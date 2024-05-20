starttime <- Sys.time()

suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))

# Receive the parameters from the website
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
}

# Position of count files
count_files <- args[1]
# Position of meta files
meta_files <- args[2]
# Minimum ratio of gene expressed in a cell type
probs <- as.numeric(args[5])
# Number of statistical iterations
myPermutation <- as.integer(args[6])
# P-Value
probs <- 1-probs

if(args[8]=='human'){
  LR_files <- "human_mr.txt"
  enzyme_files <-  "human_enzyme.txt"
  transporter_files <-  "human_transporter.txt"
}else{
  LR_files <- "mouse_mr.txt"
  enzyme_files <-  "mouse_enzyme.txt"
  transporter_files <-  "mouse_transporter.txt"
}

myCount <- read.table(count_files, sep = "\t", row.names=1,header = T, stringsAsFactors = F)
myMeta <- read.table(meta_files, sep = "\t", header = F, stringsAsFactors = F)
colnames(myMeta) <- c('samples', "cell_type")

# Normalized
TF_IDF <- function(data){
  total_count <- colSums(data)
  
  left_part <- do.call(cbind, lapply(1:length(total_count),
                                     function(i){
                                       nrow(data)*data[,i]/total_count[i]
                                     }))
  tmp.factor <- 1/rowMeans(left_part)
  right_part <- log2(tmp.factor+1)
  
  res <- left_part*right_part
  return(as.data.frame(res))
}

getNullData <- function(mydata, meta, probs = 0.1){
  data <- mydata
  my.ncol <- ncol(data)
  cell.name <- names(table(meta$cell_type))
  cell.num <- as.numeric(table(meta$cell_type))
  
  my.colnames <- rep(cell.name, cell.num)
  new.data <- data[,sample(x = 1:my.ncol, size = my.ncol, replace = FALSE)]
  new.data <- data.frame(new.data, check.names = FALSE)
  colnames(new.data) <- my.colnames
  
  expr_mean <- matrix(nrow = nrow(new.data), ncol = length(cell.name))
  myColnames <- c()
  res.df <- do.call(cbind ,lapply(cell.name, function(x){
    myCell <- x
    myMatrix <- new.data[,colnames(new.data)==myCell,drop=F]
    
    quantil.tmp <- as.numeric(apply(myMatrix, 1, function(x){
      quantile(x, probs = probs,names=FALSE)
    }))
    mean.tmp <- rowMeans(myMatrix)
    mean.tmp[which(quantil.tmp==0)]<-0 
    myMatrix_mean <- mean.tmp
    return(myMatrix_mean)
  }))
  colnames(res.df) <- cell.name
  return(res.df)
}

getMeanData <- function(mydata, meta, probs = 0.1){
  data <- mydata
  my.ncol <- ncol(data)
  cell.name <- names(table(meta$cell_type))
  cell.num <- as.numeric(table(meta$cell_type))
  
  my.colnames<-meta[match(colnames(data),as.character(meta$samples)),2]
  new.data <- data
  new.data <- data.frame(new.data, check.names = FALSE)
  colnames(new.data) <- my.colnames
  
  expr_mean <- matrix(nrow = nrow(new.data), ncol = length(cell.name))
  myColnames <- c()
  res.df <- do.call(cbind ,lapply(cell.name, function(x){
    myCell <- x
    myMatrix <- new.data[,colnames(new.data)==myCell,drop=F]
    quantil.tmp <- as.numeric(apply(myMatrix, 1, function(x){
      quantile(x, probs = probs,names=FALSE)
    }))
    mean.tmp <- rowMeans(myMatrix)
    mean.tmp[which(quantil.tmp==0)]<-0 # probs 以上全为0的gene均值也为 0
    myMatrix_mean <- mean.tmp## 每个细胞类型的gene求均值
    return(myMatrix_mean)
  }))
  colnames(res.df) <- cell.name
  return(res.df)
}

LR <- read.table(LR_files, sep = "\t", header = T, quote="", stringsAsFactors = F)
enzyme <- read.table(enzyme_files, sep = "\t", header = T, quote="", stringsAsFactors = F)
transporter <- read.table(transporter_files, sep = "\t", header = T, quote="", stringsAsFactors = F)

gene.needed <- LR$Receptor_symbol %>% str_split(";", simplify = F) %>% unlist() %>% unique()
gene.needed <- unique(c(gene.needed,enzyme$GENE_NAME,transporter$GENE_NAME))
gene.detected <- intersect(rownames(myCount), gene.needed)

myCount <- myCount[match(gene.detected,rownames(myCount)),]
myCount <- myCount[(rowSums(myCount)!=0),(colSums(myCount)!=0)]

myCount.tmp <- TF_IDF(myCount)
colnames(myCount.tmp) <- colnames(myCount)
rownames(myCount.tmp) <- rownames(myCount)

myCount <- myCount.tmp
myMeta <- myMeta[match(intersect(colnames(myCount),as.character(myMeta$samples)),as.character(myMeta$samples)),]
myCount <- myCount[,match(intersect(colnames(myCount),as.character(myMeta$samples)),colnames(myCount))]

simuDat.mean <- getMeanData(mydata = myCount, meta = myMeta, probs = probs)
simuDat.mean <- simuDat.mean[apply(simuDat.mean, 1, sum)!=0,,drop=F]

myCount <- myCount[rownames(simuDat.mean),,drop=F]
null.list = lapply(1:myPermutation, function(i){
  # print(i)
  getNullData(mydata = myCount, meta = myMeta, probs = probs)
})

mylr <- LR %>% filter( Receptor_symbol %in% rownames(simuDat.mean)) %>% .[,c(1,2,4)]

myValue.list <- list()
myP.list <- list()
mycolnames <- c()

for(s in colnames(simuDat.mean)){
  for(r in colnames(simuDat.mean)){
    sender <- s
    receiver <- r
    mycolnames <- c(mycolnames, paste(sender, receiver, sep = '-'))
    val_Pvalue <- apply(mylr, 1, function(x){
      en<-enzyme[which(enzyme$HMDB_ID==x[1]),]
      en_s<-en[which(en$enzyme_p_s=="s"),"GENE_NAME"]
      en_p<-en[which(en$enzyme_p_s=="p"),"GENE_NAME"]
    	s_gene<-intersect(en_s,rownames(simuDat.mean))
    	p_gene<-intersect(en_p,rownames(simuDat.mean))
    	if(prod(simuDat.mean[s_gene,sender])!=0 && prod(simuDat.mean[p_gene,sender])!=0){
    	    real_ligand_sval<-prod(simuDat.mean[s_gene,sender])
    	    real_ligand_pval<-prod(simuDat.mean[p_gene,sender])
    	  }else{
    	    s_simuDat.mean<-simuDat.mean[s_gene,sender]
    	    p_simuDat.mean<-simuDat.mean[p_gene,sender]
    	    real_ligand_sval<-prod(s_simuDat.mean[s_simuDat.mean!=0])
    	    real_ligand_pval<-prod(p_simuDat.mean[p_simuDat.mean!=0])
    	  }
    	E.real.val <- (real_ligand_pval^(1/length(p_gene)) - real_ligand_sval^(1/length(s_gene)))

    	tp<-transporter[which(transporter$HMDB_ID==x[1]),"GENE_NAME"]
    	tp_gene<-intersect(tp,rownames(simuDat.mean))
    	if(prod(simuDat.mean[tp_gene,sender])!=0){
    	    real_ligand_tpval<-prod(simuDat.mean[tp_gene,sender])
    	}else{
    	    tp_simuDat.mean<-simuDat.mean[tp_gene,sender]
    	    real_ligand_tpval<-prod(tp_simuDat.mean[tp_simuDat.mean!=0])
    	}
    	T.real.val <- real_ligand_tpval^(1/length(tp_gene))

      if(length(str_split(x[3], ";", simplify = F)[[1]])>1){
        i_tmp = str_split(x[3], ";", simplify = F)[[1]]
        if(sum(i_tmp %in% gene.detected) == length(i_tmp)){
          if(sum(simuDat.mean[i_tmp,receiver]!=0)==length(i_tmp)){
            real_receptor_val = mean(simuDat.mean[i_tmp,receiver])
          }else{
            return(c(0, 1))
          }
        }else{
          return(c(0, 1))
        }
      }else{
        real_receptor_val = simuDat.mean[x[3],receiver]
      }
      
      if((T.real.val==0 && E.real.val ==0) || real_receptor_val==0){
        return(c(0, 1))
      }
      
      real.val <- (T.real.val^2 + E.real.val^2 + real_receptor_val^2)^(1/2)

      null.tmp <- as.numeric(unlist(lapply(null.list, function(d){
        
        en<-enzyme[which(enzyme$HMDB_ID==x[1]),]
        en_s<-en[which(en$enzyme_p_s=="s"),"GENE_NAME"]
        en_p<-en[which(en$enzyme_p_s=="p"),"GENE_NAME"]
	  s_gene<-intersect(en_s,rownames(d))
	  p_gene<-intersect(en_p,rownames(d))
	  if(prod(d[s_gene,sender])!=0 && prod(d[p_gene,sender])!=0){
	    real_ligand_sval<-prod(d[s_gene,sender])
	    real_ligand_pval<-prod(d[p_gene,sender])
	  }else{
	    s_simuDat.mean<-d[s_gene,sender]
	    p_simuDat.mean<-d[p_gene,sender]
	    real_ligand_sval<-prod(s_simuDat.mean[s_simuDat.mean!=0])
	    real_ligand_pval<-prod(p_simuDat.mean[p_simuDat.mean!=0])
	  }
	E.real.val <- (real_ligand_pval^(1/length(p_gene)) - real_ligand_sval^(1/length(s_gene)))

	tp<-transporter[which(transporter$HMDB_ID==x[1]),"GENE_NAME"]
	tp_gene<-intersect(tp,rownames(d))
	if(prod(d[tp_gene,sender])!=0){
	    real_ligand_tpval<-prod(d[tp_gene,sender])
	}else{
	    tp_simuDat.mean<-d[tp_gene,sender]
	    real_ligand_tpval<-prod(tp_simuDat.mean[tp_simuDat.mean!=0])
	}
	T.real.val <- real_ligand_tpval^(1/length(tp_gene))

        
        if(length(str_split(x[3], ";", simplify = F)[[1]])>1){
          i_tmp = str_split(x[3], ";", simplify = F)[[1]]
          
          if(sum(d[i_tmp,receiver]!=0)==length(i_tmp)){
            real_receptor_val = mean(d[i_tmp,receiver])
          }else{
            return(c(0))
          }
        }else{
          real_receptor_val = d[x[3],receiver]
        }
        
        if((T.real.val==0 && E.real.val ==0) || real_receptor_val==0){
          return(c(0))
        }
        real.val <- (T.real.val^2 + E.real.val^2 + real_receptor_val^2)^(1/2)
        return(real.val)
      })))

      p.val.tmp.1 <- (length(null.tmp)-sum(real.val > null.tmp))/length(null.tmp)
      if(p.val.tmp.1==0){
        p.val.tmp.1 = 1/myPermutation
      }
      return(c(real.val, p.val.tmp.1))
    })

    myValue.list <- c(myValue.list, list(val_Pvalue[1,]))
    myP.list <- c(myP.list, list(val_Pvalue[2,]))
  }
}

val.df <- do.call(cbind, myValue.list)
p.df <- do.call(cbind, myP.list)

colnames(val.df) <- mycolnames
colnames(p.df) <- mycolnames
rownames(val.df) <- paste(mylr[,2], mylr[,3], sep = " - ")
rownames(p.df) <- paste(mylr[,2], mylr[,3], sep = " - ")

logic_index <- apply(p.df, 1, function(x){
  sum(x<as.numeric(args[7]))
})>0

tmp.df <- val.df[logic_index,]
my.tmp.df <- (tmp.df-min(tmp.df))/(max(tmp.df)-min(tmp.df))

write.table(p.df[logic_index,], args[3], sep = '\t', col.names = T, row.names = T, quote = F)
write.table(my.tmp.df, args[4], sep = '\t', col.names = T, row.names = T, quote = F)

endtime <- Sys.time()
