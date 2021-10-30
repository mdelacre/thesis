# For each file in CIvs.test, 
# We compared the proportion of pvalues below 5% with the proportion of CIs around point estimates that include 0

Mainfolder="D:/Documents/CIvs.test/"
subfolder=list.files(Mainfolder)
Folder=paste0(Mainfolder,subfolder)

for (i in seq_len(length(Folder))){ 
  
  # set up empty container for all estimated parameters
  est <- c("_cohen","_glass1","_glass2","'_cohen","_shieh")
  parameters <- c("p","d","g") # p = proportion of pvalues below 5%
  # d = proportion of C.I around biased estimator excluding 0
  # g = proportion of C.I around UNbiased estimator excluding 0
  col.CI <- do.call(paste0, c(expand.grid(parameters,est)))
  
  finalres <-matrix(0,length(list.files(Folder[i])),length(c("n1","n2","m1","m2","sd1","sd2",col.CI)))
  colnames(finalres) <- c("n1","n2","m1","m2","sd1","sd2",col.CI)
  head(finalres)
  for (j in seq_len(length(list.files(Folder[i])))){
    
    filepath = paste0(Folder[i],"/",list.files(Folder[i])[j])
    file=readRDS(filepath)
    
    # Extracting population parameters values and sample sizes from file names
    param <- str_extract_all(list.files(Folder[i])[j], "[[:digit:]]+\\.*[[:digit:]]*")
    n1 <- as.numeric(param[[1]][5])
    n2 <- as.numeric(param[[1]][6])
    m1 <- as.numeric(param[[1]][7])
    m2 <- as.numeric(param[[1]][8])
    sd1 <- as.numeric(param[[1]][9])
    sd2 <- as.numeric(param[[1]][10])
    
    finalres[j,1:6] <- c(n1,n2,m1,m2,sd1,sd2)    
    # pval below 5%
    pval_cohen <- sum(file[,1]<.05)/length(file[,1])
    finalres[j,7] <- pval_cohen
    pval_glass1 <- sum(file[,2]<.05)/length(file[,2])
    finalres[j,10] <- pval_glass1
    pval_glass2 <- sum(file[,3]<.05)/length(file[,3])
    finalres[j,13] <- pval_glass2
    pval_cohenprime <- sum(file[,4]<.05)/length(file[,4])
    finalres[j,16] <- pval_cohenprime
    pval_shieh <- sum(file[,5]<.05)/length(file[,5])
    finalres[j,19] <- pval_shieh
    
    # proportion of C.I excluding 0
    
    CIincluding0_Bcohen <- 1-sum(file[,6])/length(file[,6])
    finalres[j,8] <- CIincluding0_Bcohen
    CIincluding0_Bglass1 <- 1-sum(file[,8])/length(file[,8])
    finalres[j,11] <- CIincluding0_Bglass1
    CIincluding0_Bglass2 <- 1-sum(file[,10])/length(file[,10])
    finalres[j,14] <- CIincluding0_Bglass2
    CIincluding0_Bcohenprime <- 1-sum(file[,12])/length(file[,12])
    finalres[j,17] <- CIincluding0_Bcohenprime
    CIincluding0_Bshieh <- 1-sum(file[,14])/length(file[,14])
    finalres[j,20] <- CIincluding0_Bshieh
    
    CIincluding0_Ucohen <- 1-sum(file[,16])/length(file[,16])
    finalres[j,9] <- CIincluding0_Ucohen
    CIincluding0_Uglass1 <- 1-sum(file[,18])/length(file[,18])
    finalres[j,12] <- CIincluding0_Uglass1
    CIincluding0_Uglass2 <- 1-sum(file[,20])/length(file[,20])
    finalres[j,15] <- CIincluding0_Uglass2
    CIincluding0_Ucohenprime <- 1-sum(file[,22])/length(file[,22])
    finalres[j,18] <- CIincluding0_Ucohenprime
    CIincluding0_Ushieh <- 1-sum(file[,24])/length(file[,24])
    finalres[j,21] <- CIincluding0_Ushieh
    
  }
  
  setwd("D:/Documents/Github_projects/Effect-sizes/Supplemental Material 4/.TXT files/") 
  shape <- str_extract_all(Folder[i], "[[:digit:]]+\\.*[[:digit:]]*")
  if(shape[[1]][2]==2.08){
    G1 <- -2.08  
  } else {G1 <- shape[[1]][2]}
  G2 = shape[[1]][4]
  
  fname=paste0("equiv_t_CI_G1=",G1,", G2=",G2,".txt")
  write.table(finalres,fname,sep=";",dec=",")
  
}
