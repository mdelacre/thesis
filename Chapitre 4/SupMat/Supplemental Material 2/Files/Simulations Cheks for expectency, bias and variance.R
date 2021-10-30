library(stringr)

Folder="D:/ES MEASURES/G1=0,G2=0"

# 1) Comparing the expected mean (and variance) with the empirical mean (and variance) for all simulations where the normality assumptions is met

param <- c("emp","theo") # empirical (based on simulations) vs. theoretical (based on theoretical distribution)
lev <- c("mean","bias","var") 
estimator <- c("cohen","hedge","glass1","glass2","unbiasedglass1","unbiasedglass2","cohen d'","unbiased cohen d'","shieh","unbiasedshieh")

col.res <- do.call(paste, c(expand.grid(param,lev, estimator), sep = "_"))
res<-matrix(0,length(list.files(Folder)),length(col.res)+4) # +4 for mean diff, sd-ratio,n1 and sample sizes ratio 
colnames(res) <- c("m.diff","sd.ratio","n1","n.ratio",col.res)

for (j in seq_len(length(list.files(Folder)))){ 
  
  filepath = paste0(Folder,"/",list.files(Folder)[j])
  file=readRDS(filepath)
  # Extracting population parameters values from file names
  param <- str_extract_all(list.files(Folder)[j], "[[:digit:]]+\\.*[[:digit:]]*")
  n1 <- as.numeric(param[[1]][5])
  n2 <- as.numeric(param[[1]][6])
  N <- n1+n2
  m1 <- as.numeric(param[[1]][7])
  m2 <- as.numeric(param[[1]][8])
  sd1 <- as.numeric(param[[1]][9])
  sd2 <- as.numeric(param[[1]][10])
  nratio <- n1/n2
  
  ### Mean difference, sd-ratio and sample sizes ratio
  res[j,1] <- m1-m2
  res[j,2] <- sd1/sd2
  res[j,3] <- n1
  res[j,4] <- n1/n2
  
  ### 1A) Mean, bias, relative bias and variance for Cohen's ds
  
  # Mean
  res[j,5] <- mean(file[,9]) # empirical
  df=n1+n2-2 
  cohen_delta=(m1-m2)/sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/df)
  res[j,6] <- (cohen_delta*sqrt(df/2)*gamma((df-1)/2))/gamma(df/2) # theo
  
  # bias 
  res[j,7] <-  mean(file[,9]) - cohen_delta
  res[j,8] <-  cohen_delta*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2)-1)
  
  # Variance
  res[j,9] <-  var(file[,9])
  res[j,10] <- (N*df)/(n1*n2*(df-2))+cohen_delta^2*(df/(df-2)-(sqrt(df/2)*gamma((df-1)/2)/gamma(df/2))^2)
  
  ### 1B) Mean, bias and variance for Hedge's gs
  
  # Mean
  res[j,11] <- mean(file[,10]) # empirical
  df=n1+n2-2 
  res[j,12] <- cohen_delta
  
  # bias 
  res[j,13] <-  mean(file[,10]) - cohen_delta
  res[j,14] <-  0
  
  # Variance
  res[j,15] <-  var(file[,10])
  res[j,16] <- ((N*df)/(n1*n2*(df-2))+cohen_delta^2*(df/(df-2)-(sqrt(df/2)*gamma((df-1)/2)/gamma(df/2))^2))*(gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2)))^2
  
  ### 2A) Mean, bias and variance for Glass's ds using sd1 as standardizer
  
  # Mean
  res[j,17] <- mean(file[,11]) # empirical
  glass_delta1= (m1-m2)/sd1
  df=n1-1 
  res[j,18] <- (glass_delta1*sqrt(df/2)*gamma((df-1)/2))/gamma(df/2) # theo
  
  # bias
  res[j,19] <-  mean(file[,11]) - glass_delta1 
  res[j,20] <-  glass_delta1*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2)-1)
  
  # Variance
  res[j,21] <-  var(file[,11])
  res[j,22] <-  (df/(df-2))*(1/n1+sd2^2/(n2*sd1^2))+glass_delta1^2*(df/(df-2)-((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))^2)
  
  ### 3A) Mean, bias and variance for Glass's ds using sd2 as standardizer        
  
  # Mean
  res[j,23] <- mean(file[,12]) # empirical
  glass_delta2= (m1-m2)/sd2
  df=n2-1 
  res[j,24] <- (glass_delta2*sqrt(df/2)*gamma((df-1)/2))/gamma(df/2) # theo
               
  # bias
  res[j,25] <-  mean(file[,12]) - glass_delta2
  res[j,26] <-  glass_delta2*(-1+(sqrt((n2-1)/2)*gamma((n2-2)/2))/gamma((n2-1)/2))
  
  # Variance
  res[j,27] <-  var(file[,12])
  res[j,28] <-  (df/(df-2))*(1/n2+sd1^2/(n1*sd2^2))+glass_delta2^2*(df/(df-2)-((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))^2)

  ### 2B) Mean, bias and variance for unbiased Glass's gs using sd1 as standardizer
  
  # Mean
  res[j,29] <- mean(file[,13]) # empirical
  glass_delta1= (m1-m2)/sd1
  df=n1-1 
  res[j,30] <- glass_delta1
  
  # bias
  res[j,31] <-  mean(file[,13]) - glass_delta1 
  res[j,32] <-  0
  
  # Variance
  res[j,33] <-  var(file[,13])
  res[j,34] <-  ((df/(df-2))*(1/n1+sd2^2/(n2*sd1^2))+glass_delta1^2*(df/(df-2)-((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))^2))*(gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2)))^2
  
  ### 3B) Mean, bias and variance for unbiased Glass's gs using sd2 as standardizer        
  
  # Mean
  res[j,35] <- mean(file[,14]) # empirical
  glass_delta2= (m1-m2)/sd2
  df=n2-1 
  res[j,36] <- glass_delta2 
  
  # bias
  res[j,37] <-  mean(file[,14]) - glass_delta2
  res[j,38] <-  0
  
  # Variance
  res[j,39] <-  var(file[,14])
  res[j,40] <- ((df/(df-2))*(1/n2+sd1^2/(n1*sd2^2))+glass_delta2^2*(df/(df-2)-((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))^2))*(gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2)))^2  
  
  ### 4A) Mean, bias, relative bias and variance for Cohen's d's
  # Mean 
  res[j,41] <- mean(file[,17])
  cohen_deltaprime <- (m1-m2)/sqrt((sd1^2+sd2^2)/2)    
  df <- ((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
  res[j,42]   <- cohen_deltaprime*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))
  
  # Note: the expectency of (file[,17])*sqrt((file[,3]^2+file[,4]^2)/(n2*file[,3]^2+n1*file[,4]^2))
  # equals cohen_deltaprime*sqrt((sd1^2+sd2^2)/(n2*sd1^2+n1*sd2^2))*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))
  # with df = (sd1^2/n1+sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))
  # We can deduce this from the equation of the expectation of Shieh's ds, and from the link between Shieh's ds and Cohen's d's...
  
  # bias 
  res[j,43] <-  mean(file[,17]) - cohen_deltaprime
  res[j,44] <-  cohen_deltaprime*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2)-1)
  
  # Variance
  res[j,45] <-  var(file[,17])
  res[j,46] <- df/(df-2)*(2*(sd1^2/n1+sd2^2/n2)/(sd1^2+sd2^2))+cohen_deltaprime^2*(df/(df-2)-(sqrt(df/2)*gamma((df-1)/2)/gamma(df/2))^2)
  
  ### 4B) Mean, bias, relative bias and variance for unbiased Cohen's g's
  # Mean 
  res[j,47] <- mean(file[,18])
  cohen_deltaprime <- (m1-m2)/sqrt((sd1^2+sd2^2)/2)    
  df <- ((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4)
  
  res[j,48]   <- cohen_deltaprime
  
  # bias 
  res[j,49] <-  mean(file[,18]) - cohen_deltaprime
  res[j,50] <-  0
  
  # Variance
  res[j,51] <-  var(file[,18])
  k <- sqrt(sd1^2/n1+sd2^2/n2)/sqrt((sd1^2+sd2^2)/2)
  cf <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
  res[j,52] <- ((df/(df-2)*(2*(sd1^2/n1+sd2^2/n2)/(sd1^2+sd2^2))+cohen_deltaprime^2*(df/(df-2)-(sqrt(df/2)*gamma((df-1)/2)/gamma(df/2))^2)))*(gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2)))^2
  
  ### 5A) Mean, bias and variance for Shieh's ds
  
  # Mean
  res[j,53] <- mean(file[,15]) # empirical
  q1 <- n1/(n1+n2)
  q2 <- n2/(n1+n2)
  N <- n1+n2
  shieh_delta <- (m1-m2)/sqrt(sd1^2/q1+sd2^2/q2)    
  df <- (sd1^2/n1+sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))
  res[j,54]   <- shieh_delta*(sqrt(df/2)*gamma((df-1)/2))/gamma(df/2) # theo
  
  # bias 
  res[j,55] <-  mean(file[,15]) - shieh_delta
  res[j,56] <-  shieh_delta*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2)-1) # theo
  
  # Variance
  res[j,57] <-  var(file[,15])
  res[j,58] <- df/((df-2)*N)+shieh_delta^2*(df/(df-2)-((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))^2)
  
  ### 5B) Mean, bias and variance for unbiased Shieh's gs
  # Mean
  
  res[j,59] <- mean(file[,16]) # empirical
  q1 <- n1/(n1+n2)
  q2 <- n2/(n1+n2)
  N <- n1+n2
  shieh_delta <- (m1-m2)/sqrt(sd1^2/q1+sd2^2/q2)    
  df <- (sd1^2/n1+sd2^2/n2)^2/((sd1^2/n1)^2/(n1-1)+(sd2^2/n2)^2/(n2-1))
  res[j,60]   <- shieh_delta # theo
  
  # bias 
  res[j,61] <-  mean(file[,16]) - shieh_delta
  res[j,62] <-  0
  
  # Variance
  res[j,63] <-  var(file[,16])
  res[j,64] <- (df/((df-2)*N)+shieh_delta^2*(df/(df-2)-((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))^2))*(gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2)))^2
  
}

#setwd("D:/Documents/Github_projects/Effect-sizes/Supplemental Material 3/")
#write.table(res,"res.txt",sep=";",dec=",")

homosc_eqn <- subset(res,res[,2] == 1 & res[,4]==1) # When population variances and sample sizes are equal across groups
homosc_uneqn <- subset(res,res[,2] == 1 & res[,4]!=1) # When population variances are equal across groups and sample sizes are unequal
heterosc_uneqn <- subset(res,res[,2] != 1 & res[,4]!=1) # When population variances and sample sizes are unequal across groups
heterosc_eqn <- subset(res,res[,2] != 1 & res[,4]==1) # When population variances are unequal across groups and sample sizes are equal 

# No need to present results for bias (there are identical to those for the expectency)
biasresults <- c(7,8,13,14,19,20,25,26,31,32,37,38,43,44,49,50,55,56,61,62)
homosc_eqn <- homosc_eqn[,-biasresults] 
homosc_uneqn <- homosc_uneqn[,-biasresults]
heterosc_eqn <- heterosc_eqn[,-biasresults]
heterosc_uneqn <- heterosc_uneqn[,-biasresults]

###############################################
####                                       ####
####      Summary for Biased estimators    ####
####                                       ####
###############################################

coldescr <- c(1:4)
colbiasedest <- c(5:8,13:20,29:32,37:40)

homosc_eqn_biased <- homosc_eqn[,c(coldescr,colbiasedest)]
homosc_uneqn_biased <- homosc_uneqn[,c(coldescr,colbiasedest)]
heterosc_eqn_biased <- heterosc_eqn[,c(coldescr,colbiasedest)]
heterosc_uneqn_biased <- heterosc_uneqn[,c(coldescr,colbiasedest)]

### When population variances and sample sizes are equal across groups (condition a)
exp_cohen <- abs(homosc_eqn_biased[,5]-homosc_eqn_biased[,6]) 
exp_glass1 <- abs(homosc_eqn_biased[,9]-homosc_eqn_biased[,10])
exp_glass2 <- abs(homosc_eqn_biased[,13]-homosc_eqn_biased[,14])
exp_cohenprime <- abs(homosc_eqn_biased[,17]-homosc_eqn_biased[,18])
exp_shieh <- abs(homosc_eqn_biased[,21]-homosc_eqn_biased[,22])

var_cohen <- (homosc_eqn_biased[,7]/homosc_eqn_biased[,8])
var_glass1 <- (homosc_eqn_biased[,11]/homosc_eqn_biased[,12])
var_glass2 <- (homosc_eqn_biased[,15]/homosc_eqn_biased[,16])
var_cohenprime <- (homosc_eqn_biased[,19]/homosc_eqn_biased[,20])
var_shieh <- (homosc_eqn_biased[,23]/homosc_eqn_biased[,24])

check_biased <- cbind(homosc_eqn_biased[,1:4],exp_cohen,exp_glass1,exp_glass2,exp_cohenprime,exp_shieh,var_cohen,var_glass1,var_glass2,var_cohenprime,var_shieh)

maxcol<- apply(check_biased,2,max)
mincol<- apply(check_biased,2,min)
meancol<- apply(check_biased,2,mean)
sdcol<- apply(check_biased,2,sd)

result <- round(rbind(check_biased,max=maxcol,min=mincol,mean=meancol,sd=sdcol),3)  
write.table(result,"biased_condA.txt",sep=";",dec=",")

### When population variances are equal across groups and sample sizes are unequal (condition b)

exp_cohen <- abs(homosc_uneqn_biased[,5]-homosc_uneqn_biased[,6])
exp_glass1 <- abs(homosc_uneqn_biased[,9]-homosc_uneqn_biased[,10])
exp_glass2 <- abs(homosc_uneqn_biased[,13]-homosc_uneqn_biased[,14])
exp_cohenprime <- abs(homosc_uneqn_biased[,17]-homosc_uneqn_biased[,18])
exp_shieh <- abs(homosc_uneqn_biased[,21]-homosc_uneqn_biased[,22])

var_cohen <- (homosc_uneqn_biased[,7]/homosc_uneqn_biased[,8])
var_glass1 <- (homosc_uneqn_biased[,11]/homosc_uneqn_biased[,12])
var_glass2 <- (homosc_uneqn_biased[,15]/homosc_uneqn_biased[,16])
var_cohenprime <- (homosc_uneqn_biased[,19]/homosc_uneqn_biased[,20])
var_shieh <- (homosc_uneqn_biased[,23]/homosc_uneqn_biased[,24])

check_biased <- cbind(homosc_uneqn_biased[,1:4],exp_cohen,exp_glass1,exp_glass2,exp_cohenprime,exp_shieh,var_cohen,var_glass1,var_glass2,var_cohenprime,var_shieh)

maxcol<- apply(check_biased,2,max)
mincol<- apply(check_biased,2,min)
meancol<- apply(check_biased,2,mean)
sdcol<- apply(check_biased,2,sd)

result <- round(rbind(check_biased,max=maxcol,min=mincol,mean=meancol,sd=sdcol),3)  
write.table(result,"biased_condB.txt",sep=";",dec=",")

### When population variances are unequal across groups and sample sizes are equal (condition c)

exp_cohen <- abs(heterosc_eqn_biased[,5]-heterosc_eqn_biased[,6])
exp_glass1 <- abs(heterosc_eqn_biased[,9]-heterosc_eqn_biased[,10])
exp_glass2 <- abs(heterosc_eqn_biased[,13]-heterosc_eqn_biased[,14])
exp_cohenprime <- abs(heterosc_eqn_biased[,17]-heterosc_eqn_biased[,18])
exp_shieh <- abs(heterosc_eqn_biased[,21]-heterosc_eqn_biased[,22])

var_cohen <- (heterosc_eqn_biased[,7]/heterosc_eqn_biased[,8])
var_glass1 <- (heterosc_eqn_biased[,11]/heterosc_eqn_biased[,12])
var_glass2 <- (heterosc_eqn_biased[,15]/heterosc_eqn_biased[,16])
var_cohenprime <- (heterosc_eqn_biased[,19]/heterosc_eqn_biased[,20])
var_shieh <- (heterosc_eqn_biased[,23]/heterosc_eqn_biased[,24])

check_biased <- cbind(heterosc_eqn_biased[,1:4],exp_cohen,exp_glass1,exp_glass2,exp_cohenprime,exp_shieh,var_cohen,var_glass1,var_glass2,var_cohenprime,var_shieh)

maxcol<- apply(check_biased,2,max)
mincol<- apply(check_biased,2,min)
meancol<- apply(check_biased,2,mean)
sdcol<- apply(check_biased,2,sd)

result <- round(rbind(check_biased,max=maxcol,min=mincol,mean=meancol,sd=sdcol),3)  
write.table(result,"biased_condC.txt",sep=";",dec=",")

### When population variances and sample sizes are unequal across groups (condition d)

exp_cohen <- abs(heterosc_uneqn_biased[,5]-heterosc_uneqn_biased[,6])
exp_glass1 <- abs(heterosc_uneqn_biased[,9]-heterosc_uneqn_biased[,10])
exp_glass2 <- abs(heterosc_uneqn_biased[,13]-heterosc_uneqn_biased[,14])
exp_cohenprime <- abs(heterosc_uneqn_biased[,17]-heterosc_uneqn_biased[,18])
exp_shieh <- abs(heterosc_uneqn_biased[,21]-heterosc_uneqn_biased[,22])

var_cohen <- (heterosc_uneqn_biased[,7]/heterosc_uneqn_biased[,8])
var_glass1 <- (heterosc_uneqn_biased[,11]/heterosc_uneqn_biased[,12])
var_glass2 <- (heterosc_uneqn_biased[,15]/heterosc_uneqn_biased[,16])
var_cohenprime <- (heterosc_uneqn_biased[,19]/heterosc_uneqn_biased[,20])
var_shieh <- (heterosc_uneqn_biased[,23]/heterosc_uneqn_biased[,24])

check_biased <- cbind(heterosc_uneqn_biased[,1:4],exp_cohen,exp_glass1,exp_glass2,exp_cohenprime,exp_shieh,var_cohen,var_glass1,var_glass2,var_cohenprime,var_shieh)

maxcol<- apply(check_biased,2,max)
mincol<- apply(check_biased,2,min)
meancol<- apply(check_biased,2,mean)
sdcol<- apply(check_biased,2,sd)

result <- round(rbind(check_biased,max=maxcol,min=mincol,mean=meancol,sd=sdcol),3)  
write.table(result,"biased_condD.txt",sep=";",dec=",")

#################################################
####                                         ####
####      Summary for Unbiased estimators    ####
####                                         ####
#################################################

coldescr <- c(1:4)
colunbiasedest <- c(9:12,21:28,33:36,41:44)
homosc_eqn_unbiased <- homosc_eqn[,c(coldescr,colunbiasedest)]
homosc_uneqn_unbiased <- homosc_uneqn[,c(coldescr,colunbiasedest)]
heterosc_eqn_unbiased <- heterosc_eqn[,c(coldescr,colunbiasedest)]
heterosc_uneqn_unbiased <- heterosc_uneqn[,c(coldescr,colunbiasedest)]

### When population variances and sample sizes are equal across groups (condition a)

exp_cohen <- abs(homosc_eqn_unbiased[,5]-homosc_eqn_unbiased[,6]) 
exp_glass1 <- abs(homosc_eqn_unbiased[,9]-homosc_eqn_unbiased[,10])
exp_glass2 <- abs(homosc_eqn_unbiased[,13]-homosc_eqn_unbiased[,14])
exp_cohenprime <- abs(homosc_eqn_unbiased[,17]-homosc_eqn_unbiased[,18])
exp_shieh <- abs(homosc_eqn_unbiased[,21]-homosc_eqn_unbiased[,22])

var_cohen <- homosc_eqn_unbiased[,7]/homosc_eqn_unbiased[,8] 
var_glass1 <- homosc_eqn_unbiased[,11]/homosc_eqn_unbiased[,12]
var_glass2 <- homosc_eqn_unbiased[,15]/homosc_eqn_unbiased[,16]
var_cohenprime <- homosc_eqn_unbiased[,19]/homosc_eqn_unbiased[,20]
var_shieh <- homosc_eqn_unbiased[,23]/homosc_eqn_unbiased[,24]

check_unbiased <- cbind(homosc_eqn_unbiased[,1:4],exp_cohen,exp_glass1,exp_glass2,exp_cohenprime,exp_shieh,var_cohen,var_glass1,var_glass2,var_cohenprime,var_shieh)

maxcol<- apply(check_unbiased,2,max)
mincol<- apply(check_unbiased,2,min)
meancol<- apply(check_unbiased,2,mean)
sdcol<- apply(check_unbiased,2,sd)

result <- round(rbind(check_unbiased,max=maxcol,min=mincol,mean=meancol,sd=sdcol),3)  
write.table(result,"unbiased_condA.txt",sep=";",dec=",")

### When population variances are equal across groups and sample sizes are unequal (condition b)

exp_cohen <- abs(homosc_uneqn_unbiased[,5]-homosc_uneqn_unbiased[,6]) 
exp_glass1 <- abs(homosc_uneqn_unbiased[,9]-homosc_uneqn_unbiased[,10])
exp_glass2 <- abs(homosc_uneqn_unbiased[,13]-homosc_uneqn_unbiased[,14])
exp_cohenprime <- abs(homosc_uneqn_unbiased[,17]-homosc_uneqn_unbiased[,18])
exp_shieh <- abs(homosc_uneqn_unbiased[,21]-homosc_uneqn_unbiased[,22])

var_cohen <- homosc_uneqn_unbiased[,7]/homosc_uneqn_unbiased[,8] 
var_glass1 <- homosc_uneqn_unbiased[,11]/homosc_uneqn_unbiased[,12]
var_glass2 <- homosc_uneqn_unbiased[,15]/homosc_uneqn_unbiased[,16]
var_cohenprime <- homosc_uneqn_unbiased[,19]/homosc_uneqn_unbiased[,20]
var_shieh <- homosc_uneqn_unbiased[,23]/homosc_uneqn_unbiased[,24]

check_unbiased <- cbind(homosc_uneqn_unbiased[,1:4],exp_cohen,exp_glass1,exp_glass2,exp_cohenprime,exp_shieh,var_cohen,var_glass1,var_glass2,var_cohenprime,var_shieh)

maxcol<- apply(check_unbiased,2,max)
mincol<- apply(check_unbiased,2,min)
meancol<- apply(check_unbiased,2,mean)
sdcol<- apply(check_unbiased,2,sd)

result <- round(rbind(check_unbiased,max=maxcol,min=mincol,mean=meancol,sd=sdcol),3)  
write.table(result,"unbiased_condB.txt",sep=";",dec=",")

### When population variances are unequal across groups and sample sizes are equal (condition c)

exp_cohen <- abs(heterosc_eqn_unbiased[,5]-heterosc_eqn_unbiased[,6]) 
exp_glass1 <- abs(heterosc_eqn_unbiased[,9]-heterosc_eqn_unbiased[,10])
exp_glass2 <- abs(heterosc_eqn_unbiased[,13]-heterosc_eqn_unbiased[,14])
exp_cohenprime <- abs(heterosc_eqn_unbiased[,17]-heterosc_eqn_unbiased[,18])
exp_shieh <- abs(heterosc_eqn_unbiased[,21]-heterosc_eqn_unbiased[,22])

var_cohen <- heterosc_eqn_unbiased[,7]/heterosc_eqn_unbiased[,8] 
var_glass1 <- heterosc_eqn_unbiased[,11]/heterosc_eqn_unbiased[,12]
var_glass2 <- heterosc_eqn_unbiased[,15]/heterosc_eqn_unbiased[,16]
var_cohenprime <- heterosc_eqn_unbiased[,19]/heterosc_eqn_unbiased[,20]
var_shieh <- heterosc_eqn_unbiased[,23]/heterosc_eqn_unbiased[,24]

check_unbiased <- cbind(heterosc_eqn_unbiased[,1:4],exp_cohen,exp_glass1,exp_glass2,exp_cohenprime,exp_shieh,var_cohen,var_glass1,var_glass2,var_cohenprime,var_shieh)

maxcol<- apply(check_unbiased,2,max)
mincol<- apply(check_unbiased,2,min)
meancol<- apply(check_unbiased,2,mean)
sdcol<- apply(check_unbiased,2,sd)

result <- round(rbind(check_unbiased,max=maxcol,min=mincol,mean=meancol,sd=sdcol),3)  
write.table(result,"unbiased_condC.txt",sep=";",dec=",")

### When population variances and sample sizes are unequal across groups (condition d)

exp_cohen <- abs(heterosc_uneqn_unbiased[,5]-heterosc_uneqn_unbiased[,6]) 
exp_glass1 <- abs(heterosc_uneqn_unbiased[,9]-heterosc_uneqn_unbiased[,10])
exp_glass2 <- abs(heterosc_uneqn_unbiased[,13]-heterosc_uneqn_unbiased[,14])
exp_cohenprime <- abs(heterosc_uneqn_unbiased[,17]-heterosc_uneqn_unbiased[,18])
exp_shieh <- abs(heterosc_uneqn_unbiased[,21]-heterosc_uneqn_unbiased[,22])

var_cohen <- heterosc_uneqn_unbiased[,7]/heterosc_uneqn_unbiased[,8] 
var_glass1 <- heterosc_uneqn_unbiased[,11]/heterosc_uneqn_unbiased[,12]
var_glass2 <- heterosc_uneqn_unbiased[,15]/heterosc_uneqn_unbiased[,16]
var_cohenprime <- heterosc_uneqn_unbiased[,19]/heterosc_uneqn_unbiased[,20]
var_shieh <- heterosc_uneqn_unbiased[,23]/heterosc_uneqn_unbiased[,24]

check_unbiased <- cbind(heterosc_uneqn_unbiased[,1:4],exp_cohen,exp_glass1,exp_glass2,exp_cohenprime,exp_shieh,var_cohen,var_glass1,var_glass2,var_cohenprime,var_shieh)

maxcol<- apply(check_unbiased,2,max)
mincol<- apply(check_unbiased,2,min)
meancol<- apply(check_unbiased,2,mean)
sdcol<- apply(check_unbiased,2,sd)

result <- round(rbind(check_unbiased,max=maxcol,min=mincol,mean=meancol,sd=sdcol),3)  
write.table(result,"unbiased_condD.txt",sep=";",dec=",")
