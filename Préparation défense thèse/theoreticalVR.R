for (package in c("PearsonDS","gsl")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


##### Function to load packages

# Note: when kurtosis = 3 and skewness = 0, one has a normal distribution

setwd("C:/Users/mdelacre/Dropbox/Préparation défense thèse/VR expe")

get_write <- function(object,
                      n1=40,n2=40,
                      m1=1,m2=0){
  
  # compute mean and standard deviation
  nobs <- paste0(" n=[",n1,",",n2,"]")
  means <-  paste0(" means=[",m1,",",m2,"]")
  fname <-  paste(nobs, means, sep=",")
  fname <-  paste0(fname, ".rds")
  saveRDS(object, file = fname)
  
}

##### Function to generate dataset,realize statistical test and compute (and extract) p-value

get_simu     <- function(nSims=10000,
                         n1=50,n2=50,  
                         m1=1,m2=0){

  # set up empty container for all estimated parameters
  VR1 <-rep(0,nSims) 
  VR2 <-rep(0,nSims) 
  VR3 <-rep(0,nSims) 

  for (i in 1:nSims){

    #homoscedasticity parameters
    sd1 <- 2
    sd2 <- 2

    # normal distributions
    y1 <- rnorm(n1,m1,sd1)
    y2 <- rnorm(n2,m2,sd2)

    # moderately skewed
    kurt2 <- 95.75 # when asymetry, kurtosis > 0, always!
    skew2 <- -2.08
    y3 <- rpearson(n1,moments=c(m1,sd1^2,skewness=skew2*(n1-2)/sqrt(n1*(n1-1)),kurtosis=(kurt2*(n1-2)*(n1-3)-6*(n1-1))/(n1^2-1)+3)) # Note: the function takes the variance                                                                              # as an argument (not the sd)
    y4 <- rpearson(n2,moments=c(m2,sd2^2,skewness=skew2*(n2-2)/sqrt(n2*(n2-1)),kurtosis=(kurt2*(n2-2)*(n2-3)-6*(n2-1))/(n2^2-1)+3))

    # highly skewed
    kurt3 <- 95.75 # when asymetry, kurtosis > 0, always!
    skew3 <- 6.32
    y5 <- rpearson(n1,moments=c(m1,sd1^2,skewness=skew3*(n1-2)/sqrt(n1*(n1-1)),kurtosis=(kurt3*(n1-2)*(n1-3)-6*(n1-1))/(n1^2-1)+3)) # Note: the function takes the variance                                                                              # as an argument (not the sd)
    y6 <- rpearson(n2,moments=c(m2,sd2^2,skewness=skew3*(n2-2)/sqrt(n2*(n2-1)),kurtosis=(kurt3*(n2-2)*(n2-3)-6*(n2-1))/(n2^2-1)+3))
    
    # For the explanation about the transformation of skew1,skew2,kurt1 and kurt2,
    # See Skewness-and-kurtosis-in-my-simulations,-note.pdf
    ### Descriptives
    emp_vary1 <- var(y1)
    emp_vary2 <- var(y2)
    VR1[i] <- max(c(emp_vary1,emp_vary2))/min(c(emp_vary1,emp_vary2)) #normal case

    emp_vary3 <- var(y3)
    emp_vary4 <- var(y4)
    VR2[i] <- max(c(emp_vary3,emp_vary4))/min(c(emp_vary3,emp_vary4)) # moderately skewed

    emp_vary5 <- var(y5)
    emp_vary6 <- var(y6)
    VR3[i] <- max(c(emp_vary5,emp_vary6))/min(c(emp_vary5,emp_vary6)) # highly skewed
    
  } 
  
  # Extraction of the ES matrix 
  chem <- "C:/Users/mdelacre/Dropbox/Préparation défense thèse/VR expe"
  setwd(dir=chem) # destination file  
  get_write(cbind(VR1,VR2,VR3),n1,n2,m1,m2)

}

#####################################################################################
#############                       APPLICATIONS                        #############
#####################################################################################

n1 <- c(20,50,100)
n2 <- c(20,50,100)
m1 <- seq(0,4,1)
m2 <- 0

Simu=expand.grid(n1,n2,m1,m2)
length(Simu[,1])
colnames(Simu)<-c("n1","n2","m1","m2")

# performing simulations  
for (i in seq_len(length(Simu[,1]))){
  get_simu(nSims=10000,n1=Simu[i,1],n2=Simu[i,2],m1=Simu[i,3],m2=Simu[i,4])    
}

#####################################################################################
#############                       Calculs rapides                     #############
#####################################################################################

########## Proportion de VR > 3

N1 <- NULL
N2 <- NULL
M1 <- NULL
M2 <- NULL
PROPVR1 <- NULL 
PROPVR2 <- NULL
PROPVR3 <- NULL

chem<-"C:/Users/mdelacre/Dropbox/Préparation défense thèse/VR expe/ "

for (i in 1:length(Simu[,1])){
  A<-readRDS(paste0(chem,"n=[",Simu[i,1],",",Simu[i,2],"], means=[",Simu[i,3],",",Simu[i,4],"].rds"))  
  param <- str_extract_all(list.files(chem)[i], "[[:digit:]]+\\.*[[:digit:]]*")
  n1 <- as.numeric(param[[1]][1])
  n2 <- as.numeric(param[[1]][2])
  m1 <- as.numeric(param[[1]][3])
  m2 <- as.numeric(param[[1]][4])
  prop_VR1sup3 <- sum(A[,1]>3)/length(A[,1])
  prop_VR2sup3 <- sum(A[,2]>3)/length(A[,2])
  prop_VR3sup3 <- sum(A[,3]>3)/length(A[,3])
  
  N1=c(N1,n1)
  N2=c(N2,n2)
  M1=c(M1,m1)
  M2=c(M2,m2)
  PROPVR1 <- c(PROPVR1,prop_VR1sup3)
  PROPVR2 <- c(PROPVR2,prop_VR2sup3)
  PROPVR3 <- c(PROPVR3,prop_VR3sup3)
  
}

result <-cbind(N1,N2,M1,M2,PROPVR1,PROPVR2,PROPVR3)

########## VR moyen

N1 <- NULL
N2 <- NULL
M1 <- NULL
M2 <- NULL
MEDVR1 <- NULL 
MEDVR2 <- NULL
MEDVR3 <- NULL
MADVR1 <- NULL 
MADVR2 <- NULL
MADVR3 <- NULL

chem<-"C:/Users/mdelacre/Dropbox/Préparation défense thèse/VR expe/ "

for (i in 1:length(Simu[,1])){
  A<-readRDS(paste0(chem,"n=[",Simu[i,1],",",Simu[i,2],"], means=[",Simu[i,3],",",Simu[i,4],"].rds"))  
  param <- str_extract_all(list.files(chem)[i], "[[:digit:]]+\\.*[[:digit:]]*")
  n1 <- as.numeric(param[[1]][1])
  n2 <- as.numeric(param[[1]][2])
  m1 <- as.numeric(param[[1]][3])
  m2 <- as.numeric(param[[1]][4])
  med_VR1 <- median(A[,1])
  med_VR2 <- median(A[,2])
  med_VR3 <- median(A[,3])

  mad_VR1 <- mad(A[,1])
  mad_VR2 <- mad(A[,2])
  mad_VR3 <- mad(A[,3])
  
  N1=c(N1,n1)
  N2=c(N2,n2)
  M1=c(M1,m1)
  M2=c(M2,m2)
  MEDVR1 <- c(MEDVR1,med_VR1)
  MEDVR2 <- c(MEDVR2,med_VR2)
  MEDVR3 <- c(MEDVR3,med_VR3)
  MADVR1 <- c(MADVR1,mad_VR1)
  MADVR2 <- c(MADVR2,mad_VR2)
  MADVR3 <- c(MADVR3,mad_VR3)
  
}


result <-cbind(N1,N2,M1,M2,MEDVR1,MEDVR2,MEDVR3,MADVR1,MADVR2,MADVR3)
write.table(result,"res.txt",sep=";",dec=",")

########## VR max

N1 <- NULL
N2 <- NULL
M1 <- NULL
M2 <- NULL
MAXVR1 <- NULL 
MAXVR2 <- NULL
MAXVR3 <- NULL

chem<-"C:/Users/mdelacre/Dropbox/Préparation défense thèse/VR expe/ "

for (i in 1:length(Simu[,1])){
  A<-readRDS(paste0(chem,"n=[",Simu[i,1],",",Simu[i,2],"], means=[",Simu[i,3],",",Simu[i,4],"].rds"))  
  param <- str_extract_all(list.files(chem)[i], "[[:digit:]]+\\.*[[:digit:]]*")
  n1 <- as.numeric(param[[1]][1])
  n2 <- as.numeric(param[[1]][2])
  m1 <- as.numeric(param[[1]][3])
  m2 <- as.numeric(param[[1]][4])
  max_VR1 <- max(A[,1])
  max_VR2 <- max(A[,2])
  max_VR3 <- max(A[,3])
  

  N1=c(N1,n1)
  N2=c(N2,n2)
  M1=c(M1,m1)
  M2=c(M2,m2)
  MAXVR1 <- c(MAXVR1,max_VR1)
  MAXVR2 <- c(MAXVR2,max_VR2)
  MAXVR3 <- c(MAXVR3,max_VR3)

}

resultmax <-cbind(N1,N2,M1,M2,MAXVR1,MAXVR2,MAXVR3)
write.table(resultmax,"resmax.txt",sep=";",dec=",")