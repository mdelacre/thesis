# Means vs. trimmed means

library(multicon)  

M <- NULL
tr_M <- NULL
med <- NULL

for (i in 1:20000){
  
  sample <- rnorm(20)
  M <- c(M,mean(sample))
  tr_M <- c(tr_M,mean(sample,trim=.20))
  med <- c(med,median(sample))  
  
}


par(mfrow=c(1,1))
plot(density(M),ylim=c(0,max(density(M)$y))*2)
lines(density(tr_M),lty=2)
lines(density(med),lty=3)
legend("top",lty=c(1,2),legend=c("DE of mean","DE of trimmed mean"))
