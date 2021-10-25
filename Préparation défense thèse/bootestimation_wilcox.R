# Distribution d'échantillonnage de la moyenne, quand n = 20
mean <- NULL
G1 <- NULL
G2 <- NULL
Tstat <- NULL

library(moments)
for (i in 1:20000){
  data <- rnorm(20,mean=0,sd=1)
  mean <- c(mean,mean(data))
  G1 <- c(G1,skewness(data))
  G2 <- c(G2,kurtosis(data))
  t <- (mean(data))/(sd(data)/sqrt(length(data)))
  Tstat <- c(Tstat,t)
}

# estimation par bootstrap? Je pars d'un échantillon 
#(issu d'une population normale, de taille n = 20)
# Et je ré-échantillonne pour estimer la distribution d'éch
# de la moyenne, de G1 et de G2

Sample <- rnorm(20,mean=0,sd=1)

boot_mean <- NULL
boot_G1 <- NULL
boot_G2 <- NULL
boot_Tstat <- NULL

for (i in 1:20000){
  data <- sample(Sample,size=20,replace=T)
  boot_mean <- c(boot_mean,mean(data))
  boot_G1 <- c(boot_G1,skewness(data))
  boot_G2 <- c(boot_G2,kurtosis(data))
  boot_t <- (mean(data))/(sd(data)/sqrt(length(data)))
  boot_Tstat <- c(boot_Tstat,boot_t)
  
}

# Comparaison des distr calculées (gauche, code 1) 
# et estimées par boostrap (droite, code 2)

### pour la moyenne
par(mfrow=c(1,2))
plot(density(mean),
     main="DE de la moy", # DE = distribution d'échantillonnage
     xlab=paste0("G1 = ",round(skewness(mean),3),
                 "\nG2 = ",round(kurtosis(mean),3))
    ) 
plot(density(boot_mean),
     main="bootstrap",
     xlab=paste0("G1 = ",round(skewness(boot_mean),3),
                 "\nG2 = ",round(kurtosis(boot_mean),3))
    ) 

### pour G1
par(mfrow=c(1,2))
plot(density(G1),
     main="DE de G1", # DE = distribution d'échantillonnage
     xlab=paste0("G1 = ",round(skewness(G1),3),
                 "\nG2 = ",round(kurtosis(G1),3))
) 
plot(density(boot_G1),
     main="bootstrap",
     xlab=paste0("G1 = ",round(skewness(boot_G1),3),
                 "\nG2 = ",round(kurtosis(boot_G1),3))
) 

### pour G2
par(mfrow=c(1,2))
plot(density(G2),
     main="DE de G2", # DE = distribution d'échantillonnage
     xlab=paste0("G1 = ",round(skewness(G2),3),
                 "\nG2 = ",round(kurtosis(G2),3))
)
plot(density(boot_G2),
     main="bootstrap",
     xlab=paste0("G1 = ",round(skewness(boot_G2),3),
                 "\nG2 = ",round(kurtosis(boot_G2),3))
) 

### pour la statistique t
par(mfrow=c(1,2))
plot(density(Tstat),
     main="DE de la statistique t", # DE = distribution d'échantillonnage
     xlab=paste0("G1 = ",round(skewness(Tstat),3),
                 "\nG2 = ",round(kurtosis(Tstat),3))
) 
plot(density(boot_Tstat),
     main="bootstrap",
     xlab=paste0("G1 = ",round(skewness(boot_Tstat),3),
                 "\nG2 = ",round(kurtosis(boot_Tstat),3))
) 
