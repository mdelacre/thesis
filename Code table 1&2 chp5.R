# Table 1: proportion d'itérations conduisant à NRHO
# Sans contrainte de puissance

semin <- function(n){

meandiff <- seq(.1,.5,.1)
NRH0 <- NULL


for (j in seq_len(length(meandiff))){

  mu1 <- 0
  mu2 <- meandiff[j]
  pval <- NULL
  
  for (i in seq_len(100000)){
    ech1 <- rnorm(n,mu1,1)
    ech2 <- rnorm(n,mu2,1)
    res <- t.test(ech1,ech2,var.equal=T)
    pval <- c(pval,res$p.value) # p-value
  }

  NRH0 <- c(NRH0,sum(pval>.05)/length(pval))
  
    }

 return(round(NRH0,3))

}

A<-semin(n=100)
B<-semin(n=200)
C<-semin(n=300)
D<-semin(n=400)
E<-semin(n=500)
F<-semin(n=600)
G<-semin(n=700)

setwd("C:/Users/Admin/Documents/Github projects/thesis/Chapitre 5/Illustration")
write.table(rbind(A,B,C,D,E,F,G),"table1.txt",sep=";",dec=",")

# Table 2: proportion d'itérations conduisant à NRHO
# Avec contrainte de puissance: puissance d'au moins 80% pour détecter
# Une différence de .1

semin2 <- function(n){
  
  meandiff <- seq(.1,.5,.1)
  NRH0 <- NULL
  for (j in seq_len(length(meandiff))){
    
    mu1 <- 0
    mu2 <- meandiff[j]
    equiv <- NULL
    
    for (i in seq_len(100000)){
      ech1 <- rnorm(n,mu1,1)
      ech2 <- rnorm(n,mu2,1)
      res <- t.test(ech1,ech2,var.equal=T)
      pval <- res$p.value # p-value
      S <- sqrt(((n-1)*sd(ech1)^2+(n-1)*sd(ech2)^2)/(2*n-2))

      # Estimation de la puissance      
      required_power <- .80
      alpha <- .05
      theta <- .3 # SESOI (exprimé en différence de moyennes)
      res <- pwr.t.test(n=n,sig.level=alpha,
                        power=required_power,
                        type="two.sample",
                        alternative="two.sided")
      d_min <- res$d # taille d'effet standardisée qui permet d'atteindre la
                     # puissance requise
      S_max <- theta/d_min # d_min <- theta/S_max <--> S_max = theta/d_min
      equivalence <- (pval>.05 & S <= S_max)
      equiv <- c(equiv,equivalence) 
    }

    NRH0 <- c(NRH0,sum(equiv==TRUE)/length(pval))
    
  }
  
  return(NRH0)
  
}

A<-semin2(n=100)
B<-semin2(n=200)
C<-semin2(n=300)
D<-semin2(n=400)
E<-semin2(n=500)
F<-semin2(n=600)
G<-semin2(n=700)

write.table(rbind(A,B,C,D,E,F,G),"table2.txt",sep=";",dec=",")

