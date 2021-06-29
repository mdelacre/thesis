semin <- function(n,nSim=100000){

meandiff <- seq(0,.5,.1)
NRH0_1 <- NULL
NRH0_2 <- NULL

for (j in seq_len(length(meandiff))){

  mu1 <- 0
  mu2 <- meandiff[j]
  equiv1 <- NULL
  equiv2 <- NULL
  
  for (i in seq_len(nSim)){
    ech1 <- rnorm(n,mu1,1)
    ech2 <- rnorm(n,mu2,1)
    res <- t.test(ech1,ech2,var.equal=T)
    
    # Sans contrainte
    equivalence1 <- (res$p.value>.05)
    equiv1 <- c(equiv1,equivalence1)    
    
    # Avec contrainte
    required_power <- .80
    alpha <- .05
    theta <- .3 # SESOI (exprimé en différence de moyennes)
    pow <- pwr.t.test(n=n,sig.level=alpha,
                      power=required_power,
                      type="two.sample",
                      alternative="two.sided")
    d_min <- pow$d # taille d'effet standardisée qui permet d'atteindre
                   # la puissance requise
    S_max <- theta/d_min # d_min <- theta/S_max <--> S_max = theta/d_min
    S <- sqrt(((n-1)*sd(ech1)^2+(n-1)*sd(ech2)^2)/(2*n-2))
    
    equivalence2 <- (res$p.value>.05 & S <= S_max) # Si S < S_max, alors on a assez de puissance
    equiv2 <- c(equiv2,equivalence2)    
  }

  NRH0_1 <- c(NRH0_1,sum(equiv1==TRUE)/nSim)
  NRH0_2 <- c(NRH0_2,sum(equiv2==TRUE)/nSim)
  
    }

 return(rbind(round(NRH0_1,3),round(NRH0_2,3)))

}

A<-semin(n=100)
B<-semin(n=200)
C<-semin(n=300)
D<-semin(n=400)
E<-semin(n=500)
F<-semin(n=600)
G<-semin(n=700)

setwd("C:/Users/Admin/Documents/Github projects/thesis/Chapitre 5/Illustration")
# Table 1: proportion d'itérations conduisant à NRHO
# sans contrainte
write.table(rbind(A[1,],B[1,],C[1,],D[1,],E[1,],F[1,],G[1,]),"table1.txt",sep=";",dec=",")

# Table 2: proportion d'itérations conduisant à NRHO
# avec contrainte
write.table(rbind(A[2,],B[2,],C[2,],D[2,],E[2,],F[2,],G[2,]),"table2.txt",sep=";",dec=",")


