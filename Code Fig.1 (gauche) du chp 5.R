################################################
####        APPROCHE DE LA PUISSANCE        ####
################################################

# paramètres de départ
alpha <- .05
required_power <- .80
n <- 50
df <- 2*n-2
theta <- 20   

# Différences de moyennes testées 
meandiff <- seq(-theta,theta,.0001) # différence de moyennes

# t critique (en deça duquel on va NRH0)
t_crit <- qt(1-alpha/2,df) # on conclura au NRH0 si t < 2.228 
# t = meandiff/SE <--> SE = meandiff/t 
# avec SE = S*sqrt(2/n) et S = écart-type poolé
# En fonction de meandiff, le SE au delà duquel on conclura
# au NRH0 sera calculé comme suit:
crit_SE <- abs(meandiff)/t_crit#  quand n1=n2=n

# Attention, on interprétera un NRH0 comme "preuve" d'équivalence QUE SI
# On atteind une puissance minimale donnée (cf. required_power; e.g.: 80%)

# Sachant cela, quel sera le SE maximum autorisé?

library(pwr)
res <- pwr.t.test(n=n,sig.level=alpha,
           power=required_power,
           type="two.sample",
           alternative="two.sided")
d_min <- res$d # taille d'effet standardisée qui permet d'atteindre la
               # puissance requise

S_max <- theta/d_min # d_min <- theta/S_max <--> S_max = theta/d_min
SE_max <- S_max*sqrt(2/n) # SE <- S*sqrt(2/N)

################################################
####                  TOST                  ####
################################################

# Différences de moyennes testées 
theta1 <- -theta
theta2 <- theta

# Rappel: équivalence si l'IC à 1-2alpha est entièrement inclus entre theta1 et theta2.
#t critique 
t_crit_TOST <- qt(1-alpha,df) # 1-alpha au lieu de 1-alpha/2, car on calcule l'IC à 1-2alpha%

#Rappel: IC autour de meandiff = meandiff+-t_crit*SE
# Si meandiff = 0 (le meilleur des cas), t_crit*SE ne peut pas dépasser theta2 (=20)
# c'est-à-dire theta2-abs(meandiff)
# t_crit*SE < theta2-abs(meandiff) <--> SE < [theta2-abs(meandiff)]/t_crit 
# = theta2/t_crit quand meandiff=0
SE_max_TOST <- theta2/t_crit_TOST # 11.03472

# Si -theta2 < meandiff < theta2, t_crit*SE ne peut pas dépasser theta2-abs(meandiff)
# ex.: si meandiff = 1, et que theta2 = 20, le demi intervalle de confiance ne 
# peut pas valoir plus de 19 sinon la borne sup sortira de la zone d'équivalence.
# idem si meandiff = -1.

# Pour conclure à l'équivalence
# t_crit_TOST*SE < theta2-abs(meandiff) <--> abs(meandiff) < theta2-t_crit_TOST*SE
# SE < [theta2-abs(meandiff)]/t_crit_TOST

crit_SE_TOST <- (theta2-abs(meandiff))/t_crit_TOST

# Représentation graphique:


setwd("C:/Users/Admin/OneDrive/Documents/Github projects/thesis/Chapitre 5/Illustration")
png("Fig1.png",width=1500,height=1500, units = "px", res = 300)
ylimsup <- round(SE_max*3/100,1)*100 # pour que la lim supérieure soit un multiple de 10
# valant approximativement le double du SE_max
plot(0,0,
     pch=19,cex=.01,ylim=c(0,ylimsup),xlim=c(-theta,theta),col="white",
     bty="n",
     xlab=expression(paste(bar(X[1])," - ",bar(X[2]))),
     ylab=expression(paste("s ",sqrt(paste(2,"/n"))))
)
axis(2,at=0:ylimsup,labels=rep("",length(0:ylimsup)),col.ticks="grey")
axis(2,at=seq(0,ylimsup,5),labels=rep("",length(seq(0,ylimsup,5))),col.ticks="black")

# Partie qui se rapporte à l'approche de la puissance
SE_max_data<-max(crit_SE[crit_SE<=SE_max]) # dans le vecteur crit_SE, quel est celui qui  
                                           # se rapproche le plus du SE_max?
segments(0,crit_SE[meandiff==0],meandiff[crit_SE==SE_max_data],SE_max_data)
segments(meandiff[crit_SE==SE_max_data][1],SE_max_data,
         meandiff[crit_SE==SE_max_data][2],SE_max_data)
x <- c(meandiff[crit_SE==SE_max_data][1],0,meandiff[crit_SE==SE_max_data][2])
y <- c(SE_max_data,0,SE_max_data)
polygon(x, y,col="grey")

# Partie qui se rapporte au TOST
segments(theta1,0,theta2,0)
segments(theta1,0,0,SE_max_dataTOST)
segments(0,SE_max_dataTOST,theta2,0)
x <- c(theta1,0,theta2)
y <- c(0,SE_max_dataTOST,0)
polygon(x, y,density=20)

dev.off()



#setwd("C:/Users/Admin/OneDrive/Documents/Github projects/thesis/Chapitre 5/Illustration")
plot(meandiff,crit_SE_TOST,pch=19,cex=.1)
abline(h=SE_max_TOST)

segments(0,crit_SE[meandiff==0],meandiff[crit_SE==SE_max_data],SE_max_data)
segments(meandiff[crit_SE==SE_max_data][1],SE_max_data,
         meandiff[crit_SE==SE_max_data][2],SE_max_data)
x <- c(meandiff[crit_SE==SE_max_data][1],0,meandiff[crit_SE==SE_max_data][2])
y <- c(SE_max_data,0,SE_max_data)
polygon(x, y,col="grey")
dev.off()

