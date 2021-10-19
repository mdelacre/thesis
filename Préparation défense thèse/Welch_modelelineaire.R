set.seed(5) # pour avoir toujours les mêmes valeurs
n1 <- 27
n2 <- 15
Gr1 <- rnorm(n1)
Gr2 <- rnorm(n2)
y <- c(Gr1,Gr2)
x <- c(rep(1,n1),rep(2,n2))

#Representation graphique des points + descriptives

plot(y~x)
means <- tapply(y, x, mean) # Moyenne par groupe
m0 <- means[1]
m1 <- means[2]

#################################################################
### Les valeurs de b0 et b1 via l'approche du modèle linéaire ###
#################################################################

# (a) avec le dummy coding
x <- c(rep(0,n1),rep(1,n2))
regr<-summary(lm(y~1+x)) # shorthand for y = 1*b0+x*b1

round(regr$coefficients[1],5)==round(m0,5) # y = b0+b1x --> quand x = 0, y = b0
                                                 # b0 = moyenne du groupe codé 0
round(regr$coefficient[2],5)==round(m1-m0,5) # b1 = pente = delta_y/delta_x = m1-m0

# (b) Avec des pr�dicteurs centr�s
x <- c(rep(-.5,n1),rep(.5,n2)) 
regr<-summary(lm(y~x))

round(regr$coefficients[1],5)==round((m0+m1)/2,5) # b0 est à mi-chemin entre les 2 moyennes, 
                                                  # et correspond donc à la moyenne non pondérée
                                                  # des deux moyennes

round(regr$coefficient[2],5)==round(m1-m0,5) # b1 = pente = delta_y/delta_x = (m1-m0)/1 = (m1-m0)

####################
### t de Student ###
####################

library(nlme)
regr_student<-summary(nlme::gls(y~1+x)) # shorthand for y = 1*b0+x*b1
ttest_student<-t.test(Gr1,Gr2,var.equal=TRUE)

# Via la regression et le menu du test t, 
#on obtient exactement les memes...

# ... valeurs de statistique (au signe près)
t_regr_student<-abs(regr_student$tTable["x","t-value"])
t_ttest_student<-abs(ttest_student$statistic)
round(t_regr_student,5)==round(t_ttest_student,5)

# ... degres de liberte
df_regr_student<-regr_student # impossible de l'extraire, donc lire la dernière valeur residual
df_ttest_student<-ttest_student$parameter
df_regr_student
df_ttest_student

# ... p-valeurs
p_regr_student<-regr_student$tTable["x","p-value"]
p_ttest_student<-ttest_student$p.value
round(p_regr_student,5)==round(p_ttest_student,5)

# ... et IC autour de b1
bornes_regr_student<-abs(confint(lm(y~1+x))["x",])
c(min(bornes_regr_student),max(bornes_regr_student))
bornes_ttest_student <-abs(c(ttest_student$conf.int[1],ttest_student$conf.int[2]))
c(min(bornes_ttest_student),max(bornes_ttest_student))

# Comprendre pourquoi avec nlme, j'ai des bornes differentes? 

####################
###  t de Welch  ###
####################

library(nlme)
library(lmerTest)

regr_welch<-nlme::gls(y~1+x,weights=nlme::varIdent(form=~1|x),method="REML") # shorthand for y = 1*b0+x*b1
ttest_welch<-t.test(Gr1,Gr2,var.equal=FALSE)

# Via la regression et le menu du test t, 
#on obtient exactement les memes...

# ... valeurs de statistique (au signe pres)
t_regr_welch<-abs(summary(regr_welch)$tTable["x","t-value"])
t_ttest_welch<-abs(ttest_welch$statistic)
round(t_regr_welch,5)==round(t_ttest_welch,5)

# ... degres de liberte

#library(emmeans)
#regr_welch.contrat <- contrast(emmeans(regr_welch,specs="x"),method="pairwise")

df_regr_welch<-regr_welch # impossible de l'extraire, donc lire la dernière valeur residual
df_ttest_welch<-ttest_welch$parameter
df_regr_welch
df_ttest_welch

# ... p-valeurs
p_regr_welch<-regr_welch$tTable["x","p-value"]
p_ttest_welch<-ttest_welch$p.value
round(p_regr_welch,5)==round(p_ttest_welch,5)

