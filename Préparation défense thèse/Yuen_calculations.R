
# Série de base
A<-c(1,3,5,6,7,10,11,12,75,75)
B<-c(1,4,5,8,19,20,21,21,23,24)

### Trimming
TR=0.20

### Série utilisée pour calculer la moyenne trimmée:
### suppression des 20% des scores les plus faibles
### et des 20% des scores les plus élevés

trimA <- sort(A)[(TR*length(A)+1):(length(A)-(TR*length(A)))]
# si trimming = 20% : trimA=c(5,6,7,10,11,12)
trimB<-sort(B)[(TR*length(B)+1):(length(B)-(TR*length(B)))]
# si trimming = 20% : trimB=c(5,8,19,20,21,21)

# Moyennes trimmées
trim_m1=mean(trimA)
trim_m2=mean(trimB)

# Longueur des moyennes trimmées
h1=length(trimA)
h2=length(trimB)

### Série utilisée pour calculer la variance/e.t. winsorisée:
### On remplace les scores trimmés par les plus extrêmes
### valeurs non trimmées
winsorA<- c(rep(trimA[1],length(A)*TR),
            trimA,
            rep(trimA[length(trimA)],length(A)*TR))
# si trimming = 20% : winsorA=c(5,5,5,6,7,10,11,12,12,12)
winsorB<- c(rep(trimB[1],length(B)*TR),
            trimB,
            rep(trimB[length(trimB)],length(B)*TR))
# si trimming = 20%: winsorB<-c(5,5,5,8,19,20,21,21,21,21)

# Ecart-types winsorisées
winsor_sd1=sd(winsorA) # sqrt(sum((winsorA-mean(winsorA))^2)/(n1-1))
winsor_sd2=sd(winsorB) # sqrt(sum((winsorB-mean(winsorB))^2)/(n2-1))


#library(multicon)  
# sd1 = sqrt(winvar(A,tr=TR))
# sd2 = sqrt(winvar(B,tr=TR))

# Longueur des séries de départ (et winsorisées)
n1 = length(winsorA)
n2 = length(winsorB)

#Estimation de l'erreur standard faite par Yuen
d1 = ((n1-1)*winsor_sd1^2)/(h1*(h1-1))
d2 = ((n2-1)*winsor_sd2^2)/(h2*(h2-1))

winsor_sd1^2/n1
winsor_sd2^2/n2

t = yuen.t.test(A,B, alternative = "two.sided",tr=TR)$statistic
t = (trim_m1-trim_m2)/sqrt(d1+d2)

# Note: quand TR=0, le test de Yuen = le test de Welch (cf. démo manuscripte)
# vrai aussi pour ses df, et donc forcément pour sa p-valeur.

