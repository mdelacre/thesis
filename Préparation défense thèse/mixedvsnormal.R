n=100
p1=.9
mean1=0
sd1=1

p2=.1
mean2=0
sd2=10

A=rmixnorm(n, mean1 = mean1, sd1 = sd1, mean2 = mean2, sd2 = sd2, p = p1)
B=rnorm(n)

#par(mfrow=c(1,2))
#plot(density(A))
#lines(density(B),col="red")

par(mfrow=c(1,2))
qqnorm(A,main="mixed normal")
lines(seq(-10,10,.01),seq(-10,10,.01),lty=2,col="grey")
qqnorm(B,main="normal")
lines(seq(-10,10,.01),seq(-10,10,.01),lty=2,col="grey")
