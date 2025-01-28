# Make panels of Figs S5.1 for manuscript entitled 
# "What is demographic lability and when might we expect to see it?"

rm(list=ls(all=TRUE))
graphics.off()

# fig S5.1 A

r0=0.5
a=0
b=.4
z=seq(-12,12,.1)
r=1/(1+exp(-a-b*z))
windows()
plot(z,r,xlab="Environment, z",ylab="Adult survival",
     type="l",lwd=4,cex.axis=2)
lines(c(-5,-5),c(0,.12),lty='dashed',lwd=4,col='red')
lines(c(-12,-5),c(.12,.12),lty='dashed',lwd=4,col='red')
lines(c(6,6),c(0,.92),lty='dashed',lwd=4,col='black')
lines(c(-12,6),c(.92,.92),lty='dashed',lwd=4,col='black')

# fig S5.1 B

windows()
zf=seq(-7,-3,.1)
rf=1/(1+exp(-a-b*zf))
zs=seq(4,8,.1)
rs=1/(1+exp(-a-b*zs))
x=seq(-2,2,.1)
matplot(x,cbind(rf,rs),xlab="Environment, z",ylab="Adult survival",
        type="l",lwd=4,cex.axis=2,lty='solid',col=c('red','black'))

# fig S5.1 C

a=log(2/3)
b=.4
z=seq(-12,12,.1)
r=1/(1+exp(-a-b*z))
r1=r
r2=2*r
r3=3*r
windows()
matplot(z,cbind(r1,r3),xlab="Environment, z",ylab="Fertility",
        type="l",lwd=4,cex.axis=2,col=c('black','red'),lty='solid')
abline(v=0,lty='dashed',lwd=3)

# fig S5.1 D

z=seq(-2,2,.1)
r=1/(1+exp(-a-b*z))
windows()
matplot(z,cbind(r),xlab="Environment, z",ylab="Fertility",
        type="l",lwd=4,cex.axis=2)
abline(v=0,lty='dashed',lwd=3)


