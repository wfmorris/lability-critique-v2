# Make panels of Fig 2 for manuscript entitled 
# "What is demographic lability and when might we expect to see it?"

# **********************************************************************
# Make fig to show how a convex VRDR might evolve

rm(list=ls(all=TRUE))
graphics.off()

z=seq(-3,3,.1)
a=1
b=0.35
R1=exp(a+b*z)

int=7
slope=1.2
offset=0
asymp=6
R2=exp(log(int)+slope*(z+offset))
R2=asymp*R2/(1+R2)

# make panel A
windows()
matplot(z,cbind(R1,R2),type='l',xlab='',ylab='',   #xlab='Environment, z',ylab='Recruitment',
        col=c('red','blue'),lwd=4,lty='solid',cex.axis=2,ylim=c(0,8),xaxt='n')
axis(side=1, labels=FALSE)
lines(c(0,0),c(0,8),lwd=3,lty='dashed')

M=0
S=1
windows()
plot(z,dnorm(z,M,S),type='l',xlab='',ylab='',   #xlab='Environment, z',ylab='Density',
     col='black',lwd=4,lty='solid',cex.axis=2,ylim=c(0,0.5))
lines(c(0,0),c(0,0.5),lwd=3,lty='dashed')

zs=rnorm(1E6,M,S)
R1s=exp(1+b*zs)
R2s=exp(log(int)+slope*(zs+offset))
R2s=asymp*R2s/(1+R2s)

xmax=6.5  # max of x axis for pdf fig

d1=density(R1s)
M1a=exp(1) # R1 at mean of driver
M1b=mean(R1s)
x1=d1$x
y1=d1$y
y1=y1[x1<xmax]
x1=x1[x1<xmax]

d2=density(R2s)
M2a=exp(log(int)+slope*offset)
M2a=asymp*M2a/(1+M2a)

M2b=mean(R2s)
x2=d2$x
xplus=seq(max(x2),xmax,.01)
x2=c(x2,xplus)
y2=d2$y
y2=c(y2,rep(0,length(xplus)))

windows()
plot(x1,y1,type='l',xlab='',ylab='',yaxt='n',col='red',lwd=4,
     xlim=c(0,xmax),ylim=c(0,0.9),cex.axis=2)
axis(side=2, labels=FALSE)
lines(x2,y2,lwd=4,col='blue')
lines(c(M1a,M1a),c(0,8),lwd=2.5,lty='dashed',col='red')
lines(c(M1b,M1b),c(0,8),lwd=2.5,lty='solid',col='red')
lines(c(M2a,M2a),c(0,8),lwd=2.5,lty='dashed',col='blue')
lines(c(M2b,M2b),c(0,8),lwd=2.5,lty='solid',col='blue')