# make panels of Figs S7.1 and S7.2 for manuscript entitled 
# "What is demographic lability and when might we expect to see it?"

rm(list=ls(all=TRUE))
graphics.off()

# Fig S7.1

# slow life history; adult survival varies (quadratic function) with beta environmental stochasticity

sj=.1
f=.5
sa=.95

# quadratic functions
b=.06 # slope parameter
cs=c(-.06,-.03,0,.03,.06)  # curvature: convex if c>0, concave if c<0
z=seq(-.5,.5,.01)
z2=z^2

sas=matrix(0,length(z),length(cs))
for(i in 1:length(cs)){
  sas[,i]= sa+b*z+cs[i]*z2  # quadratic function
}

# fig S7.1 A

windows()
matplot(z,sas,type="l",ylim=c(0.9,1.0),ylab="s.a",xlab="Environment, z",lty="solid",lwd=3,cex.axis=2,col='black')

# simulate log lambda.g for different values of b, beta environment

sj=.1
f=.5
sa=.95

m=0.5
v=.05 # max v = m*(1-m)
aa=m*(m*(1-m)/v - 1)
bb=(1-m)*(m*(1-m)/v - 1)

# fig S7.1 B

windows()
z=seq(-.5,.5,.01)
plot(z,dbeta(z+0.5,aa,bb),xlab="Environment, z",ylab="Density",type="l",lwd=4,cex.axis=2)


tmax=50000
zs=rbeta(tmax,aa,bb)-0.5  # z's centered on zero
zs2=zs^2
n0=matrix(c(0.5,0.5),2,1)

# sa varies, f and sj constant

# quadratic parameters
b=.06 # slope parameter
cs=seq(-.06,.06,0.002)  # curvature: convex if c>0, concave if c<0

LLg.c=array(0,length(cs))

for(i in 1:length(cs)){
  
  lambdas=array(0,tmax)
  n=n0
  
  sas = sa+b*zs+cs[i]*zs2  # quadratic function
  
  A=matrix(c(0,sj,f,0),2,2)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  LLg.c[i]=mean(log(lambdas))
  
}

# Fig S7.1 C

windows()
plot(cs,LLg.c,xlab="Curvature parameter, c",ylab="log(geometric mean)",
     type="l",lwd=3,cex.axis=2)
abline(0,0,lty='dashed',lwd=3)
abline(v=0,lty='dotted',lwd=3)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Fig S7.2

# Fig S7.2 A

# plot power functions
cs=c(.75,1,2,4,8)  # convex if c>1, concave if c<1
as=((2^cs)*sa-1)/((2^cs)-1)
z=seq(0,1,.01)

sas=matrix(0,length(z),length(cs))
for(i in 1:length(cs)){
  sas[,i]= as[i]+(1-as[i])*(z^cs[i])  # power function
}
windows()
matplot(z,sas,type="l",ylab="s.a",xlab="Environment, z",lty="solid",lwd=3,cex.axis=2,col='black')


# simulate log lambda.g for different values of c, Beta environment

sj=.1
f=.5
sa=.95

tmax=50000

m=0.5
v=.05 # max v = m*(1-m)
aa=m*(m*(1-m)/v - 1)
bb=(1-m)*(m*(1-m)/v - 1)

# Fig S7.2 B - same as Fig. S7.1 B

# plot beta pdf
windows()
plot(z,dbeta(z,aa,bb),xlab="Environment, z",ylab="Density",type="l",lwd=3,cex.axis=2)

zs=rbeta(tmax,aa,bb)
n0=matrix(c(0.5,0.5),2,1)

# sa varies, f and sj constant

cs=seq(.8,6,.2)  # convex if c>1, concave if c<1
as=((2^cs)*sa-1)/((2^cs)-1)

LLg.c=array(0,length(cs))

for(i in 1:length(cs)){
  
  lambdas=array(0,tmax)
  n=n0
  sas = as[i]+(1-as[i])*(zs^cs[i])  # power function
  A=matrix(c(0,sj,f,0),2,2)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  LLg.c[i]=mean(log(lambdas))
  
}

# Fig S7.2 C

windows()
plot(cs,LLg.c,xlab="Power function exponent, c",ylab="log(geometric mean)",
     type="l",lwd=3,cex.axis=2)
abline(0,0,lty='dashed',lwd=3)
abline(v=1,lty='dotted',lwd=3)
