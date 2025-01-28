# make panels of Fig 3 for manuscript entitled 
# "What is demographic lability and when might we expect to see it?"

rm(list=ls(all=TRUE))
graphics.off()

# Plot unimodal fertility and adult survival for fast and slow life histories

f.max=2.5
ff=1/f.max
a.f=log(ff/(1-ff))
b=-0.5
z=seq(-3,3,.1)  # vector of environments

# fast life history vital rates
sj=.1
f=9.9
sa=.01

a.sa=log(sa/(1-sa))
sa.z=1/(1+exp(-a.sa-b*z^2))

f.z=f.max*f/(1+exp(-a.f-b*z^2))

sa.z.fast=sa.z
f.z.fast=f.z

# slow life history
sj=.1
f=0.5
sa=.95

A=matrix(c(0,sj,f,sa),2,2)
max(eigen(A)$val)

a.sa=log(sa/(1-sa))
sa.z=1/(1+exp(-a.sa-b*z^2))

f.z=f.max*f/(1+exp(-a.f-b*z^2))

sa.z.slow=sa.z
f.z.slow=f.z

# make Fig 3A,B
windows()
matplot(z,cbind(f.z.fast,f.z.slow),type="l",
        main="Fertility",ylab="Fertility",
        lwd=4,col=c('red','black'),cex.axis=2,lty='dashed')
windows()
matplot(z,log(cbind(sa.z.fast,sa.z.slow)),type="l",
        main="Adult Survival",ylab="log(adult survival)",
        lwd=4,col=c('red','black'),cex.axis=2,lty='solid')


# compute intercept for fertility using the approach of Le Coeur et al.
f.max=2.5
ff=1/f.max
a.f=log(ff/(1-ff))

# simulate log lambda G under various scenarios, but unimodal vr's

b=-0.5 # curvature parameter

tmax=50000
n0=matrix(c(0.5,0.5),2,1)  # starting population vector
sigmas=c(.5,1)   # values of sd(z) to use
logLg=logLa=array(0,length(sigmas))

# *************************************************************

# FAST LIFE HISTORY: fertility varies, survival rates constant
sj=.1
f=9.9
sa=.01
f.max=2.5*f  # following Le Coeur et al.

for(s in 1:length(sigmas)){
  
  zs=rnorm(tmax,0,sigmas[s])  # sample environments with mean 0
  lambdas=array(0,tmax)
  n=n0  
  A=matrix(c(0,sj,0,sa),2,2)  # constant portion of projection matrix
  fs=f.max/(1+exp(-a.f-b*zs^2)) # fertility unimodal
  for(t in 1:tmax){
    A[1,2]=fs[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  logLa[s]=log(mean(lambdas))
  logLg[s]=mean(log(lambdas))
  
}
LLa.f.fast=logLa
LLg.f.fast=logLg

# SLOW LIFE HISTORY: fert varies, s's constant
sj=.1
f=0.5
sa=.95
f.max=2.5*f
# a.f=log(f/(1-f)) as above

for(s in 1:length(sigmas)){
  
  zs=rnorm(tmax,0,sigmas[s])
  lambdas=array(0,tmax)
  n=n0
  A=matrix(c(0,sj,0,sa),2,2)
  for(t in 1:tmax){
    A[1,2]=f.max/(1+exp(-a.f-b*zs[t]^2)) # fertility unimodal
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  logLa[s]=log(mean(lambdas))
  logLg[s]=mean(log(lambdas))
}
LLa.f.slow=logLa
LLg.f.slow=logLg

# FAST LIFE HISTORY: sa varies, fert and sj constant
sj=.1
f=9.9
sa=.01
a.sa=log(sa/(1-sa))

for(s in 1:length(sigmas)){
  
  zs=rnorm(tmax,0,sigmas[s])
  lambdas=array(0,tmax)
  n=n0  
  A=matrix(c(0,sj,f,0),2,2)
  for(t in 1:tmax){
    A[2,2]=1/(1+exp(-a.sa-b*zs[t]^2)) # sa unimodal
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  logLa[s]=log(mean(lambdas))
  logLg[s]=mean(log(lambdas))
}
LLa.sa.fast=logLa
LLg.sa.fast=logLg

# SLOW LIFE HISTORY: sa varies, fert and sj constant
sj=.1
f=0.5
sa=.95
a.sa=log(sa/(1-sa))

for(s in 1:length(sigmas)){
  
  zs=rnorm(tmax,0,sigmas[s])
  lambdas=array(0,tmax)
  n=n0
  A=matrix(c(0,sj,f,0),2,2)
  for(t in 1:tmax){
    A[2,2]=1/(1+exp(-a.sa-b*zs[t]^2)) # sa unimodal
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  logLa[s]=log(mean(lambdas))
  logLg[s]=mean(log(lambdas))
}
LLa.sa.slow=logLa
LLg.sa.slow=logLg

# make Fig 3C

# plot logLg's vs sigma

heights=rbind(LLg.f.fast,LLg.f.slow,LLg.sa.fast,LLg.sa.slow)
windows()
barplot(heights,beside=T,density=c(10,10,-1,-1),
        cex.axis=2,col=c('red','black','red','black'),border=T,
        legend.text=c('Fecundity/Fast','Fecundity/Slow','Adult Survival/Fast','Adult Survival/Slow'),
        args.legend = list(x = "bottomleft", cex=1.3))