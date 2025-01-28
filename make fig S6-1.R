# Make panels of Fig S6.1 for manuscript entitled 
# "What is demographic lability and when might we expect to see it?"

rm(list=ls(all=TRUE))
graphics.off()

# *************************************************************
# EXPLORE EFFECT OF INCREASING ENVTAL VARIATION CENTERED ON DIFFERENT 
# PORTIONS OF THE UNIMODAL (quadratic logistic) SURVIVAL FUNCTION FOR 
# THE SLOW LIFE HISTORY: BUT ADJUST OTHER 
# VITAL RATES TO GET LAM.1=1 IN A CONSTANT ENVIRONMENT
# Here we e

# mean vital rates for the slow life history at mean(z)=0
sj=.1
sa=.95
f=0.5
a.sa=log(sa/(1-sa)) # adult survival intercept at mean(z)=0

b=-0.5  # curvature parameter for quadratic logistic

Mz=fs=matrix(0,1,3) # array to store environmental means

# code to find values of fertility that produce determinstic lambda=1 at the mean 
# environment as the mean environment changes but the survival function is fixed 

Mz[1]=-3
(sa.low=1/(1+exp(-a.sa-b*(Mz[1])^2))) 
fs[1]=8.257; A=matrix(c(0,sj,fs[1],sa.low),2,2); (lam1=eigen(A)$val[1])

Mz[2]=-1.5
(sa.low=1/(1+exp(-a.sa-b*(Mz[2])^2))) # sa unimodal
fs[2]=1.395; A=matrix(c(0,sj,fs[2],sa.low),2,2); (lam1=eigen(A)$val[1])

Mz[3]=0
fs[3]=f # f already yields lambda=1 when mean(z)=0

# Fig S6.1A

z=seq(-5,5,.1)
sa.z=1/(1+exp(-a.sa-b*z^2))

windows()
matplot(z,cbind(sa.z),type="l",
        #ylab="Adult survival",xlab='Environment, z',
        ylab="",xlab="",
        lwd=4,col=c('black'),cex.axis=2,lty='solid')

# Fig S6.1B

z=seq(-5,5,.0005)
zd1a=dnorm(z,mean=Mz[1],sd=.1)
zd1b=dnorm(z,mean=Mz[1],sd=.5)
zd2a=dnorm(z,mean=Mz[2],sd=.1)
zd2b=dnorm(z,mean=Mz[2],sd=.5)
zd3a=dnorm(z,mean=Mz[3],sd=.1)
zd3b=dnorm(z,mean=Mz[3],sd=.5)

cbbPalette = c("#E69F00","#E69F00","#56B4E9","#56B4E9","#009E73","#009E73")
windows()
matplot(z,cbind(zd1a,zd1b,zd2a,zd2b,zd3a,zd3b),type='l',
        #xlab='Environment, z',ylab='Density',
        xlab="",ylab="",
        lwd=4,col=cbbPalette,cex.axis=2,lty='solid')


# Get stochastic lambda by simulation for each choice of mean(z)

tmax=5E5
n0=matrix(0.5,2,1)  # starting population vector
sds=seq(0,.5,0.02)  # sd(z) values to use
lam.s=sa.mean=matrix(0,length(sds),length(Mz))

for(m in 1:length(Mz)){
  
  A=matrix(c(0,sj,fs[m],0),2,2); # fill matrix with constant vr's
  
  for(i in 1:length(sds)){
    
    zs=rnorm(tmax,mean=Mz[m],sd=sds[i])  # generate random envts. with current mean and sd
    
    sas=1/(1+exp(-a.sa-b*(zs)^2)) # sa unimodal - get all adult survival rates
    sa.mean[i,m]=mean(sas)  # save arithmetic mean adult survival
    
    n=n0
    lambdas=matrix(0,1,tmax)
    for(t in 1:tmax){
      A[2,2]=sas[t]   # put variable adult survival into matrix
      n=A %*% n
      lambdas[t]=sum(n)
      n=n/sum(n)
    }
    
    lam.s[i,m]=exp(mean(log(lambdas)))
  }
}

# Fig S6.1C

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73")
windows()
matplot(sds,sa.mean,type='l',
        #xlab='SD(environment)',ylab='Arithmetic mean adult survival',
        xlab='',ylab='',
        lwd=4,col=cbbPalette,cex.axis=2,lty='solid')

# Fig S6.1D

windows()
matplot(sds,lam.s,type='l',
        #xlab='SD(environment)',ylab='Stochastic lambda',
        xlab='',ylab='',
        lwd=4,col=cbbPalette,cex.axis=2,lty='solid')

# PUSH THE EXAMPLE ABOVE WITH STRONGEST ADAPTIVE LABILITY WITH LOW SD(Z)
# TO HIGHER SD(Z)

tmax=5E5
n0=matrix(0.5,2,1)
sds=seq(0,10,0.1)  # wider range of SD(z)
lam.s=sa.mean=matrix(0,length(sds))

m=1 # this is the case where M(z)=-3 and f=8.257

A=matrix(c(0,sj,fs[m],0),2,2); # fill matrix with constant vr's

for(i in 1:length(sds)){
  
  zs=rnorm(tmax,mean=Mz[m],sd=sds[i])  # generate random envts. with current mean and sd
  
  sas=1/(1+exp(-a.sa-b*(zs)^2)) # sa unimodal - get all adult survival rates
  sa.mean[i]=mean(sas)  # save arithmetic mean adult survival
  
  n=n0
  lambdas=matrix(0,1,tmax)
  for(t in 1:tmax){
    A[2,2]=sas[t]   # but variable adult survival into matrix
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  
  lam.s[i]=exp(mean(log(lambdas)))
}

# Fig S6.1E

windows()
plot(sds,lam.s,type='l',
     #xlab='SD(environment)',ylab='Stochastic lambda',
     xlab='',ylab='',
    lwd=4,col="#E69F00",cex.axis=2,lty='solid')
lines(c(0,10),c(1,1),lty='dashed',lwd=4,col='black')