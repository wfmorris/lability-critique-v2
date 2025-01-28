# Make panels of Fig. S1.2 for manuscript entitled 
# "What is demographic lability and when might we expect to see it?"

# Fig plots the combinations of curvature of a vital rate function and environmental
# driver variance (equivalently, skewness) at which the stochastic growth rate
# exceeds the growth rate that would result if the vital rate were fixed at its
# value at the mean or the median of the driver. Specifically, here we consider
# the slow life history with all vital rates fixed except for adult survival, 
# which follows a power function of driver z, with z following a right-skewed Beta
# distribution with mean 0.3 - thus greater Beta variance means greater right skew.
# Here we set the parameters of the adult survival function so that its value at
# mean(z) yields a deterministic population growth rate of 1 when the driver is 
# fixed at z=mean(z)

rm(list=ls(all=TRUE))
graphics.off()

# beta mean and variance of driver z
mu=0.3  # if mu<0.5, z is right-skewed
V.max = mu*(1-mu) # constraint on V for the beta

nvars=20  # number of values of the Beta variance to use
Vs=seq(0.01*V.max,0.2*V.max,0.19*V.max/nvars) # restrict variance to range where z distribution is unimodal
#Vs=0.01*V.max # use only minimum environmental variance
nz=length(Vs) # number of distributions of z

# use desired mean mu and variance V to get shape parameters for the beta 
as=mu*(mu*(1-mu)/Vs - 1)
bs=(1-mu)*(mu*(1-mu)/Vs - 1)

dz=0.002
z=seq(0,1,dz) # driver levels

# Fig. S1.2B

# plot the Beta distributions of z
beta.pdf=matrix(0,length(z),nz)
for(col in 1:nz) beta.pdf[,col]=dbeta(z,as[col],bs[col])
windows()
matplot(z,beta.pdf,type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.5,lty='solid')

Ms=qbeta(0.5,as,bs)  # compute medians

# vital rates for slow life history
sj=0.1
f=0.5
sa.mu=0.95 # adult survival at the mean environment

# vital rate function for sa(z) is the power function sa(z)=C+(A-C)*z^B
A=1  # sa(z)=A when z=1
Bs=seq(0.4,1.1,0.1)  # power function exponent; sa(z) concave if B<1, convex if B>1, linear if B=1
nf=length(Bs)  # number of vital rate functions
C.mu=(sa.mu-A*mu^Bs)/(1-mu^Bs) # intercept of sa(z) such that sa(mu)=sa.mu, given A and B

# Fig. S1.2A

# plot adult survival functions
sa.z=matrix(0,length(z),nf)
for(col in 1:nf) sa.z[,col]=C.mu[col]+(A-C.mu[col])*z^Bs[col]
windows()
matplot(z,sa.z,type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.5,lty='solid',xaxt='n')

P=matrix(c(0,sj,f,0),2,2)  # matrix with constant vital rates

# matrices to store results; rows = sa functions, columns = z distributions
sa.M=lambda.M=lambda.s=matrix(0,length(C.mu),length(Ms))

# compute adult survival at each MEDIAN z
for(row in 1:nf) sa.M[row,]=C.mu[row]+(A-C.mu[row])*Ms^Bs[row]

# compute deterministic lambda when z is fixed at its MEDIAN, for all sa(z) functions
P.now=P
for(row in 1:nf){
  for(col in 1:nz){
    P.now[2,2]=sa.M[row,col]
    lambda.M[row,col]=eigen(P.now)$val[1]
  }
}

# compute stochastic lambda when z varies, for all z distributions and sa(z) functions
tmax=5E4 # years to compute stochastic lambda
buffer=1E3  # number of initial years to discard before convergence to stationary distribution is reached

nreps=200 # number of times to compute stochastic lambda for each distribution of z
Bcrit=matrix(0,nreps,nz) # array to store value of power function exponent where lambda.s first exceeds lambda.M
flag=array(0,nz) # array to store whether lambda.s exceeds lambda.M even at lowest B explored

for(rep in 1:nreps){
  P.now=P
  for(col in 1:nz){ # loop over distributions of driver z
    zs=rbeta(tmax,as[col],bs[col]) # generate tmax random values of z; use with each sa function
    for(row in 1:nf){  # loop over versions of the adult survival function
      lams=array(0,tmax) # array to store annual population growth rates
      n0=matrix(c(0.5,0.5),2,1)  # starting population vector - uniform stage structure
      sa.t=C.mu[row]+(A-C.mu[row])*zs^Bs[row] # adult survival probabilities for each year
      for(t in 1:tmax){
        P.now[2,2]=sa.t[t]
        n1=P.now %*% n0
        lams[t]=sum(n1)  
        n0=n1/sum(n1) # renormalize population vector each year
      }
      lambda.s[row,col]=exp(mean(log(lams[buffer:tmax]))) # compute geometric mean growth rate
    }
  }
  delta.lambda = lambda.s-lambda.M # difference between stochastic lambda and lambda at median z
  for(col in 1:nz){
    index=max(which(delta.lambda[,col]<0))  
    if(index == - Inf) {
      flag[col]=1
      index = 1  # if delta.lambda doesn't change sign, save the lowest B used
    }
    Bcrit[rep,col]=Bs[index] # find lowest B where stochastic lambda exceeds lambda.M
  }
}
meanBcrit=colMeans(Bcrit)

# Fig. S1.2C

windows()
plot(Vs,meanBcrit,type='b',pch=NA,
     xlim=c(Vs[1],Vs[nz]),ylim=c(0.4,1.2),xlab="",ylab="",
     lwd=3,cex.axis=2.5,lty='solid')
for(s in 2:length(meanBcrit)) lines(c(Vs[s-1],Vs[s]),c(meanBcrit[s-1],meanBcrit[s]),lty='solid',lwd=3)
lines(c(Vs[1],Vs[nz]),c(1,1),lty='dashed',lwd=3)