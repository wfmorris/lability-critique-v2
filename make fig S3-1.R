# Make panels of Fig. S3.1 for manuscript entitled 
# "What is demographic lability and when might we expect to see it?"

# Code computes the localized selection gradient to small perturbations to a vital  
# rate function near specific values of the environmental driver. Specifically, we 
# use a Normal distribution with narrow variance and low maximum to push up a convex
# adult survival curve near specific values of the driver z, keeping adult survival 
# unchanged at values of z away from the center of the perturbation. We compute the
# localized selection gradient as the difference in the stochastic growth rate 
# between the perturbed and original survival curves divided by the maximum 
# perturbation. Here we use the slow life history with all vital rates except
# for adult survival fixed, a quadratic adult survival curve, and a shifted
# Beta-distributed environmental distribution that is symmetrical and centered 
# on z=0. By varying the center of the normally distributed perturbation, we 
# quantify localized selection on the vital rate function across the environmental
# distribution. In the second part of the script, we examine the case of antagonistic
# pleiotropy by combining a positive perturbation across the driver range with a 
# negative perturbation of equal amplitude near the upper end of the environmental
# distribution, where the (convex) adult survival function reaches its highest values.

rm(list=ls(all=TRUE))
graphics.off()

# baseline values - SLOW LIFE HISTORY - yields lambda=1 for z fixed at 0
sa=0.95
sj=0.1
f=0.5

z=seq(0,1,.01)-.5 # vector of driver values
z2=z^2

# the vital rate function; same as solid green curve in Fig. 1A
b=0.06
c=0.06  # convex if c>1, concave if c<1
sa.z = sa+b*z+c*z2  # quadratic vital rate function

# Fig. S3.1A

# plot the vital rate function - cf. Fig. 1A, solid green curve
windows()
plot(z,sa.z,
     #xlab='Environment, z',ylab='Adult survival',
     xlab='',ylab='',xaxt='n',
     type="l",ylim=c(.9,1),
        col='green',lwd=4,lty='solid',cex.axis=2)

tmax=1e6    # time steps to compute lam.s
trash=1e4   # initial time steps to discard 
n0=matrix(c(0.5,0.5),2,1)  # starting population vector
A=matrix(c(0,sj,f,0),2,2)  # constant parts of projection matrix

# generate Beta-distributed environmental states - same as distribution in Fig. 1B

# compute Beta parameters from desired mean (m) and variance (v)
m=0.5
v=.05 # max v = m*(1-m)
aa=m*(m*(1-m)/v - 1)
bb=(1-m)*(m*(1-m)/v - 1)

# Fig. S3.1B

# plot the Beta pdf
windows()
z=seq(-.5,.5,.01)
plot(z,dbeta(z+0.5,aa,bb),
     #xlab="Environment, z",ylab="Density",
     xlab='',ylab='',xaxt='n',type="l",lwd=4,cex.axis=2)

# generate tmax random enviromental states
zs=rbeta(tmax,aa,bb)
zs=zs-0.5   # make z vary from -1/2 to +1/2`
zs2=zs^2    # only need to do this once

# allocate memory
LLp=LLap=array(0,length(z)) # arrays for log lambda.s
# arrays for checking sa remains in the range (0,1)
samin.1=samax.1=samin.2=samax.2=LLp  
  
#***** Get the baseline stochastic growth rate ******

lambdas=array(0,tmax)
n=n0
sas = sa+b*zs+c*zs2   # quadratic function - make all tmax s.a's
samin.0=min(sas)   # check that sa stays in bounds
samax.0=max(sas)
for(t in 1:tmax){
  A[2,2]=sas[t]
  n=A %*% n
  lambdas[t]=sum(n)
  n=n/sum(n)
}
lambdas=lambdas[-(1:trash)]  # discard initial lambdas - pre-stationary distribution?
LL0=mean(log(lambdas)) # compute baseline stochastic lambda

# ********* Now perturb the adult survival function locally ********************

# perturbation parameters
delta=0.005 # maximum perturbation of sa; with baseline coefficients of the sa 
            # function, and delta=0.005, perturbed sa remains at or below 1 for 
            # all values of z
sd.pert=0.01 # governs how far from the perturbation center sa is perturbed
scale=delta/dnorm(0,0,sd.pert) #scaling factor so that perturbation is delta at center of perturbation

# Fig. S3.1C

# plot the purturbations centered on 3 values of z
zz=seq(-.5,.5,.001)
p1=scale*dnorm(zz,-.45,sd.pert)
p2=scale*dnorm(zz,0,sd.pert)
p3=scale*dnorm(zz,.45,sd.pert)
windows()
matplot(zz,cbind(p1,p2,p3),
        #xlab='Environment, z',ylab='Vital rate perturbation',
        xlab='',ylab='',xaxt='n',
        type='l',col=c('red','black','blue'),
        lty='solid',lwd=3,cex.axis=2)

# loop over localized perturbations of the vital rate curve
for(i in 1:length(z)){
  zp=z[i] # center of the perturbation of the vital rate curve
  lambdas=array(0,tmax)
  n=n0
  sas = sa+b*zs+c*zs2+scale*dnorm(zs,zp,sd.pert)  # get sa's for all times  
                                                  # using perturbed vital rate function
  samin.1[i]=min(sas)   # check that sa stays in bounds
  samax.1[i]=max(sas)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  
  lambdas=lambdas[-(1:trash)]  # discard initial lambdas - pre-stationary distribution?
  LLp[i]=mean(log(lambdas)) # compute and store stochastic lambdas for perturbations
  
}

# ****** Now examine antagonistic pleiotropy  ****************

# loop over localized perturbations of the vital rate curve
for(i in 1:length(z)){   # centering perturbations from z=-.5 to z=.4
  zp=z[i] # center of the perturbation of the vital rate curve
  lambdas=array(0,tmax)
  n=n0
  # get sa's for all times using perturbed vital rate function
  # but now add a DOWNWARD perturbation centered at z=0.45; 
  # thus the upward and downward perturbations cancel at z=0.45
  sas = sa+b*zs+c*zs2+scale*(dnorm(zs,zp,sd.pert)-dnorm(zs,0.45,sd.pert))

  samin.2[i]=min(sas)   # check that sa stays in bounds
  samax.2[i]=max(sas)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  
  lambdas=lambdas[-(1:trash)]  # discard initial lambdas - pre-stationary distribution?
  LLap[i]=mean(log(lambdas)) # compute and store stochastic lambdas for perturbations
  
}

# compute and plot the localized selection gradients with and without antagonistic pleiotropy 

dLsdv=(LLp-LL0)/delta  # selection gradients WITHOUT antagonistic pleiotropy
dLsdv.ap=(LLap-LL0)/delta # selection gradients WITH antagonistic pleiotrophy

# Fig. S3.1D

windows()
matplot(z,cbind(dLsdv,dLsdv.ap),type='l',
        #xlab="Environment, z",ylab='Selection gradient',
        xlab='',ylab='',
        col=c('red','black'),lty=c('solid','dashed'),lwd=4,cex.axis=2)
lines(c(.45,.45),c(-.01,.4),lty='dotted',lwd=4)