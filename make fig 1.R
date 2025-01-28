# Make panels of Fig 1 for manuscript entitled 
# "What is demographic lability and when might we expect to see it?"

# Vital rate/driver relationship (VRDR) for adult survival, s.a, vs driver z
# in text: s.a(z) = s.a(0) + b*z + c*z^2
# in code: s.a(z) = sa0 + b*z + c*z^2

rm(list=ls(all=TRUE))
graphics.off()

# First make Fig 1A which plots possible changes to the initially linear
# VRDR for adult survival (sa)

z=seq(0,1,.02)-.5
z2=z^2

sa=.95
b0=.06
c=0  # convex if c>0, concave if c<0

sa.0 = sa+b0*z+c*z2  # quadratic function - linear version

c=0.06  # convex if c>1, concave if c<1
sa.c0 = sa+b0*z+c*z2  # quadratic function - no tradeoff

T1=.08 # medium tradeoff between c and sa
sa.c1 = sa-T1*c+b0*z+c*z2  # quadratic function - tradoff parameter T1

T2=.16 # high tradeoff between c and sa
sa.c2 = sa-T2*c+b0*z+c*z2  # quadratic function - tradeoff parameter T2

c=0
b=b0-.02
sa.b = sa+b*z+c*z2  # quadratic function - linear, less steep

b=b0
sa=sa+.005
sa.sa = sa+b*z+c*z2  # quadratic function - linear, increased intercept

windows()
matplot(z,cbind(sa.0,sa.c0,sa.c1,sa.c2,sa.sa,sa.b),xlab='Environment, z',
        ylab='Adult survival',type="l",ylim=c(.9,1),
        col=c('black','green','green','green','red','blue'),lwd=c(5,4,4,4,4,4),
        lty=c('solid','solid','dotted','dotted','solid','solid'),cex.axis=2)


# Fig 1B  - identical to Figs S4.1B, S6.1 B, and Fig S7.1 B

# compute Beta parameters from desired mean (m) and variance (v)
m=0.5
v=.05 # max v = m*(1-m)
aa=m*(m*(1-m)/v - 1)
bb=(1-m)*(m*(1-m)/v - 1)

windows()
z=seq(-.5,.5,.01)
plot(z,dbeta(z+0.5,aa,bb),xlab="Environment, z",ylab="Density",type="l",lwd=4,cex.axis=2)


###################################################################################
# compute sensitivities to sa, b, and c while varying initial b (overall steepness),
# then make Fig 1C

# baseline values - SLOW LIFE HISTORY
sa=0.95
sj=.1
f=.5

delta=1e-4

tmax=1e6    # time steps to compute lam.s; long times needed to 
            # compute selection gradient for b when b0 is close to zero
trash=1e4   # initial time steps to discard 
n0=matrix(c(0.5,0.5),2,1)  # starting population vector
A=matrix(c(0,sj,f,0),2,2)  # constant parts of projection matrix

# generate Beta-distributed environmental states
zs=rbeta(tmax,aa,bb)
zs=zs-0.5   # make z vary from -1/2 to +1/2`
zs2=zs^2    # only need to do this once

b0s=seq(0,.08,0.002) # initial slopes of the linear VRDR

# allocate memory
LL0=LLc0=LLc1=LLc2=LLsa=LLb=array(0,length(b0s)) # arrays for log lambda.s
samin.0=samin.c0=samin.c1=samin.c2=samin.sa=samin.b=LL0 # arrays to check sa stays in bounds
samax.0=samax.c0=samax.c1=samax.c2=samax.sa=samax.b=LL0

for(i in 1:length(b0s)){    # loop over initial slopes
  
  b0=b0s[i] # current inital slope
  b=b0

  #***** baseline ******
  c=0  # linear if c=0, convex if c>0, concave if c<0

  lambdas=array(0,tmax)
  n=n0
  sas = sa+b*zs+c*zs2   # quadratic function - make all tmax s.a's
  samin.0[i]=min(sas)   # check that sa stays in bounds
  samax.0[i]=max(sas)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  
  lambdas=lambdas[-(1:trash)]  # discard initial lambdas - pre-stationary distribution?
  LL0[i]=mean(log(lambdas))

  #***** for selection gradient to c - curvature parameter ******
  c=delta    # perturb c

  T=0 # no tradeoff between b and sa
  sa.now = sa

  lambdas=array(0,tmax)
  n=n0
  sas = sa.now+b*zs+c*zs2  # quadratic function
  samin.c0[i]=min(sas)  # check that sa stays in bounds
  samax.c0[i]=max(sas)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  
  lambdas=lambdas[-(1:trash)]
  LLc0[i]=mean(log(lambdas))
  
  T1=.08 # medium tradeoff between c and sa
  sa.now = sa - T1*c

  lambdas=array(0,tmax)
  n=n0
  sas = sa.now+b*zs+c*zs2  # quadratic function
  samin.c1[i]=min(sas)  # check that sa stays in bounds
  samax.c1[i]=max(sas)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  
  lambdas=lambdas[-(1:trash)]
  LLc1[i]=mean(log(lambdas))
  
  T2=0.16 # high tradeoff between c and sa
  sa.now = sa - T2*c

  lambdas=array(0,tmax)
  n=n0
  sas = sa.now+b*zs+c*zs2  # quadratic function
  samin.c2[i]=min(sas)  # check that sa stays in bounds
  samax.c2[i]=max(sas)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  
  lambdas=lambdas[-(1:trash)]
  LLc2[i]=mean(log(lambdas))
  
  #***** for selection gradient to sa ******
  c=0           # revert to baseline c
  s2=sa+delta   # now perturb sa

  lambdas=array(0,tmax)
  n=n0
  sas = s2+b*zs+c*zs2  # quadratic function, using s2
  samin.sa[i]=min(sas) # check that sa stays in bounds 
  samax.sa[i]=max(sas)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  
  lambdas=lambdas[-(1:trash)]
  LLsa[i]=mean(log(lambdas))
  
  #***** for selection gradient to b - slope parameter ****** 
  b=b0+delta       # perturb b
  
  lambdas=array(0,tmax)
  n=n0
  sas = sa+b*zs+c*zs2  # quadratic function
  samin.b[i]=min(sas)  # check that sa stays in bounds
  samax.b[i]=max(sas)
  for(t in 1:tmax){
    A[2,2]=sas[t]
    n=A %*% n
    lambdas[t]=sum(n)
    n=n/sum(n)
  }
  
  lambdas=lambdas[-(1:trash)]
  LLb[i]=mean(log(lambdas))

} # for(i in ...

# check samin and samax arrays to be sure that sa stayed in bounds

lam.s = exp(LL0)  # baseline log lambda.s

# compute selection gradients
S.c0=(exp(LLc0)-lam.s)/delta
S.c1=(exp(LLc1)-lam.s)/delta
S.c2=(exp(LLc2)-lam.s)/delta
S.sa=(exp(LLsa)-lam.s)/delta
S.b=(exp(LLb)-lam.s)/delta

Sens=(cbind(S.c0,S.c1,S.c2,S.sa,S.b))

# make Fig 1C
windows()
matplot(b0s,abs(Sens),log='y',xlab='Initial steepness, c0',ylab='|Selection gradient|',
        type='l',col=c('green','green','green','red','blue'),
        lty=c('solid','dashed','dashed','solid','dashed'),lwd=4,cex.axis=2)
