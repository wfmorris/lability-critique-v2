# This code accompanies a manuscript entitled 
# "What is demographic lability and when might we expect to see it?"
# It uses data obtained by digitizing Fig. 2 in Barraquand et al., 2014, 
# Journal of Animal Ecology. The data (found in the file 
# "skua data from Barraquand.csv") describes the relationship between fecundity
# ("productivity") of long-tailed skua and abundance of lemmings at two sites
# in Greenland. This script fits 3 numerical response curves (linear, Monod, and
# logistic) using maximum likelihood (assuming Normally-distributed errors) and 
# then uses AICc (Hurvitch and Tsai 1989) to identify the most parsimonious 
# curve. The 3 curves are fit to data for each site separately and for the two
# sites combined. The script also produces fig. S2.1 in appendix S2 in the manuscript.  

graphics.off()

data=read.csv("skua data from Barraquand.csv")

p=c(4,3,3,2) # number of parameters in the 3 models to fit, including residual variance
logL=AICc=array(0,4)

# First fit functions using only data from the Traill site

data.t=subset(data,site=="Traill")
attach(data.t)

q=15 # number of data points

# sum of squares function for the logistic used by Barraquand
SSlogist=function(p) sum( (fledglings - (p[1]*(1-1/(1+(lemmings/p[2])^p[3]))) )^2 ) # p is the vector of parameter values
out.log.t=nlm(SSlogist,c(1.75,6,3),gradtol=1e-7)	# perform the nonlinear minimization, using Barraquand's values as starting params 
params.log.t=out.log.t$estimate; params.log.t	# retrieve the parameter estimates
V.log.t=out.log.t$minimum/q;		# convert the minimum sum of squares to the residual variance for 
                                # the logistic model

# sum of squares function for Monod function
SSmonod=function(p) sum( (fledglings - p[1]*lemmings/(p[2]+lemmings) )^2 ) # p is the vector of parameter values
out.mon.t=nlm(SSmonod,c(1.6,3),gradtol=1e-7)	# perform the nonlinear minimization, using Barraquand's values as starting params 
params.mon.t=out.mon.t$estimate	# retrieve the parameter estimates
V.mon.t=out.mon.t$minimum/q;		# convert the minimum sum of squares to the residual variance for 
                                # the Monod model

# linear function - non-zero intercept
out.lin1.t=lm(fledglings~lemmings)
summary(out.lin1.t)
params.lin1.t=out.lin1.t$coef
V.lin1.t = sum( (fledglings - (params.lin1.t[1]+params.lin1.t[2]*lemmings))^2 )/q

# linear function - zero intercept
out.lin0.t=lm(fledglings~0+lemmings)
summary(out.lin0.t)
params.lin0.t=out.lin0.t$coef
V.lin0.t = sum( (fledglings - params.lin0.t*lemmings)^2 )/q

np=4    # logistic model has 4 parameters, including residual variance V.r
LL=-.5*q*(log(2*pi*V.log.t)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[1]=LL
AICc[1]=AIC

np=3    # Monod model has 3 parameters, including V.r
LL=-.5*q*(log(2*pi*V.mon.t)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[2]=LL
AICc[2]=AIC

np=3    # non-zero intercept linear model has 3 parameters, including V.r
LL=-.5*q*(log(2*pi*V.lin1.t)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[3]=LL
AICc[3]=AIC

np=2    # linear (zero intercept) model has 2 parameters, including V.r
LL=-.5*q*(log(2*pi*V.lin0.t)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[4]=LL
AICc[4]=AIC

# zero-intercept linear model has lowest AICc
results=data.frame(Model=c("Logistic","Monod","Linear1","Linear0"),p,logL,AICc)
results.traill = results


# Next fit functions using only data from the Zackenberg site

data.z=subset(data,site=="Zackenberg")
attach(data.z)

q=15  # number of data points

# sum of squares function for the logistic used by Barraquand
SSlogist=function(p) sum( (fledglings - (p[1]*(1-1/(1+(lemmings/p[2])^p[3]))) )^2 ) # p is the vector of parameter values
out.log.z=nlm(SSlogist,c(1.75,6,3),gradtol=1e-7)	# perform the nonlinear minimization, using Barraquand's values as starting params 
params.log.z=out.log.z$estimate; # retrieve the parameter estimates
V.log.z=out.log.z$minimum/q;		# convert the minimum sum of squares to the residual variance for 
                                # the logistic model

# sum of squares function for Monod function
SSmonod=function(p) sum( (fledglings - p[1]*lemmings/(p[2]+lemmings) )^2 ) # p is the vector of parameter values
out.mon.z=nlm(SSmonod,c(1.6,3),gradtol=1e-7)	# perform the nonlinear minimization, using Barraquand's values as starting params 
params.mon.z=out.mon.z$estimate # retrieve the parameter estimates
V.mon.z=out.mon.z$minimum/q;		# convert the minimum sum of squares to the residual variance for 
                                # the Monod model

# linear function - non-zero intercept
out.lin1.z=lm(fledglings~lemmings)
summary(out.lin1.z)
params.lin1.z=out.lin1.z$coef
V.lin1.z = sum( (fledglings - (params.lin1.z[1]+params.lin1.z[2]*lemmings))^2 )/q

# linear function - zero intercept
out.lin0.z=lm(fledglings~0+lemmings)
summary(out.lin0.z)
params.lin0.z=out.lin0.z$coef
V.lin0.z = sum( (fledglings - params.lin0.z*lemmings)^2 )/q

np=4    # logistic model has 4 parameters, including the resicual variance V.r
LL=-.5*q*(log(2*pi*V.log.z)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[1]=LL
AICc[1]=AIC

np=3    # Monod model has 3 parameters, including V.r
LL=-.5*q*(log(2*pi*V.mon.z)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[2]=LL
AICc[2]=AIC

np=3    # linear non-zero model has 3 parameters, including V.r
LL=-.5*q*(log(2*pi*V.lin1.z)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[3]=LL
AICc[3]=AIC

np=2    # linear (zero intercept) model has 2 parameters, including V.r
LL=-.5*q*(log(2*pi*V.lin0.z)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[4]=LL
AICc[4]=AIC

# zero-intercept linear model again has lowest AICc
results=data.frame(Model=c("Logistic","Monod","Linear1","Linear0"),p,logL,AICc)
results.zack = results


# Finally fit functions using data from both sites combined

attach(data)

q=30  # number of data points

# sum of squares function for the logistic used by Barraquand
SSlogist=function(p) sum( (fledglings - (p[1]*(1-1/(1+(lemmings/p[2])^p[3]))) )^2 ) # p is the vector of parameter values
out.log.a=nlm(SSlogist,c(1.75,6,3),gradtol=1e-7)	# perform the nonlinear minimization, using Barraquand's values as starting params 
params.log.a=out.log.a$estimate; # retrieve the parameter estimates
V.log.a=out.log.a$minimum/q;		# convert the minimum sum of squares to the residual variance for 
                                # the logistic model

# sum of squares function for Monod function
SSmonod=function(p) sum( (fledglings - p[1]*lemmings/(p[2]+lemmings) )^2 ) # p is the vector of parameter values
out.mon.a=nlm(SSmonod,c(1.6,3),gradtol=1e-7)	# perform the nonlinear minimization, using Barraquand's values as starting params 
params.mon.a=out.mon.a$estimate # retrieve the parameter estimates
V.mon.a=out.mon.a$minimum/q;		# convert the minimum sum of squares to the residual variance for 
                                # the Monod model

# linear function - non-zero intercept
out.lin1.a=lm(fledglings~lemmings)
summary(out.lin1.a)
params.lin1.a=out.lin1.a$coef
V.lin1.a = sum( (fledglings - (params.lin1.a[1]+params.lin1.a[2]*lemmings))^2 )/q

# linear function - zero intercept
out.lin0.a=lm(fledglings~0+lemmings)
summary(out.lin0.a)
params.lin0.a=out.lin0.a$coef
V.lin0.a = sum( (fledglings - params.lin0.a*lemmings)^2 )/q

np=4    # logistic model has 4 parameters, including the resicual variance V.r
LL=-.5*q*(log(2*pi*V.log.a)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[1]=LL
AICc[1]=AIC

np=3    # Monod model has 3 parameters, including V.r
LL=-.5*q*(log(2*pi*V.mon.a)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[2]=LL
AICc[2]=AIC

np=3    # nonzero intercept linear model has 3 parameters, including V.r
LL=-.5*q*(log(2*pi*V.lin1.a)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[3]=LL
AICc[3]=AIC

np=2    # linear (zero intercept) model has 2 parameters, including V.r
LL=-.5*q*(log(2*pi*V.lin0.a)+1)	# log likelihood
AIC=-2*LL + (2*np*q/(q-np-1))
logL[4]=LL
AICc[4]=AIC

# zero-intercept linear model again has lowest AICc
results=data.frame(Model=c("Logistic","Monod","Linear1","Linear0"),p,logL,AICc)
results.all = results

# make Fig S2.1

windows()
plot(data.t$lemmings,data.t$fledglings,type='p',ylim=c(0,1.8),xlim=c(0,15),
     xlab="Lemmings per hectare",ylab="Average fledglings per nest",
     pch=19,col='red',cex=1.5,cex.axis=1.5,cex.lab=1.5)
points(data.z$lemmings,data.z$fledglings,pch=1,col='blue',cex=1.5)
xs=seq(0,15,0.2)
ys=params.lin0.a*xs
lines(xs,ys,col='black',lwd=3)
yb=1.75*(1-1/(1+(xs/6)^3)) # line fitted by eye by Barraquand et al
lines(xs,yb,col='green',lty='dotted',lwd=3)
