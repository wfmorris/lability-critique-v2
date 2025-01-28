# make panels of Fig S7.3 for manuscript entitled 
# "What is demographic lability and when might we expect to see it?"

rm(list=ls(all=TRUE))
graphics.off()
 
z=seq(0,10,0.1)  # environmental states

# panel a

a=.1
a=log(1/a - 1)
b=.8
c=2
vr=1/(1+exp(a-b*(z-c))) # logistic VRDR

windows()
plot(z,vr,type='l',xlab='',ylab='',
     col='black',lwd=4,lty='solid',cex.axis=2,ylim=c(0,1),xaxt='n')

# panel b

M=3  # mean z
Vs=seq(0.2,4.2,2)  # values of var(z)
ps=matrix(0,length(z),3)
for(i in 1:3){
  V=Vs[i]
  alpha=M^2/V  # get gamma parameters
  sigma=V/M
  p=dgamma(z,shape=alpha,scale=sigma) # compute 3 gamma pdf's, then plot them
  ps[,i]=p
}
windows()
matplot(z,ps,type='l',xlab='',ylab='',
        col=c('blue','black','red'),lwd=4,lty='solid',cex.axis=2,ylim=c(0,1)) 

# panel c

# increase variance of gamma and plot mean (black) and var (red) of VR
M=3
Vs=seq(0.2,4.2,.1)
VR.mean=VR.sd=array(0,length(Vs))
for(i in 1:length(Vs)){
  V=Vs[i]
  alpha=M^2/V
  sigma=V/M
  zs=rgamma(1E6,shape=alpha,scale=sigma) # sample from gamma distribution of z
  vrs=1/(1+exp(a-b*(zs-c)))
  VR.mean[i]=mean(vrs)
  VR.sd[i]=sd(vrs)
}
windows()
matplot(Vs,cbind(VR.mean,VR.sd),type='l',xlab='',ylab='',
        col=c('black','red'),lwd=4,lty='solid',cex.axis=2,ylim=c(0,.25))