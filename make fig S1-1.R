# make panels of Fig S1.1 for manuscript entitled 
# "What is adaptive demographic lability and when might we expect to see it?"

graphics.off()

Mu.c=Mu.v=array(0,9)

# panel a

# beta mean and variance of driver z
M=0.7  # left-skewed
V.max = M*(1-M) # constraint on V for the beta
V=0.03

# use desired mean M and variance V to get shape parameters for the beta 
a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)

z=seq(0,1,.002) # driver levels

A=2.5
B=2  # convex VRDR
vr=A*z^B 

drivers=rbeta(1E7,a,b) # sample from driver distribution
z.med=median(drivers)  # get driver median
Mu.c[1] = A*z.med^B         # vital rate at median of driver
Mu.v[1] = mean(A*drivers^B) # mean vital rate when driver varies

windows()
matplot(z,cbind(dbeta(z,a,b),vr),type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.8,lty='solid',xaxt='n')
axis(1,labels=F)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3)

# panel b

M=0.5  # symmetrical
V.max = M*(1-M)
V=0.03

a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)

z=seq(0,1,.002)

A=2.5
B=2
vr=A*z^B

drivers=rbeta(1E7,a,b)
z.med=median(drivers)
Mu.c[2] = A*z.med^B
Mu.v[2] = mean(A*drivers^B)

windows()
matplot(z,cbind(dbeta(z,a,b),vr),type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.8,lty='solid',xaxt='n',yaxt='n')
axis(1,labels=F)
axis(2,labels=F)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3)

# panel c

M=0.3  # right-skewed
V.max = M*(1-M)
V=0.03

a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)

z=seq(0,1,.002)

A=2.5
B=2
vr=A*z^B

drivers=rbeta(1E7,a,b)
z.med=median(drivers)
Mu.c[3] = A*z.med^B
Mu.v[3] = mean(A*drivers^B)

windows()
matplot(z,cbind(dbeta(z,a,b),vr),type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.8,lty='solid',xaxt='n',yaxt='n')
axis(1,labels=F)
axis(2,labels=F)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3)

# panel d

M=0.7  # left-skewed
V.max = M*(1-M)
V=0.03

a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)

z=seq(0,1,.002)

A=2.5
B=1  # VRDR linear
vr=A*z^B

drivers=rbeta(1E7,a,b)
z.med=median(drivers)
Mu.c[4] = A*z.med^B
Mu.v[4] = mean(A*drivers^B)

windows()
matplot(z,cbind(dbeta(z,a,b),vr),type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.8,lty='solid',xaxt='n')
axis(1,labels=F)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3)

# panel e

M=0.5  # symmetrical
V.max = M*(1-M)
V=0.03

a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)

z=seq(0,1,.002)

A=2.5
B=1
vr=A*z^B

drivers=rbeta(1E7,a,b)
z.med=median(drivers)
Mu.c[5] = A*z.med^B
Mu.v[5] = mean(A*drivers^B)

windows()
matplot(z,cbind(dbeta(z,a,b),vr),type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.8,lty='solid',xaxt='n',yaxt='n')
axis(1,labels=F)
axis(2,labels=F)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3)

# panel f

M=0.3  # right-skewed
V.max = M*(1-M)
V=0.03

a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)

z=seq(0,1,.002)

A=2.5
B=1
vr=A*z^B

drivers=rbeta(1E7,a,b)
z.med=median(drivers)
Mu.c[6] = A*z.med^B
Mu.v[6] = mean(A*drivers^B)

windows()
matplot(z,cbind(dbeta(z,a,b),vr),type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.8,lty='solid',xaxt='n',yaxt='n')
axis(1,labels=F)
axis(2,labels=F)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3)

# panel g

M=0.7  # left-skewed
V.max = M*(1-M)
V=0.03

a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)

z=seq(0,1,.002)

A=2.5
B=0.6  # VRDR concave
vr=A*z^B

drivers=rbeta(1E7,a,b)
z.med=median(drivers)
Mu.c[7] = A*z.med^B
Mu.v[7] = mean(A*drivers^B)

windows()
matplot(z,cbind(dbeta(z,a,b),vr),type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.8,lty='solid')
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3)

# panel h

M=0.5  # symmetrical
V.max = M*(1-M)
V=0.03

a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)

z=seq(0,1,.002)

A=2.5
B=0.6
vr=A*z^B

drivers=rbeta(1E7,a,b)
z.med=median(drivers)
Mu.c[8] = A*z.med^B
Mu.v[8] = mean(A*drivers^B)

windows()
matplot(z,cbind(dbeta(z,a,b),vr),type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.8,lty='solid',yaxt='n')
axis(2,labels=F)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3)

# panel i

M=0.3 # right-skewed
V.max = M*(1-M)
V=0.03

a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)

z=seq(0,1,.002)

A=2.5
B=0.6
vr=A*z^B

drivers=rbeta(1E7,a,b)
z.med=median(drivers)
Mu.c[9] = A*z.med^B
Mu.v[9] = mean(A*drivers^B)

windows()
matplot(z,cbind(dbeta(z,a,b),vr),type='l',
        xlab='',ylab='',lwd=3,cex.axis=2.8,lty='solid',yaxt='n')
axis(2,labels=F)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3)

# compute lability index for all panels
J=(Mu.v-Mu.c)/Mu.c
