
rm(list=ls())
setwd('/path/to/file/')
regression=T

filez=rev(paste0('',rep(rep(c(''),each=4),times=3),rep(LETTERS[1:3],each=4),rep(1:4,times=3),'_RESULTS.RData'))

bDATA=list(NULL)
for (i in 1:length(filez)){
  load(filez[i])
  bDATA[[i]]=bind.data
}

N.data=matrix(NA,length(bDATA[[1]]$N),length(bDATA))
R.data=N.data; Nd.m.data=N.data; Nd.sd.data=N.data; Nem.data=N.data

for (i in 1:length(bDATA)){
  N.data[,i]=bDATA[[i]]$N
  R.data[,i]=bDATA[[i]]$R
  Nd.m.data[,i]=bDATA[[i]]$Nd.m
  Nd.sd.data[,i]=bDATA[[i]]$Nd.sd
  Nem.data[,i]=bDATA[[i]]$Nem
}
a=1:(length(bDATA)/3)
b=(max(a)+1):(2*length(bDATA)/3)
c=(max(b)+1):length(bDATA)

par(fig=c(0.31,1,0.3,1),mar = c(4.5,4.5,3,1))
plot(NULL,NULL,xlab = 'x (nm)',ylab = 'y (nm)',xlim=c(-0.5*box.dim,0.5*box.dim),ylim=c(-0.5*box.dim,0.5*box.dim),main='Particle Simulation Example',cex.main=2,cex.lab=1.7,cex.axis=1.7)
points(R.xyz[,1],R.xyz[,2],pch="/",col='red',cex=1)
points(N.xyz[,1],N.xyz[,2],pch='--',col='green',cex=2.5)
points(E.xyz[,1],E.xyz[,2],pch='o',col='blue',cex=1)
par(fig=c(0,0.3,0.3,1),mar = rep(0,times=4),new=TRUE)
plot(NULL,NULL,main = '',xlab = '',ylab = '',xaxt='n',yaxt='n',xlim=0:1,ylim=0:1,bty='n')
legend('center',legend = c('Enzyme','RNA','Nucleosome'),col = c('blue','red','green'),fill = c('blue','red','green'),bty = 'n',cex = 2.5)

par(fig=c(0.01,0.5,0.55,1),mar = c(4.5,4.5,3,1))
plot(NULL,NULL,main = 'RNA Occupancy',xlim=c(0,t*1e3),ylim = c(0,1),xlab = 'Time (ms)',ylab = expression(B[R]*'  |Eq. 9.6.2|'),cex.main=1.6,cex.lab=1.7,cex.axis=1.7)
lines(seq(0,t,dt)*1e3,apply(R.data[,a],1,mean,na.rm=TRUE)+apply(R.data[,a],1,sd,na.rm=TRUE),col = 'red',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,a],1,mean,na.rm=TRUE)-apply(R.data[,a],1,sd,na.rm=TRUE),col = 'red',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,b],1,mean,na.rm=TRUE)+apply(R.data[,b],1,sd,na.rm=TRUE),col = 'green',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,b],1,mean,na.rm=TRUE)-apply(R.data[,b],1,sd,na.rm=TRUE),col = 'green',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,c],1,mean,na.rm=TRUE)+apply(R.data[,c],1,sd,na.rm=TRUE),col = 'blue',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,c],1,mean,na.rm=TRUE)-apply(R.data[,c],1,sd,na.rm=TRUE),col = 'blue',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,a],1,mean,na.rm=TRUE),col = 'red',lwd=2,lty='dotted')
lines(seq(0,t,dt)*1e3,apply(R.data[,b],1,mean,na.rm=TRUE),col = 'green',lwd=2,lty='dotted')
lines(seq(0,t,dt)*1e3,apply(R.data[,c],1,mean,na.rm=TRUE),col = 'blue',lwd=4,lty='dotted')
if (regression==T){
  RNA.fit1=nls(B ~ Mx*(1-exp(-k*t)),data = list(t=seq(0,t,dt),B=apply(R.data[,a],1,mean,na.rm=TRUE)),start = list(Mx=max(apply(R.data[,a],1,mean,na.rm=TRUE)),k=(apply(R.data[,a],1,mean,na.rm=TRUE)[5]-apply(R.data[,a],1,mean,na.rm=TRUE)[1])/5/dt))
  lines(seq(0,t,dt)*1e3,predict(RNA.fit1,newdata = seq(0,t,dt)),col = 'red',lwd=4,lty='solid')
}
if (regression==T){
  RNA.fit2=nls(B ~ Mx*(1-exp(-k*t)),data = list(t=seq(0,t,dt),B=apply(R.data[,b],1,mean,na.rm=TRUE)),start = list(Mx=max(apply(R.data[,b],1,mean,na.rm=TRUE)),k=(apply(R.data[,b],1,mean,na.rm=TRUE)[5]-apply(R.data[,b],1,mean,na.rm=TRUE)[1])/5/dt))
  lines(seq(0,t,dt)*1e3,predict(RNA.fit2,newdata = seq(0,t,dt)),col = 'green',lwd=4,lty='solid')
}
if (regression==0){
  RNA.fit3=nls(B ~ Mx*(1-exp(-k*t)),data = list(t=seq(0,t,dt),B=apply(R.data[,c],1,mean,na.rm=TRUE)),start = list(Mx=max(apply(R.data[,c],1,mean,na.rm=TRUE)),k=(apply(R.data[,c],1,mean,na.rm=TRUE)[5]-apply(R.data[,c],1,mean,na.rm=TRUE)[1])/5/dt))
  lines(seq(0,t,dt)*1e3,predict(RNA.fit3,newdata = seq(0,t,dt)),col = 'blue',lwd=4,lty='solid')
}

par(fig=c(0.51,1,0.55,1),mar = c(4.5,4.5,3,1),new=TRUE)
plot(NULL,NULL,main = 'Nucleosome Occupancy',xlim=c(0,t*1e3),ylim = c(0,1),xlab = 'Time (ms)',ylab = expression(B[N]*'  |Eq. 9.6.1|'),cex.main=1.6,cex.lab=1.7,cex.axis=1.7)
lines(seq(0,t,dt)*1e3,apply(N.data[,a],1,mean,na.rm=TRUE)+apply(N.data[,a],1,sd,na.rm=TRUE),col = 'red',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,a],1,mean,na.rm=TRUE)-apply(N.data[,a],1,sd,na.rm=TRUE),col = 'red',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,b],1,mean,na.rm=TRUE)+apply(N.data[,b],1,sd,na.rm=TRUE),col = 'green',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,b],1,mean,na.rm=TRUE)-apply(N.data[,b],1,sd,na.rm=TRUE),col = 'green',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,c],1,mean,na.rm=TRUE)+apply(N.data[,c],1,sd,na.rm=TRUE),col = 'blue',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,c],1,mean,na.rm=TRUE)-apply(N.data[,c],1,sd,na.rm=TRUE),col = 'blue',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,a],1,mean,na.rm=TRUE),col = 'red',lwd=2,lty='dotted')
lines(seq(0,t,dt)*1e3,apply(N.data[,b],1,mean,na.rm=TRUE),col = 'green',lwd=2,lty='dotted')
lines(seq(0,t,dt)*1e3,apply(N.data[,c],1,mean,na.rm=TRUE),col = 'blue',lwd=2,lty='dotted')
if (regression==T){
  NUC.fit1=nls(B ~ Mx*(1-exp(-k*t)),data = list(t=seq(0,t,dt),B=apply(N.data[,a],1,mean,na.rm=TRUE)),start = list(Mx=median(apply(N.data[,a],1,mean,na.rm=TRUE)),k=(apply(N.data[,a],1,mean,na.rm=TRUE)[15]-apply(N.data[,a],1,mean,na.rm=TRUE)[1])/15/dt))
  lines(seq(0,t,dt)*1e3,predict(NUC.fit1,newdata = seq(0,t,dt)),col = 'red',lwd=4,lty='solid')
}
if (regression==T){
  NUC.fit2=nls(B ~ Mx*(1-exp(-k*t)),data = list(t=seq(0,t,dt),B=apply(N.data[,b],1,mean,na.rm=TRUE)),start = list(Mx=median(apply(N.data[,b],1,mean,na.rm=TRUE)),k=(apply(N.data[,b],1,mean,na.rm=TRUE)[15]-apply(N.data[,b],1,mean,na.rm=TRUE)[1])/15/dt))
  lines(seq(0,t,dt)*1e3,predict(NUC.fit2,newdata = seq(0,t,dt)),col = 'green',lwd=4,lty='solid')
}
if (regression==T){
  NUC.fit3=nls(B ~ Mx*(1-exp(-k*t)),data = list(t=seq(0,t,dt),B=apply(N.data[,c],1,mean,na.rm=TRUE)),start = list(Mx=median(apply(N.data[,c],1,mean,na.rm=TRUE)),k=(apply(N.data[,c],1,mean,na.rm=TRUE)[15]-apply(N.data[,c],1,mean,na.rm=TRUE)[1])/15/dt))
  lines(seq(0,t,dt)*1e3,predict(NUC.fit3,newdata = seq(0,t,dt)),col = 'blue',lwd=4,lty='solid')
}

par(fig=c(0.01,0.5,0.1,0.55),mar = c(4.5,4.5,3,1),new=TRUE)
plot(NULL,NULL,main = 'Nucleosome-Protein Proximity',xlim=c(0,t*1e3),ylim = c(0,max(apply(Nd.m.data[,1:3],1,mean,na.rm=TRUE)+3*apply(Nd.m.data[,1:3],1,sd,na.rm=TRUE),na.rm = TRUE)),xlab = 'Time (ms)',ylab = expression(P[EN]*' (nm)'*'  |Eq. 9.7|'),cex.main=1.6,cex.lab=1.7,cex.axis=1.7)
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,a],1,mean,na.rm=TRUE),lty='solid',lwd=0.5,col='red')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,a],1,mean,na.rm=TRUE)+apply(Nd.m.data[,a],1,sd,na.rm=TRUE),lty='dashed',lwd=0.3,col='red')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,a],1,mean,na.rm=TRUE)-apply(Nd.m.data[,a],1,sd,na.rm=TRUE),lty='dashed',lwd=0.3,col='red')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,b],1,mean,na.rm=TRUE),lty='solid',lwd=0.5,col='green')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,b],1,mean,na.rm=TRUE)+apply(Nd.m.data[,b],1,sd,na.rm=TRUE),lty='dashed',lwd=0.3,col='green')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,b],1,mean,na.rm=TRUE)-apply(Nd.m.data[,b],1,sd,na.rm=TRUE),lty='dashed',lwd=0.3,col='green')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,c],1,mean,na.rm=TRUE),col = 'blue',lwd=0.5,lty='solid')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,c],1,mean,na.rm=TRUE)+apply(Nd.m.data[,c],1,sd,na.rm=TRUE),col = 'blue',lwd=0.3,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,c],1,mean,na.rm=TRUE)-apply(Nd.m.data[,c],1,sd,na.rm=TRUE),col = 'blue',lwd=0.3,lty='dashed')

par(fig=c(0.51,1,0.1,0.55),mar = c(4.5,4.5,3,1),new=TRUE)
plot(NULL,NULL,main = 'Effective Molarity',ylab='Relative Likelihood',xlab = expression(M[epsilon]*'  |Eq. 9.8|'),cex.main=1.6,cex.lab=1.7,cex.axis=1.7,xlim=c(0,10),ylim=c(0,1))
lines(density(Nem.data[seq(0,t*1e3,dt*1e3)>=10,a],adjust = 6,from = 0)$x,density(Nem.data[seq(0,t*1e3,dt*1e3)>=10,a],adjust = 6,from = 0)$y/max(density(Nem.data[seq(0,t*1e3,dt*1e3)>=10,a],adjust = 6,from = 0)$y),col = 'red')
lines(density(Nem.data[seq(0,t*1e3,dt*1e3)>=10,b],adjust = 6,from = 0)$x,density(Nem.data[seq(0,t*1e3,dt*1e3)>=10,b],adjust = 6,from = 0)$y/max(density(Nem.data[seq(0,t*1e3,dt*1e3)>=10,b],adjust = 6,from = 0)$y),col = 'green')
lines(density(Nem.data[seq(0,t*1e3,dt*1e3)>=10,c],adjust = 6,from = 0)$x,density(Nem.data[seq(0,t*1e3,dt*1e3)>=10,c],adjust = 6,from = 0)$y/max(density(Nem.data[seq(0,t*1e3,dt*1e3)>=10,c],adjust = 6,from = 0)$y),col = 'blue')
abline(v=mean(Nem.data[seq(0,t*1e3,dt*1e3)>=10,a]),lty='solid',lwd=3,col='red')
abline(v=mean(Nem.data[seq(0,t*1e3,dt*1e3)>=10,b]),lty='solid',lwd=3,col='green')
abline(v=mean(Nem.data[seq(0,t*1e3,dt*1e3)>=10,c]),lty='solid',lwd=3,col='blue')
if (1==2){
  par(fig=c(0.51,1,0.1,0.55),mar = c(4.5,4.5,3,1),new=TRUE)
  plot(NULL,NULL,main = 'Effective Molarity',xlim=c(0,t*1e3),ylim = c(0,0.5*max(apply(Nem.data[,1:3],1,mean,na.rm=TRUE),na.rm = TRUE)),xlab = 'Time (ms)',ylab = expression(M[epsilon]*'  |Eq. 9.8|'),cex.main=1.6,cex.lab=1.7,cex.axis=1.7)
  lines(seq(0,t,dt)*1e3,apply(Nem.data[,a],1,mean,na.rm=TRUE),lty='solid',lwd=0.15,col='red')
  lines(seq(0,t,dt)*1e3,apply(Nem.data[,b],1,mean,na.rm=TRUE),lty='solid',lwd=0.15,col='green')
  lines(seq(0,t,dt)*1e3,apply(Nem.data[,c],1,mean,na.rm=TRUE),lty='solid',lwd=0.15,col='blue')
  abline(h=1,col='grey',lwd=4)
  abline(h=mean(apply(Nem.data[,a],2,mean,na.rm=TRUE)),lty='solid',lwd=3,col='red')
  abline(h=mean(apply(Nem.data[,a],2,mean,na.rm=TRUE))+sd(apply(Nem.data[,a],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='red')
  abline(h=mean(apply(Nem.data[,a],2,mean,na.rm=TRUE))-sd(apply(Nem.data[,a],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='red')
  abline(h=mean(apply(Nem.data[,b],2,mean,na.rm=TRUE)),lty='solid',lwd=3,col='green')
  abline(h=mean(apply(Nem.data[,b],2,mean,na.rm=TRUE))+sd(apply(Nem.data[,b],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='green')
  abline(h=mean(apply(Nem.data[,b],2,mean,na.rm=TRUE))-sd(apply(Nem.data[,b],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='green')
  abline(h=mean(apply(Nem.data[,c],2,mean,na.rm=TRUE)),lty='solid',lwd=3,col='blue')
  abline(h=mean(apply(Nem.data[,c],2,mean,na.rm=TRUE))+sd(apply(Nem.data[,c],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='blue')
  abline(h=mean(apply(Nem.data[,c],2,mean,na.rm=TRUE))-sd(apply(Nem.data[,c],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='blue')
}

par(fig=c(0,0.33,0,0.1),mar = rep(0,times=4),new=TRUE)
plot(NULL,NULL,main = '',xlab = '',ylab = '',xaxt='n',yaxt='n',xlim=0:1,ylim=0:1,bty='n')
legend('center',legend = c(expression(K[dR]*' = 1 M')),col = c('blue'),fill = c('blue'),cex=2.5,bty = 'n')
par(fig=c(0.33,0.67,0,0.1),mar = rep(0,times=4),new=TRUE)
plot(NULL,NULL,main = '',xlab = '',ylab = '',xaxt='n',yaxt='n',xlim=0:1,ylim=0:1,bty='n')
legend('center',legend = c(expression(K[dR]*' = 100 nM')),col = c('green'),fill = c('green'),cex=2.5,bty = 'n')
par(fig=c(0.67,1,0,0.1),mar = rep(0,times=4),new=TRUE)
plot(NULL,NULL,main = '',xlab = '',ylab = '',xaxt='n',yaxt='n',xlim=0:1,ylim=0:1,bty='n')
legend('center',legend = c(expression(K[dR]*' = 1 nM')),col = c('red'),fill = c('red'),cex=2.5,bty = 'n')
