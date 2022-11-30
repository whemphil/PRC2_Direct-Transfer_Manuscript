
rm(list=ls())
setwd('/path/to/output/')

filez=paste0('',rep(rep(c(''),each=4),times=3),rep(LETTERS[1:3],each=4),rep(1:4,times=3),'_RESULTS.RData')

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
points(N.xyz[,1],N.xyz[,2],pch='--',col='blue',cex=2.5)
points(E.xyz[,1],E.xyz[,2],pch='o',col='black',cex=1)
par(fig=c(0,0.3,0.3,1),mar = rep(0,times=4),new=TRUE)
plot(NULL,NULL,main = '',xlab = '',ylab = '',xaxt='n',yaxt='n',xlim=0:1,ylim=0:1,bty='n')
legend('center',legend = c('Enzyme','RNA','Nucleosome'),col = c('black','red','blue'),fill = c('black','red','blue'),bty = 'n',cex = 2.5)

par(fig=c(0.51,1,0.55,1),mar = c(4.5,4.5,3,1))
plot(NULL,NULL,main = 'Nucleosome Occupancy',xlim=c(0,t*1e3),ylim = c(0,1),xlab = 'Time (ms)',ylab = 'Fraction Bound',cex.main=2,cex.lab=1.7,cex.axis=1.7)
lines(seq(0,t,dt)*1e3,apply(N.data[,a],1,mean,na.rm=TRUE),col = 'red',lwd=3,lty='solid')
lines(seq(0,t,dt)*1e3,apply(N.data[,a],1,mean,na.rm=TRUE)+apply(N.data[,a],1,sd,na.rm=TRUE),col = 'red',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,a],1,mean,na.rm=TRUE)-apply(N.data[,a],1,sd,na.rm=TRUE),col = 'red',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,b],1,mean,na.rm=TRUE),col = 'blue',lwd=3,lty='solid')
lines(seq(0,t,dt)*1e3,apply(N.data[,b],1,mean,na.rm=TRUE)+apply(N.data[,b],1,sd,na.rm=TRUE),col = 'blue',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,b],1,mean,na.rm=TRUE)-apply(N.data[,b],1,sd,na.rm=TRUE),col = 'blue',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,c],1,mean,na.rm=TRUE),col = 'green',lwd=3,lty='solid')
lines(seq(0,t,dt)*1e3,apply(N.data[,c],1,mean,na.rm=TRUE)+apply(N.data[,c],1,sd,na.rm=TRUE),col = 'green',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(N.data[,c],1,mean,na.rm=TRUE)-apply(N.data[,c],1,sd,na.rm=TRUE),col = 'green',lwd=0.5,lty='dashed')

par(fig=c(0.01,0.5,0.55,1),mar = c(4.5,4.5,3,1),new=TRUE)
plot(NULL,NULL,main = 'RNA Occupancy',xlim=c(0,t*1e3),ylim = c(0,1),xlab = 'Time (ms)',ylab = 'Fraction Bound',cex.main=2,cex.lab=1.7,cex.axis=1.7)
lines(seq(0,t,dt)*1e3,apply(R.data[,a],1,mean,na.rm=TRUE),col = 'red',lwd=3,lty='solid')
lines(seq(0,t,dt)*1e3,apply(R.data[,a],1,mean,na.rm=TRUE)+apply(R.data[,a],1,sd,na.rm=TRUE),col = 'red',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,a],1,mean,na.rm=TRUE)-apply(R.data[,a],1,sd,na.rm=TRUE),col = 'red',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,b],1,mean,na.rm=TRUE),col = 'blue',lwd=3,lty='solid')
lines(seq(0,t,dt)*1e3,apply(R.data[,b],1,mean,na.rm=TRUE)+apply(R.data[,b],1,sd,na.rm=TRUE),col = 'blue',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,b],1,mean,na.rm=TRUE)-apply(R.data[,b],1,sd,na.rm=TRUE),col = 'blue',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,c],1,mean,na.rm=TRUE),col = 'green',lwd=3,lty='solid')
lines(seq(0,t,dt)*1e3,apply(R.data[,c],1,mean,na.rm=TRUE)+apply(R.data[,c],1,sd,na.rm=TRUE),col = 'green',lwd=0.5,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(R.data[,c],1,mean,na.rm=TRUE)-apply(R.data[,c],1,sd,na.rm=TRUE),col = 'green',lwd=0.5,lty='dashed')

par(fig=c(0.01,0.5,0.1,0.55),mar = c(4.5,4.5,3,1),new=TRUE)
plot(NULL,NULL,main = 'Nucleosome-Protein Proximity',xlim=c(0,t*1e3),ylim = c(0,max(apply(Nd.m.data[,1:3],1,mean,na.rm=TRUE)+3*apply(Nd.m.data[,1:3],1,sd,na.rm=TRUE),na.rm = TRUE)),xlab = 'Time (ms)',ylab = 'Nearest Proteins (nm)',cex.main=1.8,cex.lab=1.7,cex.axis=1.7)
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,a],1,mean,na.rm=TRUE),lty='solid',lwd=0.5,col='red')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,a],1,mean,na.rm=TRUE)+apply(Nd.m.data[,a],1,sd,na.rm=TRUE),lty='dashed',lwd=0.3,col='red')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,a],1,mean,na.rm=TRUE)-apply(Nd.m.data[,a],1,sd,na.rm=TRUE),lty='dashed',lwd=0.3,col='red')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,b],1,mean,na.rm=TRUE),lty='solid',lwd=0.5,col='blue')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,b],1,mean,na.rm=TRUE)+apply(Nd.m.data[,b],1,sd,na.rm=TRUE),lty='dashed',lwd=0.3,col='blue')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,b],1,mean,na.rm=TRUE)-apply(Nd.m.data[,b],1,sd,na.rm=TRUE),lty='dashed',lwd=0.3,col='blue')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,c],1,mean,na.rm=TRUE),col = 'green',lwd=0.5,lty='solid')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,c],1,mean,na.rm=TRUE)+apply(Nd.m.data[,c],1,sd,na.rm=TRUE),col = 'green',lwd=0.3,lty='dashed')
lines(seq(0,t,dt)*1e3,apply(Nd.m.data[,c],1,mean,na.rm=TRUE)-apply(Nd.m.data[,c],1,sd,na.rm=TRUE),col = 'green',lwd=0.3,lty='dashed')

par(fig=c(0.51,1,0.1,0.55),mar = c(4.5,4.5,3,1),new=TRUE)
plot(NULL,NULL,main = 'Effective Molarity',xlim=c(0,t*1e3),ylim = c(0,0.5*max(apply(Nem.data[,1:3],1,mean,na.rm=TRUE),na.rm = TRUE)),xlab = 'Time (ms)',ylab = 'Relative [E]',cex.main=1.8,cex.lab=1.7,cex.axis=1.7)
lines(seq(0,t,dt)*1e3,apply(Nem.data[,a],1,mean,na.rm=TRUE),lty='solid',lwd=0.15,col='red')
lines(seq(0,t,dt)*1e3,apply(Nem.data[,b],1,mean,na.rm=TRUE),lty='solid',lwd=0.15,col='blue')
lines(seq(0,t,dt)*1e3,apply(Nem.data[,c],1,mean,na.rm=TRUE),lty='solid',lwd=0.15,col='green')
abline(h=1,col='grey',lwd=4)
abline(h=mean(apply(Nem.data[,a],2,mean,na.rm=TRUE)),lty='solid',lwd=3,col='red')
abline(h=mean(apply(Nem.data[,a],2,mean,na.rm=TRUE))+sd(apply(Nem.data[,a],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='red')
abline(h=mean(apply(Nem.data[,a],2,mean,na.rm=TRUE))-sd(apply(Nem.data[,a],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='red')
abline(h=mean(apply(Nem.data[,b],2,mean,na.rm=TRUE)),lty='solid',lwd=3,col='blue')
abline(h=mean(apply(Nem.data[,b],2,mean,na.rm=TRUE))+sd(apply(Nem.data[,b],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='blue')
abline(h=mean(apply(Nem.data[,b],2,mean,na.rm=TRUE))-sd(apply(Nem.data[,b],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='blue')
abline(h=mean(apply(Nem.data[,c],2,mean,na.rm=TRUE)),lty='solid',lwd=3,col='green')
abline(h=mean(apply(Nem.data[,c],2,mean,na.rm=TRUE))+sd(apply(Nem.data[,c],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='green')
abline(h=mean(apply(Nem.data[,c],2,mean,na.rm=TRUE))-sd(apply(Nem.data[,c],2,mean,na.rm=TRUE)),lty='dashed',lwd=1,col='green')

par(fig=c(0,0.33,0,0.1),mar = rep(0,times=4),new=TRUE)
plot(NULL,NULL,main = '',xlab = '',ylab = '',xaxt='n',yaxt='n',xlim=0:1,ylim=0:1,bty='n')
legend('center',legend = c(expression(K[dR]*' = 5 M')),col = c('blue'),fill = c('blue'),cex=2.5,bty = 'n')
par(fig=c(0.33,0.67,0,0.1),mar = rep(0,times=4),new=TRUE)
plot(NULL,NULL,main = '',xlab = '',ylab = '',xaxt='n',yaxt='n',xlim=0:1,ylim=0:1,bty='n')
legend('center',legend = c(expression(K[dR]*' = 500 nM')),col = c('green'),fill = c('green'),cex=2.5,bty = 'n')
par(fig=c(0.67,1,0,0.1),mar = rep(0,times=4),new=TRUE)
plot(NULL,NULL,main = '',xlab = '',ylab = '',xaxt='n',yaxt='n',xlim=0:1,ylim=0:1,bty='n')
legend('center',legend = c(expression(K[dR]*' = 5 nM')),col = c('red'),fill = c('red'),cex=2.5,bty = 'n')
