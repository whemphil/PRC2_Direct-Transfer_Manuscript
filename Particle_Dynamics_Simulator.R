##### SCRIPT INFORMATION

# Particle_Dynamics_Simulator.R

# A script to simulate particle dynamics to test the 'kinetic sink' hypothesis

# Associated with manuscript: "PRC2 'Hand-Off' From G-quadruplex RNA to dsDNA: Implications For RNA-Binding Chromatin Modifiers"

##########################

##### DEVELOPER NOTES

# The manuscript this code is associated with only utilizes variation in the 'Environment' and 'Key variables' sections of SCRIPT USER INPUT

##########################

##### SCRIPT USER INPUT

rm(list = ls())

### Environment

setwd('~/path/to/output/')

### Key Variables

# simulation parameters
t=50e-3 # simulation length; s
tau=1e-8 # time step; s
box.dim=1e3 # length/width/height of simulation box; nm

# reaction parameters
KdN=150e-9 # nucleosome dissociation constant; M
KdR=2.5e-9 # RNA dissociation constant; M
E=300 # enzyme concentration; nM
N=15 # nucleosome concentration; nM
RpN=8 # RNAs tethered per nucleosome

### Extraneous Variables

# simulation parameters
dtt=1e-4 # backup frequency; s
dt=1e-5 # save frequency; s

# reaction parameters
D=100e6 # diffusion coefficient; nm^2/s
p.d=5 # particle diameter; nm
AvN=6.022e23 # Avogadro's number

# save parameters
sim.name='SIM-null' # simulation prefix to use for save files

# simulation parameters
start.file=NULL # name of a save file to restart a simulation from

##########################

##### BEGIN SCRIPT

### Semi-Autonomous Script Functions

# additional parameter calculation
D.step=sqrt(6*tau*D) # average diffusion radius per time step; nm
k1=4e-24*pi*D*p.d*AvN # particle association rate constant; 1/M/s
kn1N=k1*KdN # nucleosome dissociation rate constant; 1/s
kn1R=k1*KdR # RNA dissociation rate constant; 1/s

# save parameters
parameters=mget(ls())

# custom functions
diffuse <- function(cart.coord,mean.walk){
  n=length(cart.coord)/3
  theta=runif(n,0,2*pi)
  phi=runif(n,0,2*pi)
  r=rexp(n,1/mean.walk)
  x=r*sin(phi)*cos(theta)
  y=r*sin(phi)*sin(theta)
  z=r*cos(phi)
  new.coord=cart.coord+cbind(x,y,z)
  return(new.coord)
}
wrap.box <- function(box.dim,cart.coord){
  boundary=0.5*box.dim
  breaches=(abs(cart.coord)>(boundary))
  wrapped=cart.coord
  wrapped[breaches]=cart.coord[breaches]-(boundary*cart.coord[breaches]/abs(cart.coord[breaches]))+(-1*boundary*cart.coord[breaches]/abs(cart.coord[breaches]))
  return(wrapped)
}
r.calc <- function(ref,E.xyz){
  temp.1=(E.xyz-matrix(rep(ref,times=nrow(E.xyz)),byrow = TRUE,ncol = 3))^2
  r=sqrt(apply(temp.1,1,sum))
  return(r)
}
bind.check <- function(ref,E.xyz,p.d){
  temp.1=(E.xyz-matrix(rep(ref,times=nrow(E.xyz)),byrow = TRUE,ncol = 3))^2
  r=sqrt(apply(temp.1,1,sum))
  bound=(r==min(r) & r<=p.d)
  if (max(bound)==1){
    id=which.max(bound)
  } else {
    id=NA
  }
  return(id)
}
binding.update <- function(E.xyz,N.xyz,R.xyz,bind.counter,kn1N,kn1R,tau,p.d,ex.bind=TRUE){
  new.count=bind.counter
  R.id=apply(R.xyz,1,bind.check,E.xyz=E.xyz,p.d = p.d)
  if (ex.bind==TRUE){
    R.id[(R.id %in% which(new.count[['E']]!=0))]=NA
  }
  if (sum(is.na(R.id)==FALSE & bind.counter[['R']]==0)>0){
    R.res=round(rexp(sum(is.na(R.id)==FALSE & bind.counter[['R']]==0),kn1R)/tau)
    new.count[['R']][(is.na(R.id)==FALSE & bind.counter[['R']]==0)]=R.res
    new.count[['E']][(R.id[((is.na(R.id)==FALSE) & (bind.counter[['R']]==0))])]=R.res
  }
  N.id=apply(N.xyz,1,bind.check,E.xyz=E.xyz,p.d = p.d)
  if (ex.bind==TRUE){
    N.id[(N.id %in% which(new.count[['E']]!=0))]=NA
  }
  if (sum(is.na(N.id)==FALSE & bind.counter[['N']]==0)>0){
    N.res=round(rexp(sum(is.na(N.id)==FALSE & bind.counter[['N']]==0),kn1N)/tau)
    new.count[['N']][(is.na(N.id)==FALSE & bind.counter[['N']]==0)]=N.res
    new.count[['E']][(N.id[((is.na(N.id)==FALSE) & (bind.counter[['N']]==0))])]=N.res
  }
  for (i in 1:length(new.count)){
    new.count[[i]][new.count[[i]]>0]=new.count[[i]][new.count[[i]]>0]-1
  }
  return(new.count)
}
skim.r <- function(data,RpN){
  temp.1=data[order(data)]
  temp.2=temp.1[1:round(RpN/2)]
  temp.3=mean(temp.2)
  return(temp.3)
}
em.calc <- function(data,bind.counter,N.n,p.d,box.dim){
  t.1=sum(data<=(10*p.d))/N.n
  t.2=4/3*pi*(10*p.d)^3
  t.3=box.dim^3
  t.4a=t.1/t.2
  t.4b=sum(bind.counter[['E']]==0)/t.3
  EM=t.4a/t.4b
  return(EM)
}
ETA.calc <- function(EST,progress){
  ETA = (Sys.time()-EST)*(1-progress)/progress + Sys.time()
  return(ETA)
}

# accessory packages
library(lubridate)

### Autonomous Script Functions

show(paste0('Start Time: ',Sys.time()))

# prepare initial conditions
timer=0
E.n=round(AvN*(E*1e-9)*(box.dim*1e-8)^3)
N.n=round(AvN*(N*1e-9)*(box.dim*1e-8)^3)
N.xyz=array(runif(N.n*3,-0.45*box.dim,0.45*box.dim),c(N.n,3))
R.xyz=matrix(rep(c(t(N.xyz)),times=RpN),ncol = 3,byrow = TRUE);R.xyz[,1]=R.xyz[,1]+rep(c((-1*round(RpN/2)):-1,1:(RpN-round(RpN/2)))*(2*p.d),each=N.n)
R.xyz=wrap.box(box.dim,R.xyz)
N.check=(sum(apply(rbind(N.xyz,R.xyz),1,r.calc,E.xyz=rbind(N.xyz,R.xyz))<(2*p.d) & apply(rbind(N.xyz,R.xyz),1,r.calc,E.xyz=rbind(N.xyz,R.xyz))>0)>0);while (N.check==TRUE){N.xyz=array(runif(N.n*3,-0.5*box.dim,0.5*box.dim),c(N.n,3));R.xyz=matrix(rep(c(t(N.xyz)),times=RpN),ncol = 3,byrow = TRUE);R.xyz[,1]=R.xyz[,1]+rep(c((-1*round(RpN/2)):-1,1:(RpN-round(RpN/2)))*(p.d+1),each=N.n);R.xyz=wrap.box(box.dim,R.xyz);N.check=(sum(apply(rbind(N.xyz,R.xyz),1,r.calc,E.xyz=rbind(N.xyz,R.xyz))<=p.d & apply(rbind(N.xyz,R.xyz),1,r.calc,E.xyz=rbind(N.xyz,R.xyz))>0)>0)};rm(N.check)
E.xyz=array(runif(E.n*3,-0.5*box.dim,0.5*box.dim),c(E.n,3))
bind.counter=list('E'=rep(0,times=as.character(E.n)),'N'=rep(0,times=as.character(N.n)),'R'=rep(0,times=as.character(N.n*RpN)))
bind.counter=binding.update(E.xyz,N.xyz,R.xyz,bind.counter,kn1N,kn1R,tau,p.d)
bind.data=data.frame('N'=rep(NA,times=as.character(t/dt+1)),'R'=rep(NA,times=as.character(t/dt+1)),'Nd.m'=rep(NA,times=as.character(t/dt+1)),'Nd.sd'=rep(NA,times=as.character(t/dt+1)),'Nem'=rep(NA,times=as.character(t/dt+1)))
bind.data$N[1]=sum(bind.counter[['N']]>0)/N.n
bind.data$R[1]=sum(bind.counter[['R']]>0)/N.n/RpN
bind.data$Nd.m[1]=mean(apply(apply(N.xyz,1,r.calc,E.xyz=E.xyz[bind.counter[['E']]==0,]),2,skim.r,RpN=RpN))
bind.data$Nd.sd[1]=sd(apply(apply(N.xyz,1,r.calc,E.xyz=E.xyz[bind.counter[['E']]==0,]),2,skim.r,RpN=RpN))
bind.data$Nem[1]=em.calc(apply(N.xyz,1,r.calc,E.xyz=E.xyz[bind.counter[['E']]==0,]),bind.counter,N.n,p.d,box.dim)

# start simulation
COUNTER.dt=0
COUNTER.dtt=0
data.i=1
loop.i=1
est = Sys.time() - seconds(15)
if (is.null(start.file)==F){
  load(start.file)
  est=Sys.time()-ete
}
while (loop.i < Inf){
  COUNTER.dt=COUNTER.dt+1
  COUNTER.dtt=COUNTER.dtt+1
  timer=loop.i*tau
  E.xyz[bind.counter[['E']]==0,]=diffuse(E.xyz[bind.counter[['E']]==0,],D.step)
  E.xyz=wrap.box(box.dim,E.xyz)
  if(round(COUNTER.dt)==round(dt/tau)){
    COUNTER.dt=0
    data.i=data.i+1
    bind.data$N[data.i]=sum(bind.counter[['N']]>0)/N.n
    bind.data$R[data.i]=sum(bind.counter[['R']]>0)/N.n/RpN
    bind.data$Nd.m[data.i]=mean(apply(apply(N.xyz,1,r.calc,E.xyz=E.xyz[bind.counter[['E']]==0,]),2,skim.r,RpN=RpN))
    bind.data$Nd.sd[data.i]=sd(apply(apply(N.xyz,1,r.calc,E.xyz=E.xyz[bind.counter[['E']]==0,]),2,skim.r,RpN=RpN))
    bind.data$Nem[data.i]=em.calc(apply(N.xyz,1,r.calc,E.xyz=E.xyz[bind.counter[['E']]==0,]),bind.counter,N.n,p.d,box.dim)
    show(paste0('Progress = ',round(100*data.i/(t/dt+1),3),' %; ETA = ',ETA.calc(est,data.i/(t/dt+1))))
  }
  if(round(loop.i)==round(t/tau)){
    RESULTS=list('parameters'=parameters,'data'=bind.data,'final.system'=list("E.xyz"=E.xyz,"N.xyz"=N.xyz,"R.xyz"=R.xyz,'time'=timer,'states'=bind.counter))
    par(fig=c(0.01,0.5,0.5,1),mar = c(4.5,4.5,3,1))
    plot(NULL,NULL,main = 'Nucleosome Occupancy',xlim=c(0,t*1e3),ylim = c(0,1),xlab = 'Time (ms)',ylab = 'Fraction Bound',cex.main=2,cex.lab=1.7,cex.axis=1.7)
    lines(seq(0,t,dt)*1e3,bind.data$N,col = 'black',lwd=3)
    par(fig=c(0.51,1,0.5,1),mar = c(4.5,4.5,3,1),new=TRUE)
    plot(NULL,NULL,main = 'RNA Occupancy',xlim=c(0,t*1e3),ylim = c(0,1),xlab = 'Time (ms)',ylab = 'Fraction Bound',cex.main=2,cex.lab=1.7,cex.axis=1.7)
    lines(seq(0,t,dt)*1e3,bind.data$R,col = 'black',lwd=3)
    par(fig=c(0.01,0.5,0,0.5),mar = c(4.5,4.5,3,1),new=TRUE)
    plot(NULL,NULL,main = 'Nucleosome-Enzyme Proximity',xlim=c(0,t*1e3),ylim = c(0,max(bind.data$Nd.m+3*bind.data$Nd.sd)),xlab = 'Time (ms)',ylab = 'Nearest Enzymes (nm)',cex.main=1.8,cex.lab=1.7,cex.axis=1.7)
    lines(seq(0,t,dt)*1e3,bind.data$Nd.m,lty='solid',lwd=5)
    lines(seq(0,t,dt)*1e3,bind.data$Nd.m+bind.data$Nd.sd,lty='dashed',lwd=2)
    lines(seq(0,t,dt)*1e3,bind.data$Nd.m-bind.data$Nd.sd,lty='dashed',lwd=2)
    par(fig=c(0.51,1,0,0.5),mar = c(4.5,4.5,3,1),new=TRUE)
    plot(NULL,NULL,main = 'Effective Molarity',xlim=c(0,t*1e3),ylim = c(0,max(bind.data$Nem)),xlab = 'Time (ms)',ylab = 'Relative [E]',cex.main=1.8,cex.lab=1.7,cex.axis=1.7)
    abline(h=1,col='grey',lwd=2)
    abline(h=mean(bind.data$Nem),col='red',lty='dashed',lwd=3)
    lines(seq(0,t,dt)*1e3,bind.data$Nem,lty='solid',lwd=3)
    save.image(paste0(sim.name,'_','RESULTS.RData'))
    show(paste0('Finish Time: ',Sys.time()))
    stop('Simulation Completed!')
  }
  bind.counter=binding.update(E.xyz,N.xyz,R.xyz,bind.counter,kn1N,kn1R,tau,p.d)
  loop.i=loop.i+1
  if(round(COUNTER.dtt)==round(dtt/tau)){
    COUNTER.dtt=0
    if(file.exists(paste0(sim.name,'_','restart-old.RData'))){
      file.remove(paste0(sim.name,'_','restart-old.RData'))
    }
    if(file.exists(paste0(sim.name,'_','restart.RData'))){
      file.rename(paste0(sim.name,'_','restart.RData'),paste0(sim.name,'_','restart-old.RData'))
    }
    ete=Sys.time()-est
    save.image(paste0(sim.name,'_','restart.RData'))
  }
}

##### END SCRIPT









