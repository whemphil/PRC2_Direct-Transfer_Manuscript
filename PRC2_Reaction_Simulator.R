##### SCRIPT INFORMATION

# PRC2_Reaction_Simulator.R

# A script to simulate HMTase reactions under a "hand-off" model for PRC2

# Associated with manuscript: "PRC2 'Hand-Off' From G-quadruplex RNA to dsDNA: Implications For RNA-Binding Chromatin Modifiers"

##########################

##### DEVELOPER NOTES

# The manuscript this code is associated with only utilizes variation in the 'Environment' and 'Key variables' sections of SCRIPT USER INPUT

##########################

##### SCRIPT USER INPUT

### Environment

rm(list = ls())
setwd('~/path/to/output/')

# Accessory packages
library(deSolve)

### Key variables

# Model constants
Kd.N=5.1e-9 # nucleosome dissociation constant, M
Kd.R=2.3e-9 # decoy dissociation constant, M
kn1.N=5.3e-4 # nucleosome dissociation rate constant, 1/s
kn1.R=8.3e-4 # decoy dissociation rate constant, 1/s
ktheta.N=100 # nucleosome displacement-transfer rate constant, 1/M/s
ktheta.R=67 # decoy displacement-transfer rate constant, 1/M/s
ktheta.NN=150 # nucleosome-nucleosome displacement-transfer rate constant, 1/M/s
kcat=1e0 # enzyme catalysis rate constant, 1/s
alpha=rev(c(2000/4^c(0:2),1)) # effective molarity adjustment parameter

# Initial values
N.T=5e-9 # total nucleosome, M
R.T=c(0,N.T*2^c(0:3)) # total RNA, M
E.T=Kd.N*2^c(-3:1) # total enzyme, M

### Extraneous variables

# Simulation parameters
time=60*60*1 # total reaction time, s
dt=1e0 # integration time-step, s

# Reaction parameters
pre.equilibrated=F # start from equilibrium binding state before initiating HMTase assay, T/F

# Saved data features
save.id='null' # identifier for save files
console.output=T # should console summary be output to file, T/F
plotting=T # should plots be produced, T/F
legend.location='topleft' # location of plot legend in legend() syntax
save.pdf='no' # decision to save pdf of graphs, yes/no

##########################

##### BEGIN SCRIPT

### Semi-autonomous script functions

# Remaining model constant calculations
k1.N=kn1.N/Kd.N # prey association rate constant, 1/M/s
k1.R=kn1.R/Kd.R # decoy association rate constant, 1/M/s

### Autonomous script functions

# Save simulation parameters
parameters=mget(ls())

# Plot bins
if (save.pdf=='yes'){
  pdf(file = paste(save.id,'.pdf',sep=''))
}
if (length(alpha)>1 & plotting==T){
  par(mfrow=c(length(E.T),length(alpha)),mar=c(4.5,5,2,1))
}

# Simulation results storage
SIM.Results=list()

# Record Equations
Equations=function(t,initials,constants){
  with(as.list(c(initials,constants)),{
    
    dE=kn1.N*(EN+ENm)+kn1.R*ER-E*(k1.N*(N+Nm)+k1.R*R)
    dN=EN*(kn1.N+ktheta.NN*Nm+alpha*ktheta.R*R)-N*(k1.N*E+ktheta.NN*ENm+alpha*ktheta.N*ER)
    dR=ER*(kn1.R+alpha*ktheta.N*(N+Nm))-R*(k1.R*E+alpha*ktheta.R*(EN+ENm))
    dNm=ENm*(ktheta.NN*N+alpha*ktheta.R*R+kn1.N)-Nm*(ktheta.NN*EN+k1.N*E+alpha*ktheta.N*ER)
    dEN=N*(k1.N*E+ktheta.NN*ENm+alpha*ktheta.N*ER)-EN*(kcat+kn1.N+ktheta.NN*Nm+alpha*ktheta.R*R)
    dER=R*(k1.R*E+alpha*ktheta.R*(EN+ENm))-ER*(alpha*ktheta.N*(N+Nm)+kn1.R)
    dENm=kcat*EN+Nm*(ktheta.NN*EN+k1.N*E+alpha*ktheta.N*ER)-ENm*(ktheta.NN*N+alpha*ktheta.R*R+kn1.N)
    
    list(c(dE,dN,dR,dNm,dEN,dER,dENm))
  })
}

# Enzyme range
for (p in 1:length(E.T)){
  
  # Alpha range
  for (j in 1:length(alpha)){
  
    # Plot frame
    if (plotting==T){
      plot(NULL,NULL,main = paste0('α = ',alpha[j],'; E:KdN(R) = ',signif(E.T[p]/Kd.N,2),'(',signif(E.T[p]/Kd.R,2),')'),xlim = c(0,time/60),ylim = c(0,N.T),xlab = 'Time (min)',ylab = '[H3K27me3] (M)',cex.main=1.3,cex.axis=1.5,cex.lab=1.5)
    }
    
    # RNA range
    for (k in 1:length(R.T)){
    
      if (pre.equilibrated==F){
        # Record times
        t=seq(0,time,dt)
        
        # Record constants
        constants=c(k1.N=k1.N,k1.R=k1.R,kn1.N=kn1.N,kn1.R=kn1.R,ktheta.N=ktheta.N,ktheta.R=ktheta.R,ktheta.NN=ktheta.NN,kcat=kcat,alpha=alpha[j])
        
        # Record initial values
        initials=c(E=E.T[p],N=N.T,R=R.T[k],Nm=0,EN=0,ER=0,ENm=0)
      }
      if (pre.equilibrated==T){
        t=seq(0,60*60*1,dt)
        constants=c(k1.N=k1.N,k1.R=k1.R,kn1.N=kn1.N,kn1.R=kn1.R,ktheta.N=ktheta.N,ktheta.R=ktheta.R,ktheta.NN=ktheta.NN,kcat=0,alpha=alpha[j])
        initials=c(E=E.T[p],N=N.T,R=R.T[k],Nm=0,EN=0,ER=0,ENm=0)
        sim.sim=as.data.frame(ode(y=initials,times=t,func=Equations,parms=constants))
        
        # Record times
        t=seq(0,time,dt)
        
        # Record constants
        constants=c(k1.N=k1.N,k1.R=k1.R,kn1.N=kn1.N,kn1.R=kn1.R,ktheta.N=ktheta.N,ktheta.R=ktheta.R,ktheta.NN=ktheta.NN,kcat=kcat,alpha=alpha[j])
        
        # Record initial values
        initials=c(E=sim.sim$E[nrow(sim.sim)],N=sim.sim$N[nrow(sim.sim)],R=sim.sim$R[nrow(sim.sim)],Nm=sim.sim$Nm[nrow(sim.sim)],EN=sim.sim$EN[nrow(sim.sim)],ER=sim.sim$ER[nrow(sim.sim)],ENm=sim.sim$ENm[nrow(sim.sim)])
      }
      
      # Numerical evaluation
      sim.sim=as.data.frame(ode(y=initials,times=t,func=Equations,parms=constants))
      
      # Clean up results
      RXN.Results=list('t'=sim.sim$time,'E'=sim.sim$E,'N'=sim.sim$N,'R'=sim.sim$R,'Nm'=sim.sim$Nm,'EN'=sim.sim$EN,'ER'=sim.sim$ER,'ENm'=sim.sim$ENm,'mT'=sim.sim$ENm+sim.sim$Nm)
      SIM.Results[[paste(R.T[k],alpha[j],E.T[p],sep = '_')]]=RXN.Results
      
      # Plot results
      if (plotting==T){
        lines(RXN.Results[['t']]/60,RXN.Results[['mT']],col=k,lwd=2)
      }
      
    }
    legend(legend.location,legend = paste0('RNA:Nuc = ',R.T/N.T),col = 1:length(R.T),fill = 1:length(R.T),cex = 1)
  
  }
  
}

# Save simulation data
SIM.Results[['Parameters']]=parameters
save(SIM.Results,file=paste0(save.id,'.RData'))
if (save.pdf=='yes'){
  dev.off()
}
if (console.output==T){
  sink(paste0(save.id,'.txt'))
  show(parameters)
  sink()
}

##### END SCRIPT









