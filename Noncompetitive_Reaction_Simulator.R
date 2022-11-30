##### SCRIPT INFORMATION

# Noncompetitive_Reaction_Simulator.R

# A script to simulate catalytic reactions under a model of noncompetitive binding between nucleosomes and RNA

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
Kd.N=200e-9 # nucleosome dissociation constant, M
Kd.R=400e-9 # decoy dissociation constant, M
kn1.N=2e-2 # nucleosome dissociation rate constant, 1/s
kn1.R=4e-2 # decoy dissociation rate constant, 1/s
kcat=0 # enzyme catalysis rate constant, 1/s
alpha=rev(c(2000/4^c(0:2),1)) # effective molarity adjustment parameter
beta=c(1/3^c(0:2),0) # RNA effect on catalysis adjustment parameter
delta1.N=1 # RNA effect on nucleosome on-rate
delta2.N=1e0 # RNA effect on nucleosome off-rate
delta1.R=1 # nucleosome effect on RNA on-rate
delta2.R=1e0 # nucleosome effect on RNA off-rate

# Initial values
N.T=50e-9 # total nucleosome, M
R.T=c(0,N.T*2^c(0:3)) # total RNA, M
E.T=Kd.N/3^c(-1:3) # total enzyme, M

### Extraneous variables

# Simulation parameters
time=60*60*1 # total reaction time, s
dt=0.25 # integration time-step, s

# Reaction parameters
pre.equilibrated=F # start from equilibrium binding state before initiating HMTase assay, T/F

# Saved data features
save.id='null' # identifier for save files
console.output=T # should console summary be output to file, T/F
plotting=F # should plots be produced, T/F
plot.type=2 # what should be plotted, 1=catalysis 2=occupancy 3=RNAlessOcupancy
legend.location='topright' # location of plot legend in legend() syntax
save.plots='yes' # decision to save plots, yes/no

##########################

##### BEGIN SCRIPT

### Semi-autonomous script functions

# Remaining model constant calculations
k1.N=kn1.N/Kd.N # prey association rate constant, 1/M/s
k1.R=kn1.R/Kd.R # decoy association rate constant, 1/M/s

### Autonomous script functions

# Save simulation parameters
parameters=mget(ls())

# Simulation results storage
SIM.Results=list()

# Record Equations
Equations=function(t,initials,constants){
  with(as.list(c(initials,constants)),{
    
    dE=kn1.N*(EN + ENm) + kn1.R*ER - E*(k1.N*(N + Nm) + k1.R*R)
    dN=kn1.N*(EN + delta2.N*ENR) - k1.N*N*(E + alpha*delta1.N*ER)
    dR=kn1.R*(ER + delta2.R*ENR + ENmR) - k1.R*R*(E + alpha*delta1.R*(EN + ENm))
    dNm=kn1.N*(ENm + delta2.N*ENmR) - k1.N*Nm*(E + alpha*delta1.N*ER)
    dEN=k1.N*E*N + delta2.R*kn1.R*ENR - EN*(kcat + alpha*delta1.R*k1.R*R + kn1.N)
    dER=k1.R*E*R + delta2.N*kn1.N*(ENR + ENmR) - ER*(kn1.R + alpha*delta1.N*k1.N*(N + Nm))
    dENm=kcat*EN + k1.N*E*Nm + delta2.R*kn1.R*ENmR - ENm*(kn1.N + alpha*delta1.R*k1.R*R)
    dENR=alpha*(delta1.R*k1.R*EN*R + delta1.N*k1.N*ER*N) - ENR*(delta2.R*kn1.R + beta*kcat + delta2.N*kn1.N)
    dENmR=beta*kcat*ENR + alpha*(delta1.R*k1.R*ENm*R + delta1.N*k1.N*ER*Nm) - ENmR*(delta2.R*kn1.R + delta2.N*kn1.N)
    
    list(c(dE,dN,dR,dNm,dEN,dER,dENm,dENR,dENmR))
  })
}

# Enzyme range
for (q in 1:length(E.T)){
  
  # Plot bins
  if (plotting==T){
    par(mfrow=c(length(beta),length(alpha)),mar=c(4.5,5,2,1),oma = c(0,0,5,0))
  }
  
  # Beta range
  for (p in 1:length(beta)){
    
    # Alpha range
    for (j in 1:length(alpha)){
    
      # Plot frame
      if (plotting==T & plot.type==1){
        plot(NULL,NULL,main = paste0('α=',round(alpha[j],3),'; β=',round(beta[p],3)),xlim = c(0,time/60),ylim = c(0,N.T),xlab = 'Time (min)',ylab = expression('[m'['T']*'] (M)'),cex.main=1.3,cex.axis=1.5,cex.lab=1.5)
      }
      if (plotting==T & plot.type==2){
        plot(NULL,NULL,main = paste0('α=',round(alpha[j],3),'; β=',round(beta[p],3)),xlim = c(0,time/60),ylim = c(0,1),xlab = 'Time (min)',ylab = expression('([EN]+[ENR])/[N'['T']*']'),cex.main=1.3,cex.axis=1.5,cex.lab=1.5)
      }
      if (plotting==T & plot.type==3){
        plot(NULL,NULL,main = paste0('α=',round(alpha[j],3),'; β=',round(beta[p],3)),xlim = c(0,time/60),ylim = c(0,1),xlab = 'Time (min)',ylab = expression('[EN]/[N'['T']*']'),cex.main=1.3,cex.axis=1.5,cex.lab=1.5)
      }
      
      # RNA range
      for (k in 1:length(R.T)){
      
        if (pre.equilibrated==F){
          # Record times
          t=seq(0,time,dt)
          
          # Record constants
          constants=c(k1.N=k1.N,k1.R=k1.R,kn1.N=kn1.N,kn1.R=kn1.R,kcat=kcat,alpha=alpha[j],beta=beta[p],delta1.N=delta1.N,delta2.N=delta2.N,delta1.R=delta1.R,delta2.R=delta2.R)
          
          # Record initial values
          initials=c(E=E.T[q],N=N.T,R=R.T[k],Nm=0,EN=0,ER=0,ENm=0,ENR=0,ENmR=0)
        }
        if (pre.equilibrated==T){
          t=seq(0,60*60*1,dt)
          constants=c(k1.N=k1.N,k1.R=k1.R,kn1.N=kn1.N,kn1.R=kn1.R,kcat=0,alpha=alpha[j],beta=beta[p],delta1.N=delta1.N,delta2.N=delta2.N,delta1.R=delta1.R,delta2.R=delta2.R)
          initials=c(E=E.T[q],N=N.T,R=R.T[k],Nm=0,EN=0,ER=0,ENm=0,ENR=0,ENmR=0)
          sim.sim=as.data.frame(ode(y=initials,times=t,func=Equations,parms=constants))
          
          # Record times
          t=seq(0,time,dt)
          
          # Record constants
          constants=c(k1.N=k1.N,k1.R=k1.R,kn1.N=kn1.N,kn1.R=kn1.R,kcat=kcat,alpha=alpha[j],beta=beta[p],delta1.N=delta1.N,delta2.N=delta2.N,delta1.R=delta1.R,delta2.R=delta2.R)
          
          # Record initial values
          initials=c(E=sim.sim$E[nrow(sim.sim)],N=sim.sim$N[nrow(sim.sim)],R=sim.sim$R[nrow(sim.sim)],Nm=sim.sim$Nm[nrow(sim.sim)],EN=sim.sim$EN[nrow(sim.sim)],ER=sim.sim$ER[nrow(sim.sim)],ENm=sim.sim$ENm[nrow(sim.sim)],ENR=sim.sim$ENR[nrow(sim.sim)],ENmR=sim.sim$ENmR[nrow(sim.sim)])
        }
        
        # Numerical evaluation
        sim.sim=as.data.frame(ode(y=initials,times=t,func=Equations,parms=constants))
        
        # Clean up results
        RXN.Results=list('t'=sim.sim$time,'E'=sim.sim$E,'N'=sim.sim$N,'R'=sim.sim$R,'Nm'=sim.sim$Nm,'EN'=sim.sim$EN,'ER'=sim.sim$ER,'ENm'=sim.sim$ENm,'ENR'=sim.sim$ENR,'ENmR'=sim.sim$ENmR,'mT'=sim.sim$ENmR+sim.sim$ENm+sim.sim$Nm,'rE'=E.T[q]/Kd.N,'rR'=R.T[k]/N.T,'alpha'=alpha[j],'beta'=beta[p],'On'=(sim.sim$ENmR+sim.sim$ENm+sim.sim$EN+sim.sim$ENR)/N.T,'Onx'=(sim.sim$ENm+sim.sim$EN)/N.T)
        SIM.Results[[paste(R.T[k],alpha[j],beta[p],E.T[q],sep = '_')]]=RXN.Results
        
        # Plot results
        if (plotting==T & plot.type==1){
          lines(RXN.Results[['t']]/60,RXN.Results[['mT']],col=k,lwd=2)
        }
        if (plotting==T & plot.type==2){
          lines(RXN.Results[['t']]/60,RXN.Results[['On']],col=k,lwd=2)
        }
        if (plotting==T & plot.type==3){
          lines(RXN.Results[['t']]/60,RXN.Results[['Onx']],col=k,lwd=2)
        }
        
      }
      
      # Add legend to plot results
      if (plotting==T){
        legend(legend.location,legend = paste0('RNA:Nuc = ',round(R.T/N.T,3)),col = 1:length(R.T),fill = 1:length(R.T),cex = 1)
        if (p==1 & j==1){
          # Panel labels
          title(paste0('E:KdN = ',round(E.T[q]/Kd.N,3)),line = 1.5,outer = TRUE,cex.main = 4)
        }
      }
      
    }
    
  }
  
  # Save plots
  if (plotting==T & save.plots=='yes'){
    rstudioapi::savePlotAsImage(paste0(save.id,'_plot-',q,'.png'),format ='png',width=1037,height=1000)
    dev.off()
  }
  
}

# Generate extrapolated data
ext.data=data.frame('Vo'=rep(NA,times=length(SIM.Results)),'On.eq'=rep(NA,times=length(SIM.Results)),'Onx.eq'=rep(NA,times=length(SIM.Results)),'rR'=rep(NA,times=length(SIM.Results)),'alpha'=rep(NA,times=length(SIM.Results)),'beta'=rep(NA,times=length(SIM.Results)),'rE'=rep(NA,times=length(SIM.Results)))
for (i in 1:length(SIM.Results)){
  ext.data$Vo[i]=SIM.Results[[i]][['mT']][SIM.Results[[i]][['t']]==(time/10)]/time*10
  ext.data$On.eq[i]=SIM.Results[[i]][['On']][SIM.Results[[i]][['t']]==time]
  ext.data$Onx.eq[i]=SIM.Results[[i]][['Onx']][SIM.Results[[i]][['t']]==time]
  ext.data$rR[i]=SIM.Results[[i]][['rR']]
  ext.data$alpha[i]=SIM.Results[[i]][['alpha']]
  ext.data$beta[i]=SIM.Results[[i]][['beta']]
  ext.data$rE[i]=SIM.Results[[i]][['rE']]
}

# Plot bins
par(mfrow=c(length(E.T),length(beta)),mar=c(4.5,5,2,1),oma = c(0,0,5,0))

# Enzyme range
for (q in 1:length(E.T)){
  
  # Beta range
  for (p in 1:length(beta)){
    
    # Plot frame
    if (plot.type==1){
      plot(NULL,NULL,main = paste0('β=',round(beta[p],3),'; E:KdN=',round(E.T[q]/Kd.N,3)),xlim = range(ext.data$rR),ylim = c(0,max(ext.data$Vo)),xlab = expression('[R'['T']*']'*'/[N'['T']*']'),ylab = expression(V[0]*' (M s'^-1*')'),cex.main=1.3,cex.axis=1.5,cex.lab=1.5)
    }
    if (plot.type==2){
      plot(NULL,NULL,main = paste0('β=',round(beta[p],3),'; E:KdN=',round(E.T[q]/Kd.N,3)),xlim = range(ext.data$rR),ylim = c(0,1),xlab = expression('[R'['T']*']'*'/[N'['T']*']'),ylab = expression('([EN]'[eq]*'+[ENR]'[eq]*')/[N'['T']*']'),cex.main=1.3,cex.axis=1.5,cex.lab=1.5)
    }
    if (plot.type==3){
      plot(NULL,NULL,main = paste0('β=',round(beta[p],3),'; E:KdN=',round(E.T[q]/Kd.N,3)),xlim = range(ext.data$rR),ylim = c(0,1),xlab = expression('[R'['T']*']'*'/[N'['T']*']'),ylab = expression('[EN]'[eq]*'/[N'['T']*']'),cex.main=1.3,cex.axis=1.5,cex.lab=1.5)
    }
    
    # Alpha range
    for (j in 1:length(alpha)){
      
      # Plot results
      if (plot.type==1){
        lines(ext.data$rR[ext.data$rE==(E.T[q]/Kd.N) & ext.data$beta==beta[p] & ext.data$alpha==alpha[j]],ext.data$Vo[ext.data$rE==(E.T[q]/Kd.N) & ext.data$beta==beta[p] & ext.data$alpha==alpha[j]],col=j,lwd=2)
      }
      if (plot.type==2){
        lines(ext.data$rR[ext.data$rE==(E.T[q]/Kd.N) & ext.data$beta==beta[p] & ext.data$alpha==alpha[j]],ext.data$On.eq[ext.data$rE==(E.T[q]/Kd.N) & ext.data$beta==beta[p] & ext.data$alpha==alpha[j]],col=j,lwd=2)
      }
      if (plot.type==3){
        lines(ext.data$rR[ext.data$rE==(E.T[q]/Kd.N) & ext.data$beta==beta[p] & ext.data$alpha==alpha[j]],ext.data$Onx.eq[ext.data$rE==(E.T[q]/Kd.N) & ext.data$beta==beta[p] & ext.data$alpha==alpha[j]],col=j,lwd=2)
      }
      
    }
      
    # Add legend to plot results
    legend(legend.location,legend = paste0('α = ',round(alpha,3)),col = 1:length(alpha),fill = 1:length(alpha),cex = 1)
    if (p==1 & q==1){
      # Panel labels
      title('All Data',line = 1.5,outer = TRUE,cex.main = 4)
    }
    
  }
    
}
  
# Save plots
if (save.plots=='yes'){
  rstudioapi::savePlotAsImage(paste0(save.id,'_plot-ALL.png'),format ='png',width=1037,height=1000)
  dev.off()
}

# Save simulation data
SIM.Results[['Parameters']]=parameters
save(SIM.Results,file=paste0(save.id,'.RData'))
if (console.output==T){
  sink(paste0(save.id,'.txt'))
  show(parameters)
  sink()
}

##### END SCRIPT









