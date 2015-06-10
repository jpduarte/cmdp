#UFCM charge model it is based on; Unified FinFET Compact Model: Modelling Trapezoidal Triple-Gate FinFETs
#this includes quantum mechanical bias dependent effects
#Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

from numpy import exp,log,sqrt

def unified_charge_model(self,Vg,Vch,vt,phi_gate,QMF,nss) :
  vt = vt
  q = self.q
  k = self.k
  T = self.T
  eo = self.eo
  eins = self.eins
  ech = self.ech
  Eg = self.Eg
  Nc = self.Nc
  Nv = self.Nv
  ni = self.ni
  phi_substrate = self.phi_substrate
  alpha_MI = self.alpha_MI
  Cins = self.Cins
  Ach = self.Ach
  Weff = self.Weff
  Nch = self.Nch

  rc  = (2.0*Cins/(Weff**2.0*ech/Ach))
  qdep  = (-q*Nch*Ach)/((vt/nss)*Cins)
  vfb_n = (phi_gate - phi_substrate -Eg/2.0-vt*log(Nch/ni))/vt
  vth_fixed_factor_SI = vfb_n+log(Cins*vt/(q*ni**2.0*2.0*Ach/Nch))
  vth_fixed_factor_Sub = log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))+vth_fixed_factor_SI#updated for low doping, log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))
  vth_N_Sub = -qdep+vth_fixed_factor_Sub+QMF*((-qdep)**(2.0/3.0))
  vth_N_SI  = -qdep+vth_fixed_factor_SI+QMF*((-qdep)**(2.0/3.0))
  Vg_local_N  = (Vg-Vch)*vt**(-1.0)
  F           = -Vg_local_N+vth_N_SI
  Vov         = (Vg_local_N-vth_N_Sub)*0.5/nss
  if (Vov>60):
    qm      = -Vov*2.0*nss
    qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
    x0      = qtrc/(exp(qtrc)-qtrc-1.0)
    x1      = qtrc*x0
    
    f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
    f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
    f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
    
    f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
    f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
    f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))

    f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
    f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
    f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))    
  else:
    qm      = exp((Vg_local_N-vth_N_Sub)*0.5/nss)
    if(qm>1.0e-7):#original 1.0e-7
      qm      = 2.0*(1.0-sqrt(1.0+(log(1.0+qm))**2.0))
      qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
      x0      = qtrc/(exp(qtrc)-qtrc-1.0)
      x1      = qtrc*x0
      
      f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
      f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
      f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
        
      f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
      f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
      f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0)) 
        
      f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
      f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
      f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))    
    else:
      qm      = -qm**2.0
  return qm  
  

