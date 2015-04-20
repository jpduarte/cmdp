#UFCM charge

from numpy import exp,log,sqrt

def unified_charge_model(Vg,Vch,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch,QMF) :
  rc  = (2.0*Cins/(Weff**2.0*ech/Ach))
  qdep  = (-q*Nch*Ach)/(vt*Cins)
  vfb_n = (phi_gate - phi_substrate -Eg/2.0-vt*log(Nch/ni))/vt
  vth_fixed_factor_SI = vfb_n+log(Cins*vt/(q*ni**2.0*2.0*Ach/Nch))
  vth_fixed_factor_Sub = log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))+vth_fixed_factor_SI#updated for low doping, log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))
  vth_N_Sub = -qdep+vth_fixed_factor_Sub+QMF*((-qdep)**(2.0/3.0))
  vth_N_SI  = -qdep+vth_fixed_factor_SI+QMF*((-qdep)**(2.0/3.0))
  Vg_local_N  = (Vg-Vch)*vt**(-1.0)
  F           = -Vg_local_N+vth_N_SI
  Vov         = (Vg_local_N-vth_N_Sub)*0.5
  if (Vov>60):
    qm      = -Vov*2.0
    qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
    x0      = qtrc/(exp(qtrc)-qtrc-1.0)
    x1      = qtrc*x0
    
    f0      = F-qm+log(-qm)+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
    f1      = -1.0+qm**(-1.0)+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
    f2      = -(qm**2.0)**(-1.0)-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
    
    f0      = F-qm+log(-qm)+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
    f1      = -1.0+qm**(-1.0)+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
    f2      = -(qm**2.0)**(-1.0)-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))  

    f0      = F-qm+log(-qm)+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
    f1      = -1.0+qm**(-1.0)+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
    f2      = -(qm**2.0)**(-1.0)-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))       
  else:
    qm      = exp((Vg_local_N-vth_N_Sub)*0.5)
    if(qm>1.0e-7):#original 1.0e-7
      qm      = 2.0*(1.0-sqrt(1.0+(log(1.0+qm))**2.0))
      qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
      x0      = qtrc/(exp(qtrc)-qtrc-1.0)
      x1      = qtrc*x0
      
      f0      = F-qm+log(-qm)+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
      f1      = -1.0+qm**(-1.0)+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
      f2      = -(qm**2.0)**(-1.0)-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
      
      f0      = F-qm+log(-qm)+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
      f1      = -1.0+qm**(-1.0)+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
      f2      = -(qm**2.0)**(-1.0)-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))   
      
      f0      = F-qm+log(-qm)+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
      f1      = -1.0+qm**(-1.0)+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
      f2      = -(qm**2.0)**(-1.0)-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))          
    else:
      qm      = -qm**2.0
  return qm  
