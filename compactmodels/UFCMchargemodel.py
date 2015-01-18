#UFCM charge

from numpy import exp,log,sqrt

def unified_charge_model(Vg,Vch,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch) :
  rc  = (2*Cins/(Weff**2*ech/Ach))
  qdep  = (-q*Nch*Ach)/(vt*Cins)
  vfb_n = (phi_gate - phi_substrate -Eg/2-vt*log(Nch/ni))/vt
  vth_fixed_factor_SI = vfb_n+log(Cins*vt/(q*ni**2*2*Ach/Nch))
  vth_fixed_factor_Sub = log((qdep*rc)**2/(exp(qdep*rc)-qdep*rc-1))+vth_fixed_factor_SI#updated for low doping, log((qdep*rc)**2/(exp(qdep*rc)-qdep*rc-1))
  vth_N_Sub = -qdep+vth_fixed_factor_Sub
  vth_N_SI  = -qdep+vth_fixed_factor_SI
  Vg_local_N  = (Vg-Vch)*vt**(-1)
  F           = -Vg_local_N+vth_N_SI
  Vov         = (Vg_local_N-vth_N_Sub)*0.5
  if (Vov>60):
    qm      = -Vov*2
    qtrc    = (qm*alpha_MI**(-1)+qdep)*rc
    x0      = qtrc*(exp(qtrc)-qtrc-1)**(-1)
    x1      = qtrc*x0
    f0      = F-qm+log(-qm)+log(x1)
    f1      = -1+qm**(-1)+(2*qtrc**(-1)-x0-1)*rc
    f2      = -(qm**2)**(-1)
    qm      = qm-(f0*f1**(-1))*(1+(f0*f2)*(2*f1**2)**(-1))
  else:
    qm      = exp((Vg_local_N-vth_N_Sub)*0.5)
    if(qm>1e-3):#original 1e-7
      qm      = 2*(1-sqrt(1+(log(1+qm))**2))
      qtrc    = (qm*alpha_MI**(-1)+qdep)*rc
      x0      = qtrc*(exp(qtrc)-qtrc-1)**(-1)
      x1      = qtrc*x0
      f0      = F-qm+log(-qm)+log(x1)
      f1      = -1+qm**(-1)+(2*qtrc**(-1)-x0-1)*rc
      f2      = -(qm**2)**(-1)
      qm      = qm-(f0*f1**(-1))*(1+(f0*f2)*(2*f1**2)**(-1))
    else:
      qm      = -qm**2
  return qm  
