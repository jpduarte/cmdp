#UFCM charge

from numpy import exp,log,sqrt

def vdsat(Vg,qm,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch,QMF) :
  rc  = (2.0*Cins/(Weff**2.0*ech/Ach))
  qdep  = (-q*Nch*Ach)/(vt*Cins)
  vfb_n = (phi_gate - phi_substrate -Eg/2.0-vt*log(Nch/ni))/vt
  vth_fixed_factor_SI = vfb_n+log(Cins*vt/(q*ni**2.0*2.0*Ach/Nch))
  vth_fixed_factor_Sub = log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))+vth_fixed_factor_SI#updated for low doping, log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))
  vth_N_Sub = -qdep+vth_fixed_factor_Sub+QMF*((-qdep)**(2.0/3.0))
  vth_N_SI  = -qdep+vth_fixed_factor_SI+QMF*((-qdep)**(2.0/3.0))
  
  qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
  x0      = qtrc/(exp(qtrc)-qtrc-1.0)
  x1      = qtrc*x0  
  Vch  =Vg-vt*(vth_N_SI-qm+log(-qm)+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0)))
  return Vch  
