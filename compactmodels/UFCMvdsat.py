#This function is to calculate the saturation drain voltage for a given saturation drain charge
#Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

from numpy import exp,log,sqrt

def vdsat(self,Vg,qm,vt,phi_gate,QMF,nss) :

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
  qdep  = (-q*Nch*Ach)/(vt*Cins)
  vfb_n = (phi_gate - phi_substrate -Eg/2.0-vt*log(Nch/ni))/vt
  vth_fixed_factor_SI = vfb_n+log(Cins*vt/(q*ni**2.0*2.0*Ach/Nch))
  vth_fixed_factor_Sub = log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))+vth_fixed_factor_SI#updated for low doping, log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))
  vth_N_Sub = -qdep+vth_fixed_factor_Sub+QMF*((-qdep)**(2.0/3.0))
  vth_N_SI  = -qdep+vth_fixed_factor_SI+QMF*((-qdep)**(2.0/3.0))
  
  qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
  x0      = qtrc/(exp(qtrc)-qtrc-1.0)
  x1      = qtrc*x0  
  Vch  =Vg-vt*(vth_N_SI-qm+nss*log(-qm)+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0)))
  return Vch  
