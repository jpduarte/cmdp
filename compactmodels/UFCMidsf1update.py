#derivate for drain current update, this calculate the derivate (f1) used to update the drain-source current in the presence of source/drain resistances
#Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

from numpy import exp,log,sqrt

from scipy.integrate import quad
from scipy.integrate import quadrature
from scipy.special.orthogonal import p_roots

import mobilitymodels

import UFCMvdsat
import UFCMchargemodel

def f1(self,qs,qd,vt,phi_gate,Rs,Rd):

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
  IDSMOD = self.IDSMOD
  DEVTYPE =self.DEVTYPE
  Lg = self.Lg
  
  qm = qs
  Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
  t = (Ft*Ach/Weff+vt*(1.0-exp((Ft/vt)*Ach/Weff)))/(Ft*(1.0-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
  c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
  mulow1 = mobilitymodels.mobility(self,DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))

  qm = qd
  Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
  t = (Ft*Ach/Weff+vt*(1.0-exp((Ft/vt)*Ach/Weff)))/(Ft*(1.0-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
  c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
  mulow2 = mobilitymodels.mobility(self,DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))

  #calculation of derivative
  return 1-1.0/Lg*(mulow1*qs*Cins*vt*Rs+mulow2*qd*Cins*vt*Rd)
