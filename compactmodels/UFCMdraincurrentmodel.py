#UFCM drain current model, it has two modes: 
#IDSMOD=0: uses effective mobility and does not account for drain-source current saturation 
#IDSMOD=1: it include mobility dependance on the channel charge and current saturation effects
#Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

from numpy import exp,log,sqrt

from scipy.integrate import quad
from scipy.integrate import quadrature
from scipy.special.orthogonal import p_roots

import mobilitymodels

import UFCMvdsat
import UFCMchargemodel

def unified_normilized_ids(self,qms,vt,phi_gate,Vd,Vs,Vg,QMFACTORCVfinal,deltaVth,SS) :
  
  #parameter definitions
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
  
  Vdseff = 0
  if (IDSMOD==0):
    #normalized drain current model ref: Unified FinFET Compact Model: Modelling Trapezoidal Triple-Gate FinFETs 
    Vch = Vd
    qmd = UFCMchargemodel.unified_charge_model(self,Vg-deltaVth,Vch,vt,phi_gate,QMFACTORCVfinal,SS)
  
    rc  = (2*Cins/(Weff**2*ech/Ach))#TODO: check factor of 2 in rc?
    qdep  = (-q*Nch*Ach)/(vt*Cins)
    qh = 1/rc-qdep
    qm = qms
    idss = qm**2/2-2*qm-qh*log(1-qm/qh)
    qm = qmd
    idsd = qm**2/2-2*qm-qh*log(1-qm/qh)
    #effective mobility calculation, averaging source/drain charges 
    qm = (qms+qmd)/2
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1-exp((Ft/vt)*Ach/Weff)))/(Ft*(1-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mu = mobilitymodels.mobility(self,DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))   
    #drain-source current calculation
    ids =  mu*(idss-idsd)
  
  if (IDSMOD==1):
    #include mobility dependance accross the channel using Gauss Quadrature effects
    
    #Gauss Quadrature coefficients
    x1 = -0.57735027
    x2 = 0.57735027
    w1 = 1.0
    w2 = 1.0    
    
    #guess for q1 charge so mobility can be approximated
    q1 = (0.0-qms)*(x1+1)/2.0 + qms
    
    qm = q1
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1.0-exp((Ft/vt)*Ach/Weff)))/(Ft*(1.0-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mulow1 = mobilitymodels.mobility(self,DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    
    mulow2 =mulow1
    ####################################################################################################
    # calculation considering current saturation effects
    #################################################################################################### 
    beta = self.betavsat # 2.0 
    vsat = self.vsat # m/s 
    
    #drain saturation charge calculation: qd,sat            
    qdsat = (2.0*Lg*vsat + mulow1*qms*vt*x1 + mulow1*vt + mulow2*qms*vt*x2 + mulow2*vt - sqrt(4.0*Lg**2*vsat**2 + 4.0*Lg*mulow1*qms*vsat*vt*x1 + 4.0*Lg*mulow1*vsat*vt + 4.0*Lg*mulow2*qms*vsat*vt*x2 + 4.0*Lg*mulow2*vsat*vt + mulow1**2*qms**2*vt**2 - 2.0*mulow1**2*qms*vt**2 + mulow1**2*vt**2 + 2.0*mulow1*mulow2*qms**2*vt**2 - 4.0*mulow1*mulow2*qms*vt**2 + 2.0*mulow1*mulow2*vt**2 + mulow2**2*qms**2*vt**2 - 2.0*mulow2**2*qms*vt**2 + mulow2**2*vt**2))/(mulow1*vt*x1 + mulow1*vt + mulow2*vt*x2 + mulow2*vt)
    
    #Vdsat calculation 
    Vdsat = UFCMvdsat.vdsat(self,Vg-deltaVth,qdsat,vt,phi_gate,QMFACTORCVfinal,SS)-Vs 

    #Vdseff calculation
    Vds = Vd-Vs
    if (Vds>1e-8):
      Vdseff = 1.0/((1.0/Vdsat)**beta+(1.0/Vds)**beta)**(1.0/beta)
    else:
      Vdseff = Vds
    
    #calculation of mobile charge at drain-side, using Vdseff
    #Note that Vs is added to account for right voltage reference  
    qmd = UFCMchargemodel.unified_charge_model(self,Vg-deltaVth,Vdseff+Vs,vt,phi_gate,QMFACTORCVfinal,SS)
    
    #current calculation using Gauss Quadrature technique
    q1 = (qmd-qms)*(x1+1)/2.0 + qms
    q2 = (qmd-qms)*(x2+1)/2.0 + qms
    
    qm = q1
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1.0-exp((Ft/vt)*Ach/Weff)))/(Ft*(1.0-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mulow1 = mobilitymodels.mobility(self,DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    ids1 = -qm*(1.0-1.0/qm)*mulow1*w1
    
    qm = q2
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1.0-exp((Ft/vt)*Ach/Weff)))/(Ft*(1.0-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mulow2 = mobilitymodels.mobility(self,DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    ids2 = -qm*(1.0-1.0/qm)*mulow2*w2   
   
    ids = 0.5*(qmd-qms)*(ids1+ids2)    

    mu  = mulow1
    
  return ids,mu,Vdseff,qmd
   
