#UFCM: Unified FinFET Compact Model
#author: Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

import UFCMchargemodel
import UFCMdraincurrentmodel
import mobilitymodels

from numpy import sqrt, exp

class compactmodel:
  def __init__(self,name):
    self.name = name
    self.version = 'v1'
    self.nodenames = ['vd','vg','vs','vb']
    
    #TODO: add units to each quantity   
    self.q       = 1.6e-19 #C
    self.k       = 1.380650e-23 
    self.T       = 300 #K
    self.HBAR    = 1.05457e-34 # untis are Joule-sec
    self.MEL     = 9.11e-31 # kg
    self.eo      = 8.854000e-12 # F / m
    self.eins    = 3.9*self.eo # F / m 
    self.ech     = 11.7*self.eo # F / m 
    self.Eg 	   = 1.169640 #eV
    self.Nc      = 2.890000e+25 #1/m^3
    self.Nv 	   = 3.140000e+25 #1/m^3 
    self.u0      = 40e-3 # m^2V^-1s^-1
    #device dimensions and parameters
    self.HFIN             = 20e-9 #m
    self.TFIN             = 10e-9 #m   
    self.TFIN_TOP         = 10e-9 #m
    self.TFIN_BASE        = 10e-9 #m  
    self.TFIN_TOP_L       = 5e-9 #m
    self.TFIN_BASE_L      = 5e-9 #m   
    self.TFIN_TOP_R       = 5e-9/2 #m
    self.TFIN_BASE_R      = 5e-9 #m          
    self.tins             = 1e-9 #m
    self.Nch              = 1.0e21 #1/m^3
    self.phi_substrate    = 4.05 #eV
    self.PHIG 	          = 4.5 #eV
    self.alpha_MI         = 4
    self.Lg               = 1e-6
    self.returnvar        = ['Ids']
    self.QMFACTORCV       = 0
    self.NFIN             = 1
    #UFCM parameters
    self.GEOMOD           = 0
    self.IDSMOD           = 0
    self.DEVTYPE          = 0
    
    self.MODmudop         = 1
    self.MODmuac          = 1
    self.MODmusr          = 1
    self.MODmurcs         = 1
    self.MODmuc           = 1
    
    #short-channel effect parameters for Vth roll-off, and DIBL, and Subthreshol-Swing
    self.vthrolloff = 0.0
    self.vthdibl    = 0.0
    self.SSrolloff  = 0.0
    self.SSdibl     = 0.0
    
    #current saturation
    self.betavsat = 2.0 
    self.vsat = 8.37e6*1e-2# m/s   
    
    self.ul = 470.5 #cm^2/Vs
    self.xi = 2.2  
    
    #Masetti Model: doping dependent
    self.masetti_umin1 = 44.9
    self.masetti_umin2 = 0.0
    self.masetti_u1    = 29.0
    self.masetti_Pc    = 9.23e16
    self.masetti_Cr    = 2.23e17
    self.masetti_Cs    = 6.1e20
    self.masetti_alpha = 0.719
    self.masetti_beta  = 2.0    
    
    self.lombardi_B = 9.925e6
    self.lombardi_C = 15e3#
    self.lombardi_N0 = 1.0
    self.lombardi_N2 = 1.0
    self.lombardi_lambdam = 0.0317
    self.lombardi_k = 1.0
    self.lombardi_delta = 0.7e14#
    self.lombardi_eta = 0.1e30#
    self.lombardi_lcrit = 1e-6 #cm
    self.lombardi_Astar = 2.0
    self.lombardi_Fref = 1    
    
    self.coulomb_D1inv = 135.0 #cm^2/Vs
    self.coulomb_alpha1inv = 0.0

    self.coulomb_cother = 0.0
    self.coulomb_v0inv = 1.5

    self.coulomb_v1inv = 2.0
    self.coulomb_D2inv =  40.0 #cm^2/Vs
    self.coulomb_alpha2inv = 0.0
    self.coulomb_v2inv = 0.5
   
    self.coulomb_t1 = 0.0003 #um
    self.coulomb_tcoulomb = 0.0 #um
    
    self.coulomb_S = 0.3042 #K(cm/V)^(2/3)
    self.coulomb_St = 0.0 #K(um)^2
    self.coulomb_t0 = 0.0005 #um
    
    self.coulomb_m_over_m0 = 1
    
    self.coulomb_mumax = 470.5
    self.coulomb_mumin = 44.9
    self.coulomb_alpha = 0.719
    
    self.coulomb_Nref = 2.23e17 #cm^-3
    self.coulomb_p = 4.0
    
    self.rcoulomb_murcs0  =  149.0 # cm^2/Vs
    self.rcoulomb_gamma1  = -0.23187
    self.rcoulomb_gamma2  = 2.1
    self.rcoulomb_gamma3  = 0.4
    self.rcoulomb_gamma4  = 0.05
    self.rcoulomb_gamma5  = 1.0
    self.rcoulomb_s       = 0.1
    self.rcoulomb_c0      = 3.0e16 #cm^-3
    self.rcoulomb_dcrit   = 0.0 # cm
    self.rcoulomb_lcrit   = 1e-6 # cm
    self.rcoulomb_lcrithk = 1e-6 # cm
    self.rcoulomb_xi      = 1.3042e7 # V-1cm-1
    
  def analog(self,*args):
    Vdi,Vgi,Vsi,Vbi = args
    """this function returns the drain current or other variables for given bias"""
    
    ###################################################
    self.vt      = self.k*self.T/self.q #V 
    self.ni      = sqrt(self.Nc*self.Nv)*exp(-self.Eg/(2*self.vt)) #1/m^3
    #geometrie Unified FinFET model parameter calculations
    if (self.GEOMOD==0): #double-gate, rec
      self.Cins   = (2.0*self.HFIN)*self.eins/self.tins
      self.Ach    = self.TFIN*self.HFIN
      self.Weff   = 2.0*self.HFIN
    if (self.GEOMOD==1): #triple-gate
      self.Cins   = (2.0*self.HFIN+self.TFIN)*self.eins/self.tins
      self.Ach    = self.TFIN*self.HFIN
      self.Weff   = (2.0*self.HFIN+self.TFIN)    
    if (self.GEOMOD==5): #trapezoidal
      self.Weff   = sqrt((self.TFIN_TOP_L-self.TFIN_BASE_L)**2+self.HFIN**2)+sqrt((self.TFIN_TOP_R-self.TFIN_BASE_R)**2+self.HFIN**2)+self.TFIN_TOP_R+self.TFIN_TOP_L
      self.Cins   = self.Weff*self.eins/self.tins
      self.Ach    = self.HFIN*(self.TFIN_TOP_R+self.TFIN_TOP_L + self.TFIN_BASE_L+self.TFIN_BASE_R)/2
 
    ##################################################### 
    #this check if this is PMOS or NMOS    
    if (self.DEVTYPE==-1): #PMOS
      flagdevtype = -1.0
      PHIG = -self.PHIG+2*self.phi_substrate+self.Eg
    else: #NMOS
      flagdevtype = 1.0
      PHIG = self.PHIG

    Vg = flagdevtype*(Vgi-Vbi)
    Vs = flagdevtype*(Vsi-Vbi)
    Vd = flagdevtype*(Vdi-Vbi)
    Vb = 0.0
    flagsweep = 0.0  
    #sweep voltage reference for drain-source
    if (Vd<Vs):
      Vs = Vd
      Vd = flagdevtype*(Vsi-Vbi)
      flagsweep = 1.0
      #print "sweep"
   
    #short channel effect calculations
    deltaVth = flagdevtype*(self.vthrolloff+self.vthdibl* (Vd-Vs))
    nVtm = self.vt*(1.0+self.SSrolloff+self.SSdibl* (Vd-Vs))
 
    #quantum mechanical effects - bias dependence parameter calculations
    mx =  0.916 * self.MEL
    fieldnormalizationfactor  = nVtm*self.Cins/(self.Weff*self.ech)
    auxQMfact  = ((3.0/4.0)*3*self.HBAR*2.0*3.141516*self.q/(4*sqrt(2*mx)))**(2.0/3.0)
    QMf   = self.QMFACTORCV*auxQMfact*(fieldnormalizationfactor)**(2.0/3.0)*(1/(self.q*nVtm))       
       
    #source side evaluation for charge  
    Vch = Vs
    qs = UFCMchargemodel.unified_charge_model(self,Vg-deltaVth,Vch,nVtm,PHIG,QMf)
   
    #drain-source current model (normalized)
    ids0,mu,vdsat,qd = UFCMdraincurrentmodel.unified_normilized_ids(self,qs,nVtm,PHIG,Vd,Vs,Vg,QMf,deltaVth)

    #drain-source current in Ampere
    idsfactor = (nVtm**2*self.Cins)/self.Lg
    
    #source-drain sweep in case Vd<Vs
    if flagsweep ==1:
      ids0=-ids0
    
    #total current counting number of fins NFIN  
    Ids = flagdevtype*ids0*idsfactor*self.NFIN
    
    #gate charge  
    Qg = -((qs+qd)*0.5+(qs-qd)**2/(6*(2*(qs+qd)+1)))*self.Cins*nVtm*self.Lg*self.NFIN*flagdevtype
    
    #attach values of variables to return
    variablesvalues = []
    for var in self.returnvar:
      exec 'variablesvalues.append('+var+')'
    #print variablesvalues
     
    return  variablesvalues,self.returnvar

