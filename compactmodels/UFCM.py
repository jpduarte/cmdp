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
    
    #short-channel effect parameters for Vth roll-off, and DIBL, and Subthreshol-Swing
    self.vthrolloff = 0.0
    self.vthdibl    = 0.0
    self.SSrolloff  = 0.0
    self.SSdibl     = 0.0
    
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
      print "sweep"
   
    #short channel effect calculations
    deltaVth = flagdevtype*(self.vthrolloff+self.vthdibl* (Vd-Vs))
    nVtm = self.vt*(1.0+self.SSrolloff+self.SSdibl* (Vd-Vs))
 
    #quantum mechanical effects - bias dependence parameter calculations
    mx =  0.916 * self.MEL
    fieldnormalizationfactor  = nVtm*self.Cins/(self.Weff*self.ech)
    auxQMfact  = ((3.0/4.0)*3*self.HBAR*2.0*3.141516*self.q/(4*sqrt(2*mx)))**(2.0/3.0)
    QMFACTORCVfinal   = self.QMFACTORCV*auxQMfact*(fieldnormalizationfactor)**(2.0/3.0)*(1/(self.q*nVtm))       
       
    #source side evaluation  
    Vch = Vs
    qs = UFCMchargemodel.unified_charge_model(Vg-deltaVth,Vch,self.q,self.k,self.T,self.eo,self.eins,self.ech,self.Eg,self.Nc,self.Nv,nVtm,self.ni,self.phi_substrate,PHIG,self.alpha_MI,self.Cins,self.Ach,self.Weff,self.Nch,QMFACTORCVfinal)
    
    #drain side evaluation  
    Vch = Vd
    qd = UFCMchargemodel.unified_charge_model(Vg-deltaVth,Vch,self.q,self.k,self.T,self.eo,self.eins,self.ech,self.Eg,self.Nc,self.Nv,nVtm,self.ni,self.phi_substrate,PHIG,self.alpha_MI,self.Cins,self.Ach,self.Weff,self.Nch,QMFACTORCVfinal)
    
    #drain-source current model (normalized)
    ids0,mu = UFCMdraincurrentmodel.unified_normilized_ids(qs,qd,self.q,self.k,self.T,self.eo,self.eins,self.ech,self.Eg,self.Nc,self.Nv,nVtm,self.ni,self.phi_substrate,PHIG,self.alpha_MI,self.Cins,self.Ach,self.Weff,self.Nch,self.IDSMOD,self.DEVTYPE)

    #drain-source current in Ampere
    idsfactor = (nVtm**2*self.Cins)/self.Lg
    
    #source-drain sweep in case Vd<Vs
    if flagsweep ==1:
      ids0=-ids0
    
    #total current counting number of fins NFIN  
    Ids = flagdevtype*ids0*idsfactor*self.NFIN
    
    #gate charge  
    Qg = -((qs+qd)*0.5+(qs-qd)**2/(6*(2*(qs+qd)+1)))*self.Cins*nVtm*self.Lg*self.NFIN*flagdevtype
    
    """qm = qs
    Ft = -((qm+(-self.q*self.Nch*self.Ach)/(nVtm*self.Cins))*self.Cins*nVtm/( self.Weff * self.ech)) #1e-2 is to transform it to V/cm
    t = (Ft*self.Ach/self.Weff+nVtm*(1-exp((Ft/nVtm)*self.Ach/self.Weff)))/(Ft*(1-exp((Ft/nVtm)*self.Ach/self.Weff))) #1e3 to transfor to um
    c = -((qm)*self.Cins*nVtm/( self.Weff * self.q*(t))) #1e-6 to transform it to cm^-3
    mu,mudop,muc,muac,musr = mobilitymodels.mobility(self.DEVTYPE,Ft*1e-2,0,self.Nch*1e-6,self.T,t*1e3,c*1e-6) """   

    #attach values of variables to return
    variablesvalues = []
    for var in self.returnvar:
      exec 'variablesvalues.append('+var+')'
    #print variablesvalues
     
    return  variablesvalues,self.returnvar

