#UFCM: Unified FinFET Compact Model
#author: Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

import UFCMchargemodel
import UFCMdraincurrentmodel

from numpy import sqrt, exp

class UFCM:
  def __init__(self,name):
    self.name = name
    self.version = 'v1'
    self.nodenames = ['vd','vg','vs','vb']
    #TODO: add units to each quantity   
    self.q       = 1.6e-19 #C
    self.vt      = 0.0259 #V
    self.k       = 1.380650e-23 
    self.T       = 300 #K
    self.eo      = 8.854000e-12 # F / m
    self.eins    = 3.453060e-11 # F / m 
    self.ech     = 1.035918e-10 # F / m 
    self.Eg 	    = 1.169640 #eV
    self.Nc      = 2.890000e+25 #1/m^3
    self.Nv 	    = 3.140000e+25 #1/m^3 
    self.ni      = 5e16 #sqrt(self.Nc*self.Nv)*exp(-self.Eg/(2*self.vt)) #1.1e16 #1/m^3
    self.u0     = 40e-3 # m^2V^-1s^-1
    #device dimensions and parameters
    self.HFIN             = 42e-9 #m   
    self.TFIN_TOP         = 7.6e-9 #m
    self.TFIN_BASE        = 2*7.6e-9 #m  
    self.TFIN_TOP_L       = 7.6e-9/2 #m
    self.TFIN_BASE_L      = 7.6e-9 #m   
    self.TFIN_TOP_R       = 7.6e-9/2 #m
    self.TFIN_BASE_R      = 7.6e-9 #m          
    self.tins             = 0.8e-9 #m
    self.Nch              = 6.0e24 #1/m^3
    self.phi_substrate    = 4.05 #eV
    self.phi_gate 	      = 4.25 #eV
    self.alpha_MI         = 0.7
    self.Lg               = 1e-6
    #UFCM parameters
    self.GEOMOD = 5
    if (self.GEOMOD==5):
      self.Weff = sqrt((self.TFIN_TOP_L-self.TFIN_BASE_L)**2+self.HFIN**2)+sqrt((self.TFIN_TOP_R-self.TFIN_BASE_R)**2+self.HFIN**2)+self.TFIN_TOP_R+self.TFIN_TOP_L
      self.Cins = self.Weff*self.eins/self.tins
      self.Ach = self.HFIN*(self.TFIN_TOP_R+self.TFIN_TOP_L + self.TFIN_BASE_L+self.TFIN_BASE_R)/2
    """self.Cins = (2*self.HFIN+self.TFIN)*self.eins/self.tins
    self.Ach = self.TFIN*self.HFIN
    self.Weff = 2*self.TFIN"""
  def updateparameter(self,name,value):
    #this funtion update a parameter in the model
    if type(value) == type(''):
      exec "self."+name+' = '+value
    else:
      exec "self."+name+' = '+str(value)    

  def draincurrent(self,*args):
    Vdi,Vgi,Vsi,Vbi = args
    """this function returns the drain current for given bias"""
    
    #sweep voltage reference for drain-source
    #bulk reference: TODO: check if this is right
    if (Vdi<Vsi):
      Vg = Vgi-Vbi
      Vs = Vdi-Vbi
      Vd = Vsi-Vbi
      Vb = 0
      flagsweep = 1
    else:
      Vg = Vgi-Vbi
      Vs = Vsi-Vbi
      Vd = Vdi-Vbi
      Vb = 0
      flagsweep = 0
    #source side evaluation  
    Vch = Vs
    qs = UFCMchargemodel.unified_charge_model(Vg,Vch,self.q,self.k,self.T,self.eo,self.eins,self.ech,self.Eg,self.Nc,self.Nv,self.vt,self.ni,self.phi_substrate,self.phi_gate,self.alpha_MI,self.Cins,self.Ach,self.Weff,self.Nch)
    #drain side evaluation  
    Vch = Vd
    qd = UFCMchargemodel.unified_charge_model(Vg,Vch,self.q,self.k,self.T,self.eo,self.eins,self.ech,self.Eg,self.Nc,self.Nv,self.vt,self.ni,self.phi_substrate,self.phi_gate,self.alpha_MI,self.Cins,self.Ach,self.Weff,self.Nch)
    
    ids0 = UFCMdraincurrentmodel.unified_normilized_ids(qs,qd,self.q,self.k,self.T,self.eo,self.eins,self.ech,self.Eg,self.Nc,self.Nv,self.vt,self.ni,self.phi_substrate,self.phi_gate,self.alpha_MI,self.Cins,self.Ach,self.Weff,self.Nch)

    idsfactor = (self.u0*self.vt**2*self.Cins)/self.Lg
    
    if flagsweep ==1:
      ids0=-ids0
      
    return ids0*idsfactor
