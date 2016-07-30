#UFCM: Unified FinFET Compact Model
#author: Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

import UFCMchargemodel
import UFCMdraincurrentmodel
import mobilitymodels
import UFCMidsf1update
#from algopy import UTPM
from scipy.special.orthogonal import p_roots

from numpy import sqrt, exp, log, cosh, pi, tan, sin, cos, arccos, sign, isnan

def cubicdepressedsol(a,b,c,d):
  #a*x**3+b*x**2+c*x+d=0
  if (True):#(abs(a)>1e-9):
    delta = sqrt(-c/(3*a))
    h = 2*a*delta**3
    delta2 = -c/(3*a)
    d2 = d**2
    h2 = 4*a**2*delta2**3
    #print ('d2-h2:',d2-h2)
    if (d2>h2):
      t = (-d+sqrt(d2-h2))/(2.0*a)
      if t>0:
        t=(t)**(1/3)
      else:
        t=-(abs(t))**(1/3) 
      s=-delta2/t
      #print ('qcubic',t-s)
      return t-s,0,0  
    else:
      theta = (1/3)*arccos(-d/h)
      x1 = 2*delta*cos(2*pi/3+theta)
      x2 = 2*delta*cos(2*pi/3-theta) #seems this is the solution, how to check this
      x3 = 2*delta*cos(theta)
      #print (x1,x2,x3)
      return x2,x1,x3
  else:
    return 0,0,0

class compactmodel:
  def __init__(self,name):
    self.name = name
    self.version = 'v1'
    self.nodenames = ['vd','vg','vs','vb']
    
    #TODO: add units to each quantity   
    self.q       = 1.6e-19 #C
    self.k       = 1.380650e-23 
    self.T       = 300 #K
    self.HBAR    = 1.05457e-34 # Joule-sec
    self.MEL     = 9.11e-31 # kg
    self.eo      = 8.854000e-12 # F / m
    self.eins    = 3.9*self.eo # F / m 
    self.ech     = 11.7*self.eo # F / m 
    self.Eg      = 1.169640 #eV
    self.Nc      = 2.890000e+25 #1/m^3
    self.Nv      = 3.140000e+25 #1/m^3 
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
    self.PHIG             = 4.5 #eV
    self.alpha_MI         = 4
    self.Lg               = 1e-6
    self.returnvar        = ['Ids']
    self.QMFACTORCV       = 0
    self.NFIN             = 1
    self.Rs = 0
    self.Rd = 0
    self.countRmodel      = 3
    #UFCM parameters
    self.GEOMOD           = 0
    self.IDSMOD           = 0
    self.DEVTYPE          = 0
    self.RDSMOD           = 0
    self.VTHROMOD         = 0
    self.VTHDIBLMOD       = 0
    self.SSMOD           = 0
    #mobility flags, 0 do not use that model, 1 for use the model
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
    
    self.VTHROFIT = 1
    self.VTHDIBLFIT = 1
    self.lambdafitRO = 1
    self.lambdafitDIBL = 1
    self.lamdafitSS = 1
    self.fitSSall = 1
    
    self.nssshift=100.0
    #mobility parameters for pmos silicon devices
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
    self.lombardi_C = 2.947e3#
    self.lombardi_N0 = 1.0
    self.lombardi_N2 = 1.0
    self.lombardi_lambdam = 0.0317
    self.lombardi_k = 1.0
    self.lombardi_delta = 2.0546e14#
    self.lombardi_eta = 2.0546e30#
    self.lombardi_lcrit = 1e-6 #cm
    self.lombardi_Astar = 2.0
    self.lombardi_Fref = 1    
    
    self.coulomb_D1inv = 135.0 #cm^2/Vs
    self.coulomb_alpha1inv = 0.0

    #self.coulomb_cother = 0.0
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
    
    #remote coulomb scattering, for high-k dielectric
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
    
    self.alpha1_P = -4.8628e+07
    self.alpha11_P = 3.1166e+08
    self.alpha111_P = 0
    self.alpha1111_P = 0
    self.t_FE = 0
    
    self.NCFETMOD = 0
    self.a0 = -1.0 #m/F
    self.b0 = 0#1.3e10 #m^5/F/coul^2
    self.c0 = 0 #m^9/F/coul^4
    
    self.cgsfe = 0
    self.gauss_n = 30
  def updateparameter(self,name,value):
    #this funtion update a parameter in the model
    exec ("self."+name+' = '+'value'   )
    
  def analog(self,*args):
    """this function returns the drain current or other variables for given bias"""  
    
    ###################################################
    ########Bias independent calculations##############
    ###################################################
    #temperature dependent vt, and ni
    self.vt      = self.k*self.T/self.q #V 
    self.ni      = 2.651e+16#sqrt(self.Nc*self.Nv)*exp(-self.Eg/(2*self.vt)) #1/m^3
    gauss_n = self.gauss_n
    
    
    #geometrie Unified FinFET model parameter calculations, TODO: add more geometries
    if (self.GEOMOD==0): #double-gate, rectangular
      self.Cins   = (2.0*self.HFIN)*self.eins/self.tins
      self.Ach    = self.TFIN*self.HFIN
      self.Weff   = 2.0*self.HFIN
    if (self.GEOMOD==1): #triple-gate, rectangular
      self.Cins   = (2.0*self.HFIN+self.TFIN)*self.eins/self.tins
      self.Ach    = self.TFIN*self.HFIN
      self.Weff   = (2.0*self.HFIN+self.TFIN)    
    if (self.GEOMOD==2): #triple-gate, trapezoidal
      self.Weff   = sqrt((self.TFIN_TOP_L-self.TFIN_BASE_L)**2+self.HFIN**2)+sqrt((self.TFIN_TOP_R-self.TFIN_BASE_R)**2+self.HFIN**2)+self.TFIN_TOP_R+self.TFIN_TOP_L
      self.Cins   = self.Weff*self.eins/self.tins
      self.Ach    = self.HFIN*(self.TFIN_TOP_R+self.TFIN_TOP_L + self.TFIN_BASE_L+self.TFIN_BASE_R)/2
    
    #this check if this is PMOS or NMOS device  
    if (self.DEVTYPE==-1): #PMOS
      flagdevtype = -1.0
      PHIG = -self.PHIG+2*self.phi_substrate+self.Eg
    else: #NMOS
      flagdevtype = 1.0
      PHIG = self.PHIG

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
    ###################################################
    ##########Bias dependent calculations##############
    ###################################################
    Vdi,Vgi,Vsi,Vbi,iterationguess = args
    Vgi=Vgi
    count_FE=0
    if True:#while (count_FE<100):
      count_FE+=1
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
        if self.RDSMOD==0:
          Rdaux = self.Rs
          Rsaux = self.Rs      
        else:
          Rdaux = self.Rs
          Rsaux = self.Rd
      else:
        if self.RDSMOD==0:
          Rdaux = self.Rs
          Rsaux = self.Rs      
        else:  
          Rdaux = self.Rd
          Rsaux = self.Rs    
     
      #short channel effect calculations: threshold voltage shift and subthreshold swing degradation
      vfb_n = (PHIG - self.phi_substrate -self.Eg/2.0-self.vt*log(self.Nch/self.ni))/self.vt
      vth_fixed_factor_SI = vfb_n+log(self.Cins*self.vt/(self.q*self.ni**2.0*2.0*self.Ach/self.Nch)) 
      rc  = (2.0*self.Cins/(self.Weff**2.0*self.ech/self.Ach))    
      qdep  = (-self.q*self.Nch*self.Ach)/((self.vt)*self.Cins)
      Vth =  vth_fixed_factor_SI*self.vt+log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))*self.vt
      self.lamda = sqrt(self.ech*self.Ach/self.Cins*(1+rc/2.0))
      
      deltaVth = 0
      #Vth roll-off
      if ( self.VTHROMOD == 0):
        vbi = self.vt*log(1e26*self.Nch/self.ni**2)
        vsl = Vth - vfb_n*self.vt-(self.q*self.Nch/self.ech)*self.lamda**2
        deltaVth = deltaVth + flagdevtype*(2*(vbi-vsl)*self.VTHROFIT/(2*cosh(self.Lg/( self.lambdafitRO*2*self.lamda))-2) )
      else:
        deltaVth = deltaVth + flagdevtype*(self.vthrolloff)

      #Vth DIBL effec
      if ( self.VTHDIBLMOD == 0):
        deltaVth = deltaVth + (-(Vd-Vs)*self.VTHDIBLFIT/(2*cosh(self.Lg/(self.lambdafitDIBL*2*self.lamda))-2) )
      else:
        deltaVth = deltaVth + flagdevtype*(self.vthdibl* (Vd-Vs))
        
      #SS calculation  
      if (self.SSMOD==0):
        vbi = self.vt*log(1e26*self.Nch/self.ni**2) #TODO: input Nd source/drain  
        Ld =   self.lamdafitSS*self.lamda*sqrt(8.0)
        #Ld = 2.0*self.Ach/(self.Weff)+2*(self.ech)*self.Cins/(self.Weff)
        #Ld = 2.0*self.Ach/self.Weff+2.0*4.0*self.tins
        yc = self.Lg/2.0+(Ld/(2*pi))*log(vbi/(vbi+(Vd-Vs)))
        SS = (1.0/0.9)/((4.0/(pi))*sin(pi*yc/self.Lg)/cosh(0.5*pi*Ld/self.Lg)+(4.0/(3.0*pi))*sin(3.0*pi*yc/self.Lg)/cosh(3.0*0.5*pi*Ld/self.Lg)+(4.0/(5.0*pi))*sin(5.0*pi*yc/self.Lg)/cosh(5.0*0.5*pi*Ld/self.Lg))
        #+(4.0/(7.0*pi))*sin(7.0*pi*yc/self.Lg)/cosh(7.0*0.5*pi*Ld/self.Lg)+(4.0/(9.0*pi))*sin(9.0*pi*yc/self.Lg)/cosh(9.0*0.5*pi*Ld/self.Lg)+(4.0/(11.0*pi))*sin(11.0*pi*yc/self.Lg)/cosh(11.0*0.5*pi*Ld/self.Lg)+(4.0/(13.0*pi))*sin(13.0*pi*yc/self.Lg)/cosh(13.0*0.5*pi*Ld/self.Lg)
        nVtm = self.vt#*(SS)
        SS = SS     
      else:  
        SS = 1.0+self.SSrolloff+self.SSdibl* (Vd-Vs)    
        nVtm = self.vt#*(SS)
        SS = SS
   
      #quantum mechanical effects - bias dependence parameter calculations
      mx =  0.916 * self.MEL
      fieldnormalizationfactor  = nVtm*self.Cins/(self.Weff*self.ech)
      auxQMfact  = ((3.0/4.0)*3*self.HBAR*2.0*3.141516*self.q/(4*sqrt(2*mx)))**(2.0/3.0)
      QMf   = self.QMFACTORCV*auxQMfact*(fieldnormalizationfactor)**(2.0/3.0)*(1/(self.q*nVtm))       
         
      #source side evaluation for charge  
      [x,w]=p_roots(gauss_n)
      qmaux = []

      j=0
      while (j<(gauss_n)): 
        if (len(iterationguess)<1):
          guess = float('nan')
        else:
          guess = iterationguess[j]      
        Vch = (Vd-Vs)*(x[j]+1.0)/2.0+Vs
        qs,vfe,qmguess = UFCMchargemodel.unified_charge_model_nc2(self,Vg,Vch,nVtm,PHIG,QMf,SS,Vs,guess)
        qmaux.append(qs)
        exec ('q'+str(j)+'=qs' )
        j=j+1
        
      ids0 = 0
      j=0
      #print (len(qmaux))
      while (j<(gauss_n)): 
        qm = qmaux[j]
        Ft = -((qm+(-q*Nch*Ach)/(self.vt*Cins))*Cins*self.vt/( Weff * ech)) #1e-2 is to transform it to V/cm
        t = (Ft*Ach/Weff+self.vt*(1.0-exp((Ft/self.vt)*Ach/Weff)))/(Ft*(1.0-exp((Ft/self.vt)*Ach/Weff))) #1e3 to transfor to um
        c = -((qm)*Cins*self.vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
        mulow2 = mobilitymodels.mobility(self,DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))      
        ids0=ids0-w[j]*qmaux[j]*mulow2
        j=j+1    
      
      ids0 = 0.5*(Vd-Vs)*ids0*self.Cins*nVtm/self.Lg/Weff/1e6
      '''ids0,mu,vdsat,qd,qdsat = UFCMdraincurrentmodel.unified_normilized_ids(self,qs,nVtm,PHIG,Vd,Vs,Vg,QMf,deltaVth,SS,flagsweep)

      #drain-source current in Ampere [C/s]
      idsfactor = (nVtm**2*self.Cins)/self.Lg
      idsfinal = ids0*idsfactor
      
      #initial guess for drain-source current with source/drain resistances
      if ((Vd-Vs)**2>1e-20):
        idsfinal = idsfinal*(Vd-Vs)/((Vd-Vs)+idsfinal*(Rdaux+Rsaux))
      else:
        idsfinal = idsfinal
      
      #iteration to solve drain-source current including source and drain resistances
      count=0'''
      
      '''while (count<self.countRmodel):  
        Vch = Vs+Rsaux*idsfinal
        qs = UFCMchargemodel.unified_charge_model(self,Vg-deltaVth,Vch,nVtm,PHIG,QMf,SS)
        ids0,mu,vdsat,qd,qdsat = UFCMdraincurrentmodel.unified_normilized_ids(self,qs,nVtm,PHIG,Vd-Rdaux*idsfinal,Vs+Rsaux*idsfinal,Vg,QMf,deltaVth,SS,flagsweep)
        
        #Newton iteration update    
        f0 = idsfinal-ids0*idsfactor
        f1 = UFCMidsf1update.f1(self,qs,qd,nVtm,PHIG,Rsaux,Rdaux)
        idsfinal = idsfinal-f0/f1
        count+=1'''
      
      #source-drain sweep in case Vd<Vs
      '''if flagsweep ==1:
        qaux = qs
        qs = qd
        qd = qs
        idsfinal=-idsfinal
      
      #total current counting number of fins NFIN  
      vedrain = -idsfinal/(qd*self.vt*self.Cins)
      vesource = -idsfinal/(qs*self.vt*self.Cins)
      Ids = flagdevtype*idsfinal*self.NFIN
      Idnorm = Ids/self.Weff*1e-6
      #gate charge, TODO: add source/drain terminal charges 
      Qg = -((qs+qd)*0.5-(qs-qd)**2/(6*(-2*(qs+qd)+1)))*self.Cins*nVtm*self.Lg*self.NFIN*flagdevtype
      
      dqd = -1/(1.0-1.0/qd)/self.vt
      dqs = -1/(1.0-1.0/qs)/self.vt
      dQg_dVg =  -((dqs+dqd)*0.5)*self.Cins*nVtm*self.Lg*self.NFIN*flagdevtype'''
    qs=-qs
    variablesvalues = []
    for var in self.returnvar:
      exec ('variablesvalues.append('+var+')' )
    
    return  variablesvalues, self.returnvar,qmaux
