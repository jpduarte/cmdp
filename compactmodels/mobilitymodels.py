#mobility models
#Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

from numpy import sqrt, exp, log

def mobility(self,DEVTYPE,Ft,Na,Nd,T,t,carrier,Ndepl):
  N_A0 = Na
  N_D0 = Nd
  x = t*1e-4 #convernt um to cm
  #########################################################################################
  #mobility due to phonon scattering
  #########################################################################################
  ul = self.ul #cm^2/Vs
  xi = self.xi
  
  muconst = ul*(T/300.0)**(-xi)

  #########################################################################################
  #Masetti Model: doping dependent
  #########################################################################################
  umin1 = self.masetti_umin1
  umin2 = self.masetti_umin2
  u1    = self.masetti_u1
  Pc    = self.masetti_Pc
  Cr    = self.masetti_Cr
  Cs    = self.masetti_Cs
  alpha = self.masetti_alpha
  beta  = self.masetti_beta

  mudop = umin1*exp(-Pc/(N_A0+N_D0))+(muconst-umin2)/(1.0+((N_A0+N_D0)/Cr)**alpha)-u1/(1.0+(Cs/(N_A0+N_D0))**beta)

  #########################################################################################
  #enhanced lombardi model
  #########################################################################################
  B = self.lombardi_B
  C = self.lombardi_C
  N0 = self.lombardi_N0
  N2 = self.lombardi_N2
  lambdam = self.lombardi_lambdam
  k = self.lombardi_k
  delta = self.lombardi_delta
  eta = self.lombardi_eta
  lcrit = self.lombardi_lcrit
  Astar = self.lombardi_Astar
  Fref = self.lombardi_Fref
  #surface contribution due to acoustic phonon scattering
  muac = (B/Ft + C*((N_A0+N_D0+N2)/N0)**lambdam/(Ft**(1.0/3.0)*(T/300.0)**k) )/exp(-x/lcrit)

  #contribution attributed to surface roughness scattering
  musr = (((Ft/Fref)**Astar/delta+Ft**3.0/eta)**(-1.0))/exp(-x/lcrit)

  #########################################################################################
  #Coulomb Scattering
  #########################################################################################
  Ninv = Nd 
  #t = #um t is the local layer thicknes  
  #c = #concentration mobile chargers per cm^3
  
  D1inv = self.coulomb_D1inv
  alpha1inv = self.coulomb_alpha1inv

  cother = self.coulomb_cother
  v0inv = self.coulomb_v0inv

  v1inv = self.coulomb_v1inv
  D2inv =  self.coulomb_D2inv
  alpha2inv = self.coulomb_alpha2inv
  v2inv = self.coulomb_v2inv
 
  t1 = self.coulomb_t1
  tcoulomb = self.coulomb_tcoulomb
  
  S = self.coulomb_S
  St = self.coulomb_St
  t0 = self.coulomb_t0
  
  m_over_m0 = self.coulomb_m_over_m0
  
  mumax = self.coulomb_mumax
  mumin = self.coulomb_mumin
  alpha = self.coulomb_alpha
  
  Nref = self.coulomb_Nref
  p = self.coulomb_p
  
  muc2dinv = sqrt((D1inv*(T/300.0)**alpha1inv*(carrier/1e18)**v0inv/(Ninv/1e18)**v1inv)**2.0+(D2inv*(T/300.0)**alpha2inv/(Ninv/1e18)**v2inv)**2.0)
  muc2d = muc2dinv#/erf((t+t1)/tcoulomb)
  
  muc = muc2d #fFt*muc3d+(1-fFt)*muc2d TODO: check this effect
  
  #########################################################################################
  #Remote Coulomb Scattering Model
  #########################################################################################
  murcs0  =  self.rcoulomb_murcs0
  gamma1  = self.rcoulomb_gamma1
  gamma2  = self.rcoulomb_gamma2
  gamma3  = self.rcoulomb_gamma3
  gamma4  = self.rcoulomb_gamma4
  gamma5  = self.rcoulomb_gamma5
  s       = self.rcoulomb_s
  c0      = self.rcoulomb_c0
  dcrit   = self.rcoulomb_dcrit
  lcrit   = self.rcoulomb_lcrit
  lcrithk = self.rcoulomb_lcrithk
  xi      = self.rcoulomb_xi
  
  dist = x
  
  fFt = 1.0-exp(-xi*Ft/Ndepl)
  gscreening = s + carrier/(c0*(Nd/3e16)**gamma5)
  Drcs = exp(-(dist+dcrit)/lcrit)
  Drcs_hk = 1.0 #exp(-disthk/lcrithk)
  murcs = (murcs0*(Nd/3e16)**gamma1*(T/300.0)**gamma2*(gscreening)**(gamma3+gamma4*log(Nd/3e16))) #/(Drcs*Drcs_hk)
  #########################################################################################
  muinv = 0.0
  if (self.MODmudop):
    muinv+= 1.0/mudop
  if (self.MODmuac):
    muinv+= 1.0/muac
  if (self.MODmusr):
    muinv+= 1.0/musr
  if (self.MODmuc==1):
    muinv+= 1.0/muc      
  if (self.MODmurcs==1):
    muinv+= 1.0/murcs            
     
  mu = (muinv)**(-1.0)
  return mu*1e-4#,mudop*1e-4,muc*1e-4,muac*1e-4,musr*1e-4,murcs*1e-4#murcs*1e-4 #*1e-4 factor is to transfor in m^2/Vs
  

