#mobility models
from numpy import sqrt, exp, log

def erf(x):
  # save the sign of x
  sign = 1 if x >= 0 else -1
  x = abs(x)

  # constants
  a1 =  0.254829592
  a2 = -0.284496736
  a3 =  1.421413741
  a4 = -1.453152027
  a5 =  1.061405429
  p  =  0.3275911

  # A&S formula 7.1.26
  t = 1.0/(1.0 + p*x)
  y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)
  return sign*y 


def mobility(DEVTYPE,Ft,Na,Nd,T,t,carrier,Ndepl):
  N_A0 = Na
  N_D0 = Nd
  x = t*1e-4 #convernt um to cm
  #########################################################################################
  #mobility due to phonon scattering
  #########################################################################################
  ul = 470.5 #cm^2/Vs
  xi = 2.2
  
  muconst = ul*(T/300.0)**(-xi)

  #########################################################################################
  #Masetti Model: doping dependent
  #########################################################################################
  umin1 = 44.9
  umin2 = 0.0
  u1    = 29.0
  Pc    = 9.23e16
  Cr    = 2.23e17
  Cs    = 6.1e20
  alpha = 0.719
  beta  = 2.0

  mudop = umin1*exp(-Pc/(N_A0+N_D0))+(muconst-umin2)/(1.0+((N_A0+N_D0)/Cr)**alpha)-u1/(1.0+(Cs/(N_A0+N_D0))**beta)

  #########################################################################################
  #enhanced lombardi model
  #########################################################################################
  B = 9.925e6
  C = 2.947e3#15e3#
  N0 = 1.0
  N2 = 1.0
  lambdam = 0.0317
  k = 1.0
  delta = 2.0546e14#0.7e14#
  eta = 2.0546e30#0.1e30#
  lcrit = 1e-6 #cm
  Astar = 2.0
  Fref = 1
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
  
  D1inv = 135.0 #cm^2/Vs
  alpha1inv = 0.0

  cother = 0.0
  v0inv = 1.5

  v1inv = 2.0
  D2inv =  40.0 #cm^2/Vs
  alpha2inv = 0.0
  v2inv = 0.5
 
  t1 = 0.0003 #um
  tcoulomb = 0.0 #um
  
  S = 0.3042 #K(cm/V)^(2/3)
  St = 0.0 #K(um)^2
  t0 = 0.0005 #um
  
  m_over_m0 = 1
  
  mumax = 470.5
  mumin = 44.9
  alpha = 0.719
  
  Nref = 2.23e17 #cm^-3
  p = 4
  
  muc2dinv = sqrt((D1inv*(T/300.0)**alpha1inv*(carrier/1e18)**v0inv/(Ninv/1e18)**v1inv)**2.0+(D2inv*(T/300.0)**alpha2inv/(Ninv/1e18)**v2inv)**2.0)
  muc2d = muc2dinv#/erf((t+t1)/tcoulomb)
  

  """Nsc = Nd+Na+cother
  P = (2.459/(3.97e13*Nsc**(-2.0/3.0))+3.828*(carrier+cother)*m_over_m0/1.36e20)**(-1.0)*(T/300.0)**2.0
  F = (0.7643*P**(0.6478)+ 2.2999 + 6.5502 *m_over_m0)/(P**(0.6478)+2.3670-0.8552*m_over_m0)
  G = 1.0-0.89233/(0.41372+P*(m_over_m0*T/300.0)**0.28227)**0.19778+0.005978/(P*(m_over_m0*T/300.0)**0.72169)**1.80618
  Nsceff = G*Ninv+cother/F 
  mun = (mumax**2/(mumax-mumin))*(T/300.0)**(alpha*3.0-1.5)
  muc = (mumax*mumin/(mumax-mumin))*(T/300.0)**(0.5)
  muc3d = mun*(Nsc/Nsceff)*(Nref/Nsc)**alpha+muc *((carrier+cother)/Nsceff)

  fFt = 1.0/(1.0+exp(S*Ft**(2.0/3.0)/T+St/((t+t0)**2.0*T)-p))"""
  muc = muc2d #fFt*muc3d+(1-fFt)*muc2d TODO: check this effect
  
  #########################################################################################
  #Remote Coulomb Scattering Model
  #########################################################################################
  murcs0  =  149.0 # cm^2/Vs
  gamma1  = -0.23187
  gamma2  = 2.1
  gamma3  = 0.4
  gamma4  = 0.05
  gamma5  = 1.0
  s       = 0.00001#0.1
  c0      = 3.0e16 #cm^-3
  dcrit   = 0.0 # cm
  lcrit   = 1e-6 # cm
  lcrithk = 1e-6 # cm
  xi      = 1.3042e7 # V-1cm-1
  
  dist = x
  
  fFt = 1.0-exp(-xi*Ft/Ndepl)
  gscreening = s + carrier/(c0*(Nd/3e16)**gamma5)
  Drcs = exp(-(dist+dcrit)/lcrit)
  Drcs_hk = 1.0 #exp(-disthk/lcrithk)
  murcs = (murcs0*(Nd/3e16)**gamma1*(T/300.0)**gamma2*(gscreening)**(gamma3+gamma4*log(Nd/3e16))) #/(Drcs*Drcs_hk)
  #########################################################################################
  mu = (1.0/mudop+1.0/muac+1.0/musr+1/murcs)**(-1.0)
  return mu*1e-4,mudop*1e-4,muc*1e-4,muac*1e-4,musr*1e-4,murcs*1e-4#murcs*1e-4 #*1e-4 factor is to transfor in m^2/Vs
  
"""qm = qs
Ft = -((qm+(-self.q*self.Nch*self.Ach)/(nVtm*self.Cins))*self.Cins*nVtm/( self.Weff * self.ech)) #1e-2 is to transform it to V/cm
t = (Ft*self.Ach/self.Weff+nVtm*(1-exp((Ft/nVtm)*self.Ach/self.Weff)))/(Ft*(1-exp((Ft/nVtm)*self.Ach/self.Weff))) #1e3 to transfor to um
c = -((qm)*self.Cins*nVtm/( self.Weff * self.q*(t))) #1e-6 to transform it to cm^-3
mu,mudop,muc,muac,musr = mobilitymodels.mobility(self.DEVTYPE,Ft*1e-2,0,self.Nch*1e-6,self.T,t*1e3,c*1e-6) """   

