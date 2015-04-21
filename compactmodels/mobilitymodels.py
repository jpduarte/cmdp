#mobility models
from numpy import sqrt, exp

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


def mobility(DEVTYPE,Ft,Na,Nd,T,t,carrier):
  N_A0 = Na
  N_D0 = Nd
  x = t*1e-4
  #mobility due to phonon scattering
  ul = 470.5 #cm^2/Vs
  xi = 2.2
  
  muconst = ul*(T/300.0)**(-xi)

  #Masetti Model: doping dependent
  umin1 = 44.9
  umin2 = 0.0
  u1 = 29.0
  Pc = 9.23e16
  Cr = 2.23e17
  Cs = 6.1e20
  alpha = 0.719
  beta = 2.0

  mudop = umin1*exp(-Pc/(N_A0+N_D0))+(muconst-umin2)/(1.0+((N_A0+N_D0)/Cr)**alpha)-u1/(1.0+(Cs/(N_A0+N_D0))**beta)

  #enhanced lombardi model
  B = 9.925e6
  C = 15e3#2.947e3
  N0 = 1.0
  N2 = 1.0
  lambdam = 0.0317
  k = 1.0
  delta = 0.7e14#2.0546e14
  eta = 0.1e30#2.0546e30
  lcrit = 1e-6 #cm
  Astar = 2.0
  Fref = 1
  #surface contribution due to acoustic phonon scattering
  muac = (B/Ft + C*((N_A0+N_D0+N2)/N0)**lambdam/(Ft**(1.0/3.0)*(T/300.0)**k) )/exp(-x/lcrit)

  #contribution attributed to surface roughness scattering
  musr = (((Ft/Fref)**Astar/delta+Ft**3.0/eta)**(-1.0))/exp(-x/lcrit)

  #Coulomb Scattering
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
  
  mu = (1.0/mudop+1.0/muac+1.0/musr)**(-1.0)
  return mu*1e-4,mudop*1e-4,muc*1e-4,muac*1e-4,musr*1e-4 #*1e-4 factor is to transfor in m^2/Vs

