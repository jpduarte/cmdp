#UFCM drain current model
#Juan Pablo Duarte, jpduarte@berkeley.edu

from numpy import exp,log,sqrt

from scipy.integrate import quad
from scipy.integrate import quadrature
from scipy.special.orthogonal import p_roots

import mobilitymodels

def fixed_quad(func,a,b,n,args=()) :
  [x,w] = p_roots(n)
  x = real(x)
  y = (b-a)*(x+1)/2.0 + a
  sum_integral = 0
  i=0
  for xi in x:
    sum_integral = sum_integral+w[i]*func((b-a)*(xi+1)/2.0 + a,*args)
    i=i+1
  return (b-a)/2.0*sum_integral

def current_integral(qm, a, b, c,eps):
     return -qm*(1-1/qm)/(1-a*qm+b*qm**2-c*(1/(eps+qm)))

def unified_normilized_ids(qms,qmd,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch,IDSMOD,DEVTYPE) :
  #normalized drain current model: Unified FinFET Compact Model: Modelling Trapezoidal Triple-Gate FinFETs
  
  if (IDSMOD==0):
    rc  = (2*Cins/(Weff**2*ech/Ach))#TODO: check factor of 2 in rc?
    qdep  = (-q*Nch*Ach)/(vt*Cins)
    qh = 1/rc-qdep
    qm = qms
    idss = qm**2/2-2*qm-qh*log(1-qm/qh)
    qm = qmd
    idsd = qm**2/2-2*qm-qh*log(1-qm/qh) 
    qm = (qms+qmd)/2
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1-exp((Ft/vt)*Ach/Weff)))/(Ft*(1-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mu,mudop,muc,muac,musr = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))   
    ids =  mu*(idss-idsd )
  
  if (IDSMOD==1):
    #ids, error = quad(current_integral, qms, qmd, args=(0.1,0.001,1,-0.1))
    #fixed_quad(current_integral, qs, qd, 3,[0.1,0.001,1,-0.1])
    x1 = -7.74596669e-01
    x2 = -4.78946310e-17
    x3 = 7.74596669e-01
    w1 = 0.55555556
    w2 = 0.88888889
    w3 = 0.55555556
    q1 = (qmd-qms)*(x1+1)/2.0 + qms
    q2 = (qmd-qms)*(x2+1)/2.0 + qms
    q3 = (qmd-qms)*(x3+1)/2.0 + qms
    
    qm = q1
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1-exp((Ft/vt)*Ach/Weff)))/(Ft*(1-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mu,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    ids1 = -qm*(1-1/qm)*mu*w1
    
    qm = q2
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1-exp((Ft/vt)*Ach/Weff)))/(Ft*(1-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mu,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    ids2 = -qm*(1-1/qm)*mu*w2    
    
    qm = q3
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1-exp((Ft/vt)*Ach/Weff)))/(Ft*(1-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mu,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6*(Weff*1e2))
    ids3 = -qm*(1-1/qm)*mu*w3     
    ids = (qmd-qms)/2.0*(ids1+ids2+ids3)

  return ids,mu
