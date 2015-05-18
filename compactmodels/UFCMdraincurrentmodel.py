#UFCM drain current model
#Juan Pablo Duarte, jpduarte@berkeley.edu

from numpy import exp,log,sqrt

from scipy.integrate import quad
from scipy.integrate import quadrature
from scipy.special.orthogonal import p_roots

import mobilitymodels

import UFCMvdsat
import UFCMchargemodel

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

def unified_normilized_ids(qms,qmd,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch,IDSMOD,DEVTYPE,Lg,Vds,Vg,QMFACTORCVfinal,deltaVth) :
  #normalized drain current model: Unified FinFET Compact Model: Modelling Trapezoidal Triple-Gate FinFETs
  Vdseff = 0
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
    mulow2,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    ids2 = -qm*(1-1/qm)*mu*w2    
    
    qm = q3
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1-exp((Ft/vt)*Ach/Weff)))/(Ft*(1-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mulow2,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6*(Weff*1e2))
    ids3 = -qm*(1-1/qm)*mu*w3     
    ids = (qmd-qms)/2.0*(ids1+ids2+ids3)
    
    Fhfs1 = Vds/Lg
    Fhfs2 = Vds/Lg
    Fhfs3 = Vds/Lg
    v1 = (ids/qms)*(vt)/Lg
    v2 = (ids/q2)*(vt)/Lg
    v3 = (ids/qmd)*(vt)/Lg    
    
  if (IDSMOD==2):
    x1 = -7.74596669e-01
    x2 = -4.78946310e-17
    x3 = 7.74596669e-01
    w1 = 0.55555556
    w2 = 0.88888889
    w3 = 0.55555556
    q1 = (qmd-qms)*(x1+1)/2.0 + qms
    q2 = (qmd-qms)*(x2+1)/2.0 + qms
    q3 = (qmd-qms)*(x3+1)/2.0 + qms
    
    me = 9.11e-31 #kg
    vsat = 8.37e6  *1e-2 # m/s 8.37e6
    taue = 0.25e-13 #energy relaxation time [s],0.25e-12
    fieldfactor = 1.0;
     
    ####################################################################################################
    # calculation without considering current saturation
    ####################################################################################################       
    qm = q1
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1-exp((Ft/vt)*Ach/Weff)))/(Ft*(1-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mulow1,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    ids1 = -qm*(1-1/qm)*mulow1*w1
    
    qm = q2
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1-exp((Ft/vt)*Ach/Weff)))/(Ft*(1-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mulow2,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    ids2 = -qm*(1-1/qm)*mulow2*w2    
    
    qm = q3
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1-exp((Ft/vt)*Ach/Weff)))/(Ft*(1-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mulow3,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6*(Weff*1e2))
    ids3 = -qm*(1-1/qm)*mulow3*w3     
    
    ids = 0.5*(qmd-qms)*(ids1+ids2+ids3)   
    ####################################################################################################
    # calculation considering current saturation
    ####################################################################################################
    beta=2.0 
    alpha = 0.0    
    vd = (ids/qmd)*(vt)/Lg
    v3 = vd
    Dsat1 = (alpha+1)/((1.0+( (alpha+1)*vd/vsat)**beta)**(1.0/beta)+alpha)
    ids = 0.5*(qmd-qms)*(ids1+ids2+ids3)*Dsat1
    vs = (ids/qms)*(vt)/Lg
    v1 = vs
    #v3 = (ids/qmd)*(vt)/Lg
    v2 = (v3+v1)*0.5
    
    Fhfs1 = v1/mulow1
    Fhfs2 = v2/mulow2
    Fhfs3 = v3/mulow3      
 
    mu  = mulow1 

  if (IDSMOD==3):
    x1 = -0.57735027
    x2 = 0.57735027
    w1 = 1.0
    w2 = 1.0    
    q1 = (qmd-qms)*(x1+1)/2.0 + qms
    q2 = (qmd-qms)*(x2+1)/2.0 + qms
    
    qm = q1
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1.0-exp((Ft/vt)*Ach/Weff)))/(Ft*(1.0-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mulow1,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    ids1 = -qm*(1.0-1.0/qm)*mulow1*w1
    
    mulow2 =mulow1
    ####################################################################################################
    # calculation considering current saturation
    #################################################################################################### 
    beta=2.0 
    vsat = 2*8.37e6*1e-2#8.37e6*1e-2 # m/s 8.37e6            
    qdsat = (2.0*Lg*vsat + mulow1*qms*vt*x1 + mulow1*vt + mulow2*qms*vt*x2 + mulow2*vt - sqrt(4.0*Lg**2*vsat**2 + 4.0*Lg*mulow1*qms*vsat*vt*x1 + 4.0*Lg*mulow1*vsat*vt + 4.0*Lg*mulow2*qms*vsat*vt*x2 + 4.0*Lg*mulow2*vsat*vt + mulow1**2*qms**2*vt**2 - 2.0*mulow1**2*qms*vt**2 + mulow1**2*vt**2 + 2.0*mulow1*mulow2*qms**2*vt**2 - 4.0*mulow1*mulow2*qms*vt**2 + 2.0*mulow1*mulow2*vt**2 + mulow2**2*qms**2*vt**2 - 2.0*mulow2**2*qms*vt**2 + mulow2**2*vt**2))/(mulow1*vt*x1 + mulow1*vt + mulow2*vt*x2 + mulow2*vt)
     
    Vdsat = UFCMvdsat.vdsat(Vg-deltaVth,qdsat,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch,QMFACTORCVfinal) 
    beta = 2.0
    if (Vds>1e-8):
      Vdseff = 1.0/((1.0/Vdsat)**beta+(1.0/Vds)**beta)**(1.0/beta)
    else:
      Vdseff = Vds
    qmd = UFCMchargemodel.unified_charge_model(Vg-deltaVth,Vdseff,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch,QMFACTORCVfinal)
    
    q1 = (qmd-qms)*(x1+1)/2.0 + qms
    q2 = (qmd-qms)*(x2+1)/2.0 + qms
    
    qm = q1
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1.0-exp((Ft/vt)*Ach/Weff)))/(Ft*(1.0-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mulow1,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    ids1 = -qm*(1.0-1.0/qm)*mulow1*w1
    
    qm = q2
    Ft = -((qm+(-q*Nch*Ach)/(vt*Cins))*Cins*vt/( Weff * ech)) #1e-2 is to transform it to V/cm
    t = (Ft*Ach/Weff+vt*(1.0-exp((Ft/vt)*Ach/Weff)))/(Ft*(1.0-exp((Ft/vt)*Ach/Weff))) #1e3 to transfor to um
    c = -((qm)*Cins*vt/( Weff * q*(t))) #1e-6 to transform it to cm^-3
    mulow2,mudop,muc,muac,musr,murcs = mobilitymodels.mobility(DEVTYPE,Ft*1e-2,0,Nch*1e-6,T,t*1e3,c*1e-6,Nch*1e-6/(Weff*1e2))
    ids2 = -qm*(1.0-1.0/qm)*mulow2*w2   
   
    ids = 0.5*(qmd-qms)*(ids1+ids2)    

    mu  = mulow1
    
  return ids,mu,Vdseff
   
