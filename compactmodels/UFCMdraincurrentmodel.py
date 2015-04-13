#UFCM drain current model
#Juan Pablo Duarte, jpduarte@berkeley.edu

from numpy import exp,log,sqrt

from scipy.integrate import quad
from scipy.integrate import quadrature
from scipy.special.orthogonal import p_roots

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

def unified_normilized_ids(qms,qmd,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch,IDSMOD) :
  #normalized drain current model: Unified FinFET Compact Model: Modelling Trapezoidal Triple-Gate FinFETs
  
  if (IDSMOD==0):
    rc  = (2*Cins/(Weff**2*ech/Ach))#TODO: check factor of 2 in rc?
    qdep  = (-q*Nch*Ach)/(vt*Cins)
    qh = 1/rc-qdep
    qm = qms
    idss = qm**2/2-2*qm-qh*log(1-qm/qh)
    qm = qmd
    idsd = qm**2/2-2*qm-qh*log(1-qm/qh) 
    ids =  idss-idsd 
  
  if (IDSMOD==1):
    ids, error = quad(current_integral, qms, qmd, args=(0.1,0.001,1,-0.1))
    #fixed_quad(current_integral, qs, qd, 3,[0.1,0.001,1,-0.1])
  return ids
