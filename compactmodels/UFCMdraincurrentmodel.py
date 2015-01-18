#UFCM drain current model
#Juan Pablo Duarte, jpduarte@berkeley.edu

from numpy import exp,log,sqrt

def unified_normilized_ids(qms,qmd,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch) :
  #normalized drain current model: Unified FinFET Compact Model: Modelling Trapezoidal Triple-Gate FinFETs
  rc  = (2*Cins/(Weff**2*ech/Ach))#TODO: check factor of 2 in rc?
  qdep  = (-q*Nch*Ach)/(vt*Cins)
  qh = 1/rc-qdep
  qm = qms
  idss = qm**2/2-2*qm-qh*log(1-qm/qh)
  qm = qmd
  idsd = qm**2/2-2*qm-qh*log(1-qm/qh)  
  return idss-idsd  
