#UFCM charge model it is based on; Unified FinFET Compact Model: Modelling Trapezoidal Triple-Gate FinFETs
#this includes quantum mechanical bias dependent effects
#Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

from numpy import exp,log,sqrt,isnan

def unified_charge_model(self,Vg,Vch,vt,phi_gate,QMF,nss) : 
  vt = vt
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


  Vg = Vg+log(self.nssshift)*(nss-1.0)*vt
  rc  = (2.0*Cins/(Weff**2.0*ech/Ach))
  qdep  = (-q*Nch*Ach)/((vt/nss)*Cins)
  vfb_n = (phi_gate - phi_substrate -Eg/2.0-vt*log(Nch/ni))/vt
  vth_fixed_factor_SI = vfb_n+log(Cins*vt/(q*ni**2.0*2.0*Ach/Nch))
  vth_fixed_factor_Sub = log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))+vth_fixed_factor_SI#updated for low doping, log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))
  vth_N_Sub = -qdep+vth_fixed_factor_Sub+QMF*((-qdep)**(2.0/3.0))
  vth_N_SI  = -qdep+vth_fixed_factor_SI+QMF*((-qdep)**(2.0/3.0))
  Vg_local_N  = (Vg-Vch)*vt**(-1.0)
  F           = -Vg_local_N+vth_N_SI
  Vov         = (Vg_local_N-vth_N_Sub)*0.5/nss
  if (Vov>60):
    qm      = -Vov*2.0*nss
    qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
    x0      = qtrc/(exp(qtrc)-qtrc-1.0)
    x1      = qtrc*x0
    
    f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
    f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
    f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
    
    f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
    f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
    f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))

    f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
    f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
    f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))       
  else:
    qm      = exp((Vg_local_N-vth_N_Sub)*0.5/nss)
    if(qm>1.0e-7):#original 1.0e-7
      qm      = 2.0*(1.0-sqrt(1.0+(log(1.0+qm))**2.0))
      qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
      x0      = qtrc/(exp(qtrc)-qtrc-1.0)
      x1      = qtrc*x0
      
      f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
      f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
      f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
        
      f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
      f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
      f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0)) 
        
      f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))
      f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))
      f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))       
    else:
      qm      = -qm**2.0
  return qm  
  

def unified_charge_model_nc(self,Vg,Vch,vt,phi_gate,QMF,nss) :
  vt = vt
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
  a0 = self.a0
  b0 = self.b0
  c0 = self.c0

  Vg = Vg+log(self.nssshift)*(nss-1.0)*vt
  rc  = (2.0*Cins/(Weff**2.0*ech/Ach))
  qdep  = (-q*Nch*Ach)/((vt/nss)*Cins)
  vfb_n = (phi_gate - phi_substrate -Eg/2.0-vt*log(Nch/ni))/vt
  vth_fixed_factor_SI = vfb_n+log(Cins*vt/(q*ni**2.0*2.0*Ach/Nch))
  vth_fixed_factor_Sub = log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))+vth_fixed_factor_SI#updated for low doping, log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))
  vth_N_Sub = -qdep+vth_fixed_factor_Sub+QMF*((-qdep)**(2.0/3.0))
  vth_N_SI  = -qdep+vth_fixed_factor_SI+QMF*((-qdep)**(2.0/3.0))
  Vg_local_N  = (Vg-Vch)*vt**(-1.0)
  F           = -Vg_local_N+vth_N_SI
  Vov         = (Vg_local_N-vth_N_Sub)*0.5/nss
  if (Vov>60):
    qm      = -Vov*2.0*nss
    qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
    x0      = qtrc/(exp(qtrc)-qtrc-1.0)
    x1      = qtrc*x0
    
    f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))+a0*qm+b0*qm**3.0+c0*qm**5.0
    f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))+a0+3.0*b0*qm**2.0+5.0*c0*qm**4.0
    f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))+6.0*b0*qm+20.0*c0*qm**3.0
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
    
    f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))+a0*qm+b0*qm**3.0+c0*qm**5.0
    f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))+a0+3.0*b0*qm**2.0+5.0*c0*qm**4.0
    f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))+6.0*b0*qm+20.0*c0*qm**3.0
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
    
    f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))+a0*qm+b0*qm**3.0+c0*qm**5.0
    f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))+a0+3.0*b0*qm**2.0+5.0*c0*qm**4.0
    f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))+6.0*b0*qm+20.0*c0*qm**3.0
    qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
  else:
    qm      = exp((Vg_local_N-vth_N_Sub)*0.5/nss)
    if(qm>1.0e-7):#original 1.0e-7
      qm      = 2.0*(1.0-sqrt(1.0+(log(1.0+qm))**2.0))
      qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
      x0      = qtrc/(exp(qtrc)-qtrc-1.0)
      x1      = qtrc*x0

      f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))+a0*qm+b0*qm**3.0+c0*qm**5.0
      f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))+a0+3.0*b0*qm**2.0+5.0*c0*qm**4.0
      f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))+6.0*b0*qm+20.0*c0*qm**3.0
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
      
      f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))+a0*qm+b0*qm**3.0+c0*qm**5.0
      f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))+a0+3.0*b0*qm**2.0+5.0*c0*qm**4.0
      f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))+6.0*b0*qm+20.0*c0*qm**3.0
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
      
      f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))+a0*qm+b0*qm**3.0+c0*qm**5.0
      f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))+a0+3.0*b0*qm**2.0+5.0*c0*qm**4.0
      f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))+6.0*b0*qm+20.0*c0*qm**3.0
      qm      = qm-(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
    else:
      qm      = -qm**2.0
  return qm  
    
def unified_charge_model_nc2(self,Vg,Vch,vt,phi_gate,QMF,nss,Vs,qguess_iteration) : 
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
        alpha1_P = self.alpha1_P 
        alpha11_P = self.alpha11_P
        alpha111_P = self.alpha111_P 
        t_FE = self.t_FE         
        a0 = 2.0*alpha1_P*t_FE*Cins/Weff;
        b0 = 4.0*alpha11_P*t_FE*(vt*Cins/Weff)**3/vt;
        c0 = 6.0*alpha111_P*t_FE*(vt*Cins/Weff)**5/vt;     

        cgsfe = self.cgsfe

        Vg = Vg+log(self.nssshift)*(nss-1.0)*vt
        rc  = (2.0*Cins/(Weff**2.0*ech/Ach))
        qdep  = (-q*Nch*Ach)/((vt/nss)*Cins)
        vfb_n = (phi_gate - phi_substrate -Eg/2.0-vt*log(Nch/ni))/vt
        vth_fixed_factor_SI = vfb_n+log(Cins*vt/(q*ni**2.0*2.0*Ach/Nch))
        vth_fixed_factor_Sub = log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))+vth_fixed_factor_SI#updated for low doping, log((qdep*rc)**2.0/(exp(qdep*rc)-qdep*rc-1.0))
        vth_N_Sub = -qdep+vth_fixed_factor_Sub+QMF*((-qdep)**(2.0/3.0))
        vth_N_SI  = -qdep+vth_fixed_factor_SI+QMF*((-qdep)**(2.0/3.0))
        Vg_local_N  = (Vg-Vch)*vt**(-1.0)
        F           = -Vg_local_N+vth_N_SI
        Vov         = (Vg_local_N-vth_N_Sub)*0.5/nss
        
        
        ########################################################Initial Guess: Begin##################################################
        qgsfe = -cgsfe*(Vg-Vs)
        auxqm = alpha1_P/(2.0*alpha11_P*(vt*Cins/Weff)**2)
        qm = -sqrt(-auxqm)
        qminter = -sqrt(-auxqm)
        qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
        x0      = qtrc/(exp(qtrc)-qtrc-1.0)
        x1      = qtrc*x0  
        F           = -Vg_local_N+vth_N_SI
        Vgn=(vth_N_SI-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0)))*vt
        dvfe = -(a0+3.0*b0*(qm+qgsfe)**2.0+5.0*c0*(qm+qgsfe)**4.0)
        ddvfe = -(6.0*b0*(qm+qgsfe)+20.0*c0*(qm+qgsfe)**3.0)
        Coxn = -1/(-1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))+dvfe)
        Coxnaux = (-(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))+ddvfe)
        Coxn2 = -Coxnaux*Coxn**2
        
        #print (a0,b0)
        qmfequesssub = exp(Vg_local_N-vth_N_Sub)
        vfes = -(a0*(-qmfequesssub+qgsfe) )   
        qmfeguess = -qm + Coxn*(Vg/vt-Vgn/vt) #+ Coxn2*(Vg/vt-Vgn/vt)**2.0/2.0 #Coxn*(Vg/vt-vthnfe)
        mtfe = 1.0
        a = qmfeguess
        b = 1.0
        c = 1.0
        d = 2.0
        qmfeguess = a+b-0.5*(a+b-sqrt((a-b)**2+c))  
        ################calculating the point of vfe max
        auxqm = alpha1_P/(2.0*alpha11_P*(vt*Cins/Weff)**2)
        qmvfemax = -sqrt(-auxqm/3.0)
        vfemax = -(a0*(qmvfemax+qgsfe)+b0*(qmvfemax+qgsfe)**3.0)
        qm = qmvfemax
        Vgn2=(vth_N_SI-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))+vfemax)*vt
        
        m = -(qminter-qmvfemax)/((Vgn-Vgn2)*(1/vt))
        vthnfe = Vgn2/vt+qmvfemax/m
        qmgues2 = m*(Vg/vt-vthnfe)
        #print (qminter,qmvfemax,Vgn,Vgn2,vthnfe*vt) 
        a = qmgues2
        b = 1.0
        c = 1.0
        d = 2.0
        qmfeguess = a+b-0.5*(a+b-sqrt((a-b)**2+c))     
        qmgues2 =     qmfeguess
               
        ###########################################3  
        qmfequesssub = qmfequesssub*exp(-vfes)
        vfes = -(a0*(-qmfequesssub+qgsfe) ) 
        qmfequesssub = exp(Vg_local_N-vth_N_Sub)*exp(-vfes)
        if qmfequesssub>1e3:
          qmfequesssub = 10000
        qmfeguess = qmfeguess*qmfequesssub/(qmfequesssub+qmfeguess)
        qmfeguessprint = qmfeguess

        ############################################################Initial Guess: End############################################    
        
        #qmfe,qmfe1,qmfe3 = cubicdepressedsol(-b0/nss,0,(-a0-1.0)/nss,-(Vg_local_N-vth_N_Sub)/nss )
        #print ('d:',-(Vg_local_N-vth_N_Sub)/nss)
        #qmlin = -(Vg_local_N-vth_N_Sub)/nss
        
        if not (isnan(qguess_iteration)):
         qm = qguess_iteration
        else:
         qm = -qmfeguess
 
        #qmguess = 0 
        if (True):
          i=1  
          while (i<20):
             i = i+1
             vfe = -(a0*qm+b0*qm**3.0+c0*qm**5.0)
             dvfe = -(a0+3.0*b0*(qm+qgsfe)**2.0+5.0*c0*(qm+qgsfe)**4.0)
             Vov         = (Vg_local_N-vfe-vth_N_Sub)*0.5/nss  
             #qm =  2.0*(1.0-sqrt(1.0+(log(1.0+exp(Vov)))**2.0))
             vfe = -(a0*(qm+qgsfe)+b0*(qm+qgsfe)**3.0+c0*(qm+qgsfe)**5.0 )
             dvfe = -(a0+3.0*b0*(qm+qgsfe)**2.0+5.0*c0*(qm+qgsfe)**4.0)
             ddvfe = -(6.0*b0*(qm+qgsfe)+20.0*c0*(qm+qgsfe)**3.0)
            
             qtrc    = (qm*alpha_MI**(-1.0)+qdep)*rc
             x0      = qtrc/(exp(qtrc)-qtrc-1.0)
             x1      = qtrc*x0          
            
             f0      = F-qm+log(-qm)*nss+log(x1)+QMF*((-(qdep+qm))**(2.0/3.0))+vfe
             f1      = -1.0+qm**(-1.0)*nss+(2.0*qtrc**(-1.0)-x0-1.0)*rc-(2.0/3.0)*QMF*((-(qdep+qm))**(-1.0/3.0))+ dvfe
             f2      = -(qm**2.0)**(-1.0)*nss-(2.0/9.0)*QMF*((-(qdep+qm))**(-4/3.0))+ddvfe
             delta = -(f0*f1**(-1.0))*(1.0+(f0*f2)*(2.0*f1**2.0)**(-1.0))
             
             #
             
             
             if (abs(delta)<10):
              qm      = qm+delta
             else:
              qm = -qmfeguess
              
             if (qm>0):
              qm = -qmfeguess
          #print (f0,f1,f2,delta,qm)
          ######################################
          vfe = -(a0*(qm+qgsfe)+b0*(qm+qgsfe)**3.0+c0*(qm+qgsfe)**5.0 )
        return qm,vfe,qmfeguess 
