#UFCM charge, qd sat
"""For a given bias and source charge, this function returns the effective
normalized drain charge including current saturation effects
"""
from numpy import exp,log,sqrt

def unified_charge_model_sat(Vg,Vch,q,k,T,eo,eins,ech,Eg,Nc,Nv,vt,ni,phi_substrate,phi_gate,alpha_MI,Cins,Ach,Weff,Nch,QMF,qs) :
"""  vt              = MODEL.PhysicalParameter.vt;
  Lg              = MODEL.DeviceDimension.Lg;
  um_fit_factor   = MODEL.FitParameter.um_fit_factor;
  um              = MODEL.PhysicalParameter.um_e*um_fit_factor;
  a               = MODEL.FitParameter.alpha_mob_fit_factor*MODEL.FitParameter.mob_fit_factor ;
  b               = MODEL.FitParameter.beta_mob_fit_factor*MODEL.FitParameter.mob_fit_factor ;
  vsat            = abs(MODEL.PhysicalParameter.vsat*MODEL.FitParameter.vsat_factor) ;
  %vsource         = MODEL.PhysicalParameter.vsource*MODEL.FitParameter.vsat_factor ;
  qdep            = MODEL.ChargeCalcParameters.qdep;
  Ids_Factor_body_effect  = MODEL.ChargeCalcParameters.Ids_Factor_body_effect;
"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Initial Guess Qs_sat%%%%%%%%%%%%%%%%%%%%%%%%%
um_eff          = um*(1+a*qs)^(-1);
qd =  2 - ((Lg^2*vsat^2 + 4*Lg*um_eff*vsat*vt + qs^2*um_eff^2*vt^2 - ...
    4*qs*um_eff^2*vt^2 + 4*um_eff^2*vt^2)^(1*2^(-1)) - Lg*vsat)*(um_eff*vt)^(-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Initial Guess Ids%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch MODEL.MODE.IV_Model
    case 0 %contant mobility, original dvch/dqm
        %Ids = int(qm*(1-1/qm-1/(qm+gama*qdep)+(rc/gama)/((qm/(gama)+qdep)*rc-((qm/(gama)+qdep)*rc)^2)) , qm)
        um_eff          = um;
    case 1 %contant mobility but with degradation, original dvch/dqm
        %Ids = int(qm*(1-1/qm-1/(qm+gama*qdep)+(rc/gama)/((qm/(gama)+qdep)*rc-((qm/(gama)+qdep)*rc)^2)) , qm)
        um_eff          = um*(1+a*(qdep+(qs+qd)*(2^(-1)))+(b*(qdep+(qs+qd)*2^(-1)))^2)^(-1);
    case 2 %mobility degradation and apprximate dvch/dqm
        %,Ids = int((1+a*(qm))^(-1)*qm*(1-1/qm) , qm)
        um_eff          = um*(1+a*qd)^(-1);
    case 3 %mobility degradation and apprximate dvch/dqm
        %Ids = int((1+a*(qm+qdep))^(-1)*qm*(1-1/qm) , qm)
        um_eff          = um*(1+a*(qd+qdep))^(-1);
    case 4 %mobility degradation and apprximate dvch/dqm
        %Ids = int((1+a*(qm)+(b*qm)^2)^(-1)*qm*(1-1/qm) , qm)
        um_eff          = um*(1+a*(qd)+(b*qd)^2)^(-1);
    case 5 %mobility degradation and apprximate dvch/dqm
        %Ids = int((1+a*(qm+qdep)+(b*(qm+qdep))^2)^(-1)*qm*(1-1/qm) , qm)
        um_eff          = um*(1+a*(qd+qdep)+(b*(qd+qdep))^2)^(-1);
    case 6 %mobility degradation and original dvch/dqm
        %Ids = int((1+a*(qm))^(-1)*qm*(1-1/qm-1/(qm+gama*qdep)+(rc/gama)/((qm/(gama)+qdep)*rc-((qm/(gama)+qdep)*rc)^2)) , qm)
        um_eff          = um*(1+a*qd)^(-1);
    case 7 %mobility degradation and original dvch/dqm
        %Ids = int((1+a*(qm+qdep))^(-1)*qm*(1-1/qm-1/(qm+gama*qdep)+(rc/gama)/((qm/(gama)+qdep)*rc-((qm/(gama)+qdep)*rc)^2)) , qm)
        um_eff          = um*(1+a*(qd+qdep))^(-1);
    otherwise
        warning('Unexpected I-V Model');
end
Ids_guess       = Ids_Factor_body_effect*(-um*vt*Lg^(-1))*current_unified_normalized(qs,qd,MODEL);
f0              = -qd*vsat-Ids_guess;
%um_eff          = um/(1+a*qd);
f1              = -vsat+Ids_Factor_body_effect*(um_eff*vt*Lg^(-1))*(qd-2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Qd_sat Final Value%%%%%%%%%%%%%%%%%%%%%%%%%%%
qd              = (qd-f0*f1^(-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Vd Efeective%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vd_sat         = Vd_SAT_model(Vg, Vd, Vb, n_SS, qd, MODEL);
Vd_eff         = Vd_eff_Model(Vg,Vd,Vs,Vb,Vd_sat,MODEL);
%%%%%%%%%%%%%Calculation of effectuve Charge Saturation at Drain%%%%%%%%%%
[qd_sat_Id, qd_sat_CV, Vth_bias]          = terminal_charge_function(Vg,Vd_eff,Vb,n_SS,MODEL);%Charge_Duarte_Universal_Explicit_2_MODEL((Vg-Vd_eff)/vt,n_SS,MODEL)/(vt*Cins);
%end

end

