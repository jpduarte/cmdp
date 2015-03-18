import numpy as np
import supportfunctions as sf
import shutil #this library is to copy file from original to AUX for simulation

"""wheretosimpath ='/users/jpduarte/research/variabilityproject/netlist_modelcards/idvgscaling/'
templatepath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/idvgnmostemplate.sp"
modelverilogpath = "\"/users/jpduarte/research/variabilityproject/model_code/code_109beta/bsimcmg.va\""
modelcardpath = "/users/jpduarte/research/variabilityproject/netlist_modelcards/templateshspice/modelcard-109template.nmos"
vgs = np.linspace(0,0.86,50)
vds = np.linspace(0.05,0.86,3)
Lparam = '20e-9' 
NFINparam = '1'"""

def preparehspiceidvg(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,Lparam,NFINparam, DEVTYPEparam):
  """*Sample netlist for BSIM-MG
  *Id-Vg Characteristics for NMOS (T = 27 C)

  .option abstol=1e-6 reltol=1e-6 post ingold
  .temp 27

  .hdl pathmodelverilog
  .include pathmodelcard

  * --- Voltage Sources ---
  vds supply  0 dc=0.05
  vgs gate  0 dc=1
  vbs bulk  0 dc=0
  * --- Transistor ---
  X1 supply gate 0 bulk nmos1 L=Lparam NFIN=NFINparam 
  * --- DC Analysis ---
  .dc vgs vgsi vgsf vgsdelta vds vds vdsi vdsf vdsdelta
  .print dc i(X1.d)
  .end"""
  #make an aux copy of hspice file to simulate
  shutil.copyfile(templatepath,wheretosimpath+'idvgaux.sp')
  #make an aux copy of modelcard file to simulate
  shutil.copyfile(modelcardpath,wheretosimpath+'modelcardaux.nmos')

  #update path of model and modelcard
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'pathmodelverilog', modelverilogpath)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'pathmodelcard', '\"modelcardaux.nmos\"')

  #bias update
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsi', str(vgs[0]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsf', str(vgs[-1]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsdelta', str(vgs[1]-vgs[0]))

  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsi', str(vds[0]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsf', str(vds[-1]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsdelta', str(vds[1]-vds[0]))

  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'Lparam', Lparam)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'NFINparam',NFINparam)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'DEVTYPEparam',DEVTYPEparam)  
def preparehspiceidvgGEO1(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,Lparam,HFINparam,TFIN_TOPparam,TFIN_BASEparam,EOTparam,NBODYparam,NFINparam):
  """*Sample netlist for BSIM-MG
  *Id-Vg Characteristics for NMOS (T = 27 C)

  .option abstol=1e-6 reltol=1e-6 post ingold
  .temp 27

  .hdl pathmodelverilog
  .include pathmodelcard

  * --- Voltage Sources ---
  vds supply  0 dc=0.05
  vgs gate  0 dc=1
  vbs bulk  0 dc=0
  * --- Transistor ---
  X1 supply gate 0 bulk nmos1 L=Lparam NFIN=NFINparam 
  * --- DC Analysis ---
  .dc vgs vgsi vgsf vgsdelta vds vds vdsi vdsf vdsdelta
  .print dc i(X1.d)
  .end"""
  #make an aux copy of hspice file to simulate
  shutil.copyfile(templatepath,wheretosimpath+'idvgaux.sp')
  #make an aux copy of modelcard file to simulate
  shutil.copyfile(modelcardpath,wheretosimpath+'modelcardaux.nmos')

  #update path of model and modelcard
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'pathmodelverilog', modelverilogpath)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'pathmodelcard', '\"modelcardaux.nmos\"')

  #bias update
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsi', str(vgs[0]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsf', str(vgs[-1]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsdelta', str(vgs[1]-vgs[0]))

  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsi', str(vds[0]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsf', str(vds[-1]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsdelta', str(vds[1]-vds[0]))

  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'Lparam', Lparam)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'HFINparam',HFINparam)  
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'TFIN_TOPparam', TFIN_TOPparam)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'TFIN_BASEparam',TFIN_BASEparam)   
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'EOTparam', EOTparam)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'NBODYparam',NBODYparam)  
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'NFINparam', NFINparam)  
  
def preparehspiceidvgGEO1v2(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,Lparam,HFINparam,TFIN_TOPparam,TFIN_BASEparam,EOTparam,NBODYparam,NFINparam,PHIGparam,RSHSparam,RSHDparam):
  """*Sample netlist for BSIM-MG
  *Id-Vg Characteristics for NMOS (T = 27 C)

  .option abstol=1e-6 reltol=1e-6 post ingold
  .temp 27

  .hdl pathmodelverilog
  .include pathmodelcard

  * --- Voltage Sources ---
  vds supply  0 dc=0.05
  vgs gate  0 dc=1
  vbs bulk  0 dc=0
  * --- Transistor ---
  X1 supply gate 0 bulk nmos1 L=Lparam NFIN=NFINparam 
  * --- DC Analysis ---
  .dc vgs vgsi vgsf vgsdelta vds vds vdsi vdsf vdsdelta
  .print dc i(X1.d)
  .end"""
  #make an aux copy of hspice file to simulate
  shutil.copyfile(templatepath,wheretosimpath+'idvgaux.sp')
  #make an aux copy of modelcard file to simulate
  shutil.copyfile(modelcardpath,wheretosimpath+'modelcardaux.nmos')

  #update path of model and modelcard
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'pathmodelverilog', modelverilogpath)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'pathmodelcard', '\"modelcardaux.nmos\"')

  #bias update
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsi', str(vgs[0]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsf', str(vgs[-1]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsdelta', str(vgs[1]-vgs[0]))

  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsi', str(vds[0]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsf', str(vds[-1]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsdelta', str(vds[1]-vds[0]))

  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'Lparam', Lparam)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'HFINparam',HFINparam)  
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'TFIN_TOPparam', TFIN_TOPparam)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'TFIN_BASEparam',TFIN_BASEparam)   
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'EOTparam', EOTparam)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'NBODYparam',NBODYparam)  
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'NFINparam', NFINparam)  
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'PHIGparam', PHIGparam)  
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'RSHSparam', RSHSparam)  
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'RSHDparam', RSHDparam)  
  
def preparehspiceidvgGEO4(wheretosimpath,templatepath,modelverilogpath,modelcardpath,vgs,vds,Lparam,Ach_UFCMparam,Cins_UFCMparam,W_UFCMparam,NBODYparam,NFINparam):
#L=Lparam Ach_UFCM=Ach_UFCMparam Cins_UFCM=Cins_UFCMparam W_UFCM=W_UFCMparam NBODY=NBODYparam NFIN=NFINparam 
  """*Sample netlist for BSIM-MG
  *Id-Vg Characteristics for NMOS (T = 27 C)

  .option abstol=1e-6 reltol=1e-6 post ingold
  .temp 27

  .hdl pathmodelverilog
  .include pathmodelcard

  * --- Voltage Sources ---
  vds supply  0 dc=0.05
  vgs gate  0 dc=1
  vbs bulk  0 dc=0
  * --- Transistor ---
  X1 supply gate 0 bulk nmos1 L=Lparam NFIN=NFINparam 
  * --- DC Analysis ---
  .dc vgs vgsi vgsf vgsdelta vds vds vdsi vdsf vdsdelta
  .print dc i(X1.d)
  .end"""
  #make an aux copy of hspice file to simulate
  shutil.copyfile(templatepath,wheretosimpath+'idvgaux.sp')
  #make an aux copy of modelcard file to simulate
  shutil.copyfile(modelcardpath,wheretosimpath+'modelcardaux.nmos')

  #update path of model and modelcard
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'pathmodelverilog', modelverilogpath)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'pathmodelcard', '\"modelcardaux.nmos\"')

  #bias update
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsi', str(vgs[0]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsf', str(vgs[-1]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vgsdelta', str(vgs[1]-vgs[0]))

  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsi', str(vds[0]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsf', str(vds[-1]))
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'vdsdelta', str(vds[1]-vds[0]))

  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'Lparam', Lparam)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'Ach_UFCMparam',Ach_UFCMparam)  
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'Cins_UFCMparam', Cins_UFCMparam)
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'W_UFCMparam',W_UFCMparam)   
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'NBODYparam',NBODYparam)  
  sf.inplace_change(wheretosimpath+'idvgaux.sp', 'NFINparam', NFINparam)    
