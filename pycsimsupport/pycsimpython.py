#pycsimpython: it is a wrapper for pycsim simulator
#author: Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley
import os, sys
import time
import supportfunctions
import imp # to import class using path directly

import numpy as np

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

class pycsimpython:
  def __init__(self,name):
    self.name = name
    self.version = 'v1'  
    self.devicesymbol = 'X' #add device symbol for Hspice, ex: sim1.updateparameter('devicesymbol','X')
    self.analysis = 'dc' #sim1.updateparameter('analysis','dc')#define type of simulations
    self.modelpath = '' #sim1.updateparameter('modelpath','/users/jpduarte/BSIM_CM_Matlab/BSIM_model_development_v2/DM_Verilog_Hspice/Models_Verilog/BSIMIMG/code/bsimimg.va')
    self.modelcardpath = '' #sim1.updateparameter('modelcardpath','/users/jpduarte/BSIM_CM_Matlab/BSIM_model_development_v2/DM_Verilog_Hspice/Models_Verilog/BSIMIMG/benchmark_tests/modelcard.nmos')#add path to model card of device under study
    self.simulationfolder = ''#sim1.updateparameter('simulationfolder',rootfolder+'/cmdp/user1/project1/hspicesimulations/idvg/')#add path to folder which will contain simulation files
    self.simfilename = 'anyname'#sim1.updateparameter('simfilename','hspicesimaux')#define simulation file name
    self.simresultfilename = 'outname.txt'#sim1.updateparameter('simresultfilename','hspicesimauxresult.txt')#define simulation final results file 
    self.nodes = ['Vd', 'Vg', 'Vs', 'Vb'] #sim1.updateparameter('nodes',['Vd', 'Vg', 'Vs', 'Vb'])#include node names in the order defined in code
    self.dcbiases = [np.linspace(0.05,1,2), np.linspace(0.0,1,20), [0], [0]] #sim1.updateparameter('dcbiases',[np.linspace(0.05,1,2), np.linspace(0.0,1,20), [0], [0]])#values for bias conditions of nodes
    self.deviceparameter = ['L']#sim1.updateparameter('deviceparameter',['L'])#device parameters defined to sweep in simulation
    self.deviceparametervalue = [[1000e-9]]#sim1.updateparameter('deviceparametervalue',[[1000e-9]])#device parameter values for simulation
    self.vartosafe = ['Ids']#sim1.updateparameter('vartosafe',['Ids','qs']) #add variables to save   
  def updateparameter(self,name,value):
    #this funtion update a parameter in the model
    exec "self."+name+' = '+'value'    
  def runsim(self):
    #check if simulation path exist, if not, it create the folder
    if not os.path.exists(self.simulationfolder):
      os.makedirs(self.simulationfolder)
 
    #import class with path to model & create device for simulation
    CURRENT_DIR = os.path.dirname(os.path.abspath(self.modelpath))
    sys.path.insert(0, CURRENT_DIR)
    deviceaux = imp.load_source('model', self.modelpath)
    device =  deviceaux.compactmodel('X1')

    #update paramter from modelcard
    modelcard = open(self.modelcardpath, 'r') 
    for line in modelcard:
      if (line.find('+')>-1):
        line = line.replace('+', '')
        line = line.replace('=', ' ')
        #TODO replace all n,u,p,f, etc
        stringinline = str.split(line)
        valueparam = stringinline[-1]
        if (isfloat(valueparam)):
          stringtoexec = 'device.'+stringinline[0]+ '='+valueparam
        else:
          stringtoexec = 'device.'+stringinline[0] +'='+'device.'+valueparam        
        exec stringtoexec     
    modelcard.close()
    
    device.returnvar = self.vartosave
    
    #################################simulation, device evaluation
    fileresult = open(self.simulationfolder+self.simresultfilename, 'w') 
    #analysis set up
    stringanalysis = ''
    for Vbias in self.nodes:
      stringanalysis = stringanalysis + Vbias + ' '
    for parameteraux in self.deviceparameter:
      stringanalysis = stringanalysis + parameteraux + ' '
    for vartoprint in self.vartosave:
      stringanalysis = stringanalysis + vartoprint + ' '      
    fileresult.write(stringanalysis+' \n')  
    #print stringanalysis
    
    #generate combination of all bias and parameter values
    count=0
    stringallarraybias =''
    for arraybias in self.dcbiases:
      stringarraybias = 'ar'+str(count) + ' = np.asarray(arraybias)'
      stringallarraybias = stringallarraybias + ' ar'+str(count)+',' 
      exec stringarraybias
      count+=1
    for arraydeviceparameter in self.deviceparametervalue:
      stringarraybias = 'ar'+str(count) + ' = np.asarray(arraydeviceparameter)'
      stringallarraybias = stringallarraybias + ' ar'+str(count)+',' 
      exec stringarraybias
      count+=1          
    stringtoeval = 'allvaldc = supportfunctions.meshgrid2('+stringallarraybias[:-1]+')'
    exec stringtoeval
    
    #this save to text file by evaluating model
    i=0
    column=len(allvaldc)
    rown=len(allvaldc[0])
    while (i<rown):
      stringtowrite = ''
      j=0
      bias = []
      while (j<column):
        stringtowrite = stringtowrite+str(allvaldc[j][i])+' '
        if (j<len(self.nodes)):
          bias.append(allvaldc[j][i])
        j+=1  
      valuesvar,namesvar = device.analog(*tuple(bias))#model evaluation
      resultsimstring = ' '.join(map(str, valuesvar)) 
      fileresult.write(stringtowrite+resultsimstring+'\n') 
      i+=1 

    fileresult.close()
