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
      else:
        stringtoexec = 'device.'+self.deviceparameter[j-len(self.nodes)]+ '='+str(allvaldc[j][i])
        exec stringtoexec
      j+=1  
    valuesvar,namesvar = device.analog(*tuple(bias))#model evaluation
    resultsimstring = ' '.join(map(str, valuesvar)) 
    fileresult.write(stringtowrite+resultsimstring+'\n') 
    i+=1 
  fileresult.close()
