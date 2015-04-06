#hspicepython: it is a class to wrap hspice from python lenguage
#author: Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley
import os
import time
import supportfunctions

import numpy as np

class hspicepython:
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
    self.vartosave = ['Ids']#sim1.updateparameter('vartosave',['Ids','qs']) #add variables to save   
  def updateparameter(self,name,value):
    #this funtion update a parameter in the model
    exec "self."+name+' = '+'value'    
  def runsim(self):
    #check if simulation path exist, if not, it create the folder
    if not os.path.exists(self.simulationfolder):
      os.makedirs(self.simulationfolder)
    #create simulation file
    simfile = open(self.simulationfolder+self.simfilename+'.sp', 'w')
    simfile.write('*script to generate hspice simulation using cmdp, Juan Duarte \n')
    simfile.write('*Date: '+ time.strftime("%m/%d/%Y")+ ', time: ' + time.strftime("%H:%M:%S") + '\n\n')
    #TODO: do not hard code the following 3 lines, use them as parameters?
    simfile.write('.option abstol=1e-6 reltol=1e-6 post ingold \n')
    simfile.write('.option measform=1 \n')
    simfile.write('.temp 27 \n')
    simfile.write('\n')
    #add model path and modelcard path
    simfile.write('.hdl "'+self.modelpath+'" \n')
    simfile.write('.include "'+self.modelcardpath+'" \n')
    simfile.write('\n')
    
    #bias and parameter set up
    for Vbias in self.nodes:
      simfile.write('.PARAM '+Vbias+'_value = 0 \n')
    count = 0
    for parameteraux in self.deviceparameter:
      valueaux= str(self.deviceparametervalue[count][0])
      simfile.write('.PARAM '+parameteraux+'_value = '+valueaux+' \n')
      count+=1
    simfile.write('\n')
    
    #voltage sources set up
    for Vbias in self.nodes:
      simfile.write(Vbias+' '+Vbias+' 0.0 dc = ' + Vbias+'_value \n')    
    simfile.write('\n')

    #device set up: X1 Vd Vg Vs Vb nmos L = 'L_value' TODO: gate length always there?
    stringdevice = self.devicesymbol+'1'
    for Vbias in self.nodes:
      stringdevice = stringdevice+' ' +Vbias    
    
    #obtain devicetype from modelcard  
    modelcardfile = open(self.modelcardpath, 'r') 
    for line in modelcardfile:
      if line.find(".model")==0:
        print line.find(".model")
        stringinline = str.split(line)
        self.devicetype=stringinline[1]
    modelcardfile.close()    
        
    stringdevice = stringdevice + ' '+self.devicetype+ ' '+ self.deviceparameter[0]+' = \''+self.deviceparameter[0]+'_value\''+'\n'
    simfile.write(stringdevice)
    simfile.write('\n')
        
    #analysis set up
    stringanalysis = '.DATA datadc'
    for Vbias in self.nodes:
      stringanalysis = stringanalysis+' ' +Vbias+'_value'
    for parameteraux in self.deviceparameter:
      stringanalysis = stringanalysis+' ' +parameteraux +'_value'
    simfile.write(stringanalysis+' \n')       
    
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

    #write to file combination of all bias and parameter values
    i=0
    column=len(allvaldc)
    rown=len(allvaldc[0])
    while (i<rown):
      stringtowrite = ''
      j=0
      while (j<column):
        stringtowrite = stringtowrite+str(allvaldc[j][i])+' '
        j+=1
      simfile.write(stringtowrite+'\n') 
      i+=1 
    #finish writing dc analysis   
    simfile.write('.ENDDATA \n.dc sweep DATA = datadc \n')  

    
    #print variable, .print dc X1:Dmob X1:U0 X1:Eeffm
    stringtowrite = '.print dc '
    for varname in self.vartosave:
      stringtowrite = stringtowrite + self.devicesymbol+'1:' +varname+' ' 
    simfile.write(stringtowrite)
    simfile.write('\n.end')
    
    simfile.close()   
    #hspice run file TODO: how to check simulation was aborted?
    os.system('hspice ' + self.simulationfolder +self.simfilename+'.sp -o ' + self.simulationfolder+self.simfilename)
   
    #parse results
    outputfiletoread = open(self.simulationfolder+self.simfilename+'.lis', 'r') 
    #state machine
    #00: searching for "x1:"
    #01: searching for x1.d and parsing parameters
    #02: search y, meanhwile parsing data id-vg
    state=0

    #TODO: parse case with many output variables, it include more than one read
    countruns=0
    for line in outputfiletoread:
      #print line
      if (state==0):
        if (line.find("x1:")>-1):
          state=1
          hspicefileresult = open(self.simulationfolder+self.simresultfilename, 'w') 
          stringoutput = ''
          for Vbias in self.nodes:
            stringoutput = stringoutput +Vbias+' '
          for parameteraux in self.deviceparameter:
            stringoutput = stringoutput+parameteraux+' '
          for varname in self.vartosave:
            stringoutput = stringoutput + varname+' '   
          hspicefileresult.write(stringoutput+' \n') 
          i=0
      elif (state==1):
        if line.find("y")>-1 :
          state=0
          hspicefileresult.close()
        else:     
          column=len(allvaldc)
          stringtowrite = ''
          j=0
          while (j<column):
            stringtowrite = stringtowrite+str(allvaldc[j][i])+' '
            j+=1
          #add results  
          stringinline = str.split(line)
          resultcount = 0
          for result in stringinline:
            if resultcount>0:
              stringtowrite = stringtowrite + result+' '
            resultcount+=1 
          hspicefileresult.write(stringtowrite+'\n') 
          i+=1        
