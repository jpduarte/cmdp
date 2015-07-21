#fit support using lsf for compact models
#Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley
import os, sys
import time
import imp # to import class using path directly
import numpy as np
from numpy import loadtxt
import scipy.optimize as optimization
import supportfunctions
import shutil #to copy a file
import fileinput
import functionauxfit

'''
arrange all data input

fit set
bias
param to fit

update
'''

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
   
class fitmodelclass:
  def __init__(self,name):
    self.name = name
    self.version = 'v1'  
    self.devicesymbol = 'X' #add device symbol for Hspice, ex: sim1.updateparameter('devicesymbol','X')
    self.analysis = 'dc' #sim1.updateparameter('analysis','dc')#define type of simulations
    self.modelpath = '' #sim1.updateparameter('modelpath','/users/jpduarte/BSIM_CM_Matlab/BSIM_model_development_v2/DM_Verilog_Hspice/Models_Verilog/BSIMIMG/code/bsimimg.va')
    self.modelcardpath = '' #sim1.updateparameter('modelcardpath','/users/jpduarte/BSIM_CM_Matlab/BSIM_model_development_v2/DM_Verilog_Hspice/Models_Verilog/BSIMIMG/benchmark_tests/modelcard.nmos')#add path to model card of device under study
    self.modelcardpathfinal = ''
    self.fitfolder = ''#sim1.updateparameter('simulationfolder',rootfolder+'/cmdp/user1/project1/hspicesimulations/idvg/')#add path to folder which will contain simulation files
    self.fitfilename = 'anyname'#sim1.updateparameter('simfilename','hspicesimaux')#define simulation file name
    self.simulationresultsfilename = 'outname.txt'#sim1.updateparameter('simresultfilename','hspicesimauxresult.txt')#define simulation final results file 
    self.nodes = ['Vd', 'Vg', 'Vs', 'Vb'] #sim1.updateparameter('nodes',['Vd', 'Vg', 'Vs', 'Vb'])#include node names
    self.deviceparameter = ['L']#sim1.updateparameter('deviceparameter',['L'])#device parameters defined to sweep in simulation     
    self.biasrange = [[0,1],[0,1],[0,1],[0,1]]#indicate range where fit has to be done
    self.deviceparameterrange = [[50e-9,100e-9]]#indicate range where fit has to be done
    self.vartofitdata = ['Id']#variable that will be use as reference for fitting, from data  
    self.vartofitmodel = ['Ids']#variable that will be use as reference for fitting, from model      
    self.vartosave = ['Ids']
    self.alldatafile = 'alldata.txt'
    self.paramtoinclude = ['Ids','Vd']#this is to include parameters in final data file which contains all the data
    self.paramtofit = ['tins']#parameters that are being fit
    self.inputfileformat = ''
  ############################################################################
  def updateparameter(self,name,value):
    #this funtion update a parameter in the model
    exec "self."+name+' = '+'value'    
  ############################################################################
  def resetdata(self):
    #this function erase file with all data so it can be created again
    if os.path.isfile(self.fitfolder+self.alldatafile): 
      os.remove(self.fitfolder+self.alldatafile)     
  ############################################################################
  def adddata(self,datapathandfile,extraparmnames,extraparmvalues):
    #create folder if does not exist to save files for fitting process
    if not os.path.exists(self.fitfolder):
      os.makedirs(self.fitfolder)  
     
    #check if file with all data for fitting exist, if it does not exist it creates one, if not it will append all new added data
    if os.path.isfile(self.fitfolder+self.alldatafile): 
      #print "file with all data will be appended"
      filealldata = open(self.fitfolder+self.alldatafile, 'a') 
    else:
      #print "file with all data will be created"    
      filealldata = open(self.fitfolder+self.alldatafile, 'w') 
      #analysis set up
      stringanalysis = ''
      for paramtoinclude in self.paramtoinclude:
        stringanalysis = stringanalysis + paramtoinclude + ' '    
      filealldata.write(stringanalysis+' \n')        
    
    target = open( datapathandfile, 'r')
    if (self.inputfileformat == ''):
      header = str.split(target.readline())
      
      #check the order of how parameters must be printed in final file  
      indexinparam = []
      for paramtoinclude in header:
        indexinparam.append( self.paramtoinclude.index(paramtoinclude))     
      for paramtoinclude in extraparmnames:
        indexinparam.append( self.paramtoinclude.index(paramtoinclude))       
      
      for line in target:
        listvalues = str.split(line)
        listvalues = listvalues+extraparmvalues     
        valuessorted = [x for y, x in sorted(zip(indexinparam,listvalues))]
        stringtoprint = ''
        for valueprint in valuessorted:
          stringtoprint+=valueprint + ' '
        filealldata.write(stringtoprint+'\n')
 
    if (self.inputfileformat == 'asifdata'):
      enterdata = 0
      
      for line in target:
        if 'DataName' in line:
          header = [x.strip() for x in line.split(',')]
          header = header[1:]
          
          #check the order of how parameters must be printed in final file  
          indexinparam = []
          for paramtoinclude in header:
            indexinparam.append( self.paramtoinclude.index(paramtoinclude))   

          for paramtoinclude in extraparmnames:
            indexinparam.append( self.paramtoinclude.index(paramtoinclude))           
          break

      for line in target:
        listvalues = [x.strip() for x in line.split(',')] 
        listvalues = listvalues[1:-1]
        listvalues = listvalues+extraparmvalues     
        valuessorted = [x for y, x in sorted(zip(indexinparam,listvalues))]
        stringtoprint = ''
        for valueprint in valuessorted:
          stringtoprint+=valueprint + ' '
        filealldata.write(stringtoprint+'\n') 
 
    filealldata.close() 
  def addalldatainfolder(self,folder,name_parameters_missing, value_parameters_missing):
    data_files = [(x[0], x[2]) for x in os.walk(folder)]
    for datafile in data_files[0][1]:
      print "\n Adding: "+datafile
      self.adddata(folder+datafile, name_parameters_missing, value_parameters_missing) 
  ############################################################################     
  ############################################################################     
  ############################################################################      
  def fitparameters(self):
    print "Starting fit preparation"
    filealldata = open(self.fitfolder+self.alldatafile, 'r')
    header = str.split(filealldata.readline())
    #get index for parameters and biases
    #print header
    indexnode= []
    for nodename in self.nodes:
      indexnode.append( header.index(nodename))  
    #print  indexnode
       
    indexparam= []
    for deviceparameter in self.deviceparameter:
      indexparam.append( header.index(deviceparameter))  
    #print  indexparam   
    
    indexy= []
    for vartofitdata in self.vartofitdata:
      indexy.append( header.index(vartofitdata))  
    #print  indexy       
    
    #check what biases and parameters fullfill fit requerements, biasparam and ydata is created 
    datalist = loadtxt(filealldata,skiprows = 1)
    filealldata.close()
    i=0
    k=0
    biasparam = []
    ydata = []
    for datalistaux in datalist[:,0]:
      j = 0
      flag = 0
      biasparamaux = []
      for index in indexnode:
        maxaux = max(self.biasrange[j])
        minaux = min(self.biasrange[j])
        #print maxaux,minaux,datalist[i,index]
        if ((datalist[i,index]>minaux) and (datalist[i,index]<maxaux)):
          flag+= 1
          biasparamaux.append(datalist[i,index])
        j+=1
        
      j = 0
      for index in indexparam:
        maxaux = max(self.deviceparameterrange[j])
        minaux = min(self.deviceparameterrange[j])
        #print maxaux,minaux,datalist[i,index]
        if ((datalist[i,index]>minaux) and (datalist[i,index]<maxaux)):
          flag+= 1
          biasparamaux.append(datalist[i,index])
        j+=1   
        
      if flag==(len(indexnode)+len(indexparam)):
        biasparam.append(biasparamaux)
        ydata.append(datalist[i,indexy[0]])
      i+=1   
    
    biasandparam = np.array(biasparam)
    ydata = np.array(ydata)
    #print ydata
    
    #import class with path to model & create device for simulation
    CURRENT_DIR = os.path.dirname(os.path.abspath(self.modelpath))
    sys.path.insert(0, CURRENT_DIR)
    deviceaux = imp.load_source('model', self.modelpath)
    device =  deviceaux.compactmodel('X1')

    #update paramter from modelcard
    #TODO: make this a function outside, because many parts use this modelcard to parameter setting
    modelcard = open(self.modelcardpath, 'r') 
    for line in modelcard:
      if (line.find('*')!=0):
        if (line.find('+')>-1):
          line = line.replace('+', '')
          line = line.replace('=', ' ')
          #TODO replace all n,u,p,f, units by number so users can write it in that way as well
          stringinline = str.split(line)
          valueparam = stringinline[-1]
          if (isfloat(valueparam)):
            stringtoexec = 'device.'+stringinline[0]+ '='+valueparam
          else:
            stringtoexec = 'device.'+stringinline[0] +'='+'device.'+valueparam        
          exec stringtoexec     
    modelcard.close()
    device.returnvar = self.vartofitmodel 
    
    #create parameterinitial initial guess, parameterinitial
    parameterinitial = []
    for paramtofit in self.paramtofit:
      exec 'parameterinitial.append(device.'+paramtofit+')' 
    parameterinitial = np.array(parameterinitial)
    #create function to fit based on parameters to fit
    functofit = functionauxfit.fungen(self,device)
    
    #curve fitting
    print "model fitting using optimization.curve_fit"
    parametersfit = optimization.curve_fit(functofit, biasandparam, ydata, parameterinitial)

    #update model card  
    modelcardfinal = open(self.modelcardpathfinal, 'w') 
    modelcard = open(self.modelcardpath, 'r') 
    for line in modelcard:
      flag=0
      if (line[0]=='+'):
        i=0
        for paramtofit in self.paramtofit:
          if paramtofit in line:
            flag=1
            print "updating modelcard with paramter " + paramtofit + ' = ' +str(parametersfit[0][i])
            modelcardfinal.write('+'+paramtofit+' = ' + str(parametersfit[0][i])+'\n')
          i+=1      
      if flag==0:
        modelcardfinal.write(line)
    modelcard.close()
    modelcardfinal.close()
  ############################################################################       
  ############################################################################ 
  ############################################################################ 
  def runsim(self,modelcardpath):
    print "simulation of model"
    #check if simulation path exist, if not, it create the folder
    if not os.path.exists(self.fitfolder):
      os.makedirs(self.fitfolder)
   
    #import class with path to model & create device for simulation
    CURRENT_DIR = os.path.dirname(os.path.abspath(self.modelpath))
    sys.path.insert(0, CURRENT_DIR)
    deviceaux = imp.load_source('model', self.modelpath)
    device =  deviceaux.compactmodel('X1')

    #update paramter from modelcard
    modelcard = open(modelcardpath, 'r') 
    for line in modelcard:
      if (line[0]=='+'):
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
    #print device.returnvar
    
    #################################simulation, device evaluation
    fileresult = open(self.fitfolder+self.simulationresultsfilename, 'w') 
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

