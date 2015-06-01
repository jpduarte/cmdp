#fit support using lsf for compact models
#Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley
import os, sys
import time
import imp # to import class using path directly
import numpy as np
from numpy import loadtxt
import scipy.optimize as optimization
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
    
def fungen1(self,device):
  #create function of x input, and 1 variable to optimize "a"
  def functofit(x,a):
    #update variable to optimize
    for paramtofit in self.paramtofit:
      device.updateparameter(paramtofit,a)   
    #evaluate x 2D array, it contains bias and parameters      
    rows = len(x[:,0])
    i=0
    y = []
    while (i<rows):
      #get biases
      bias = x[i,0:len(self.nodes)]
      j=0
      for deviceparameter in self.deviceparameter:
        device.updateparameter(deviceparameter,x[i,len(self.nodes)+j])       
        j+=1
      valuesvar,namesvar = device.analog(*tuple(bias))
      y.append(valuesvar[0])
      i+=1
    return np.array(y)
  return functofit   
   
class fitmodelclass:
  def __init__(self,name):
    self.name = name
    self.version = 'v1'  
    self.devicesymbol = 'X' #add device symbol for Hspice, ex: sim1.updateparameter('devicesymbol','X')
    self.analysis = 'dc' #sim1.updateparameter('analysis','dc')#define type of simulations
    self.modelpath = '' #sim1.updateparameter('modelpath','/users/jpduarte/BSIM_CM_Matlab/BSIM_model_development_v2/DM_Verilog_Hspice/Models_Verilog/BSIMIMG/code/bsimimg.va')
    self.fitcardpath = '' #sim1.updateparameter('modelcardpath','/users/jpduarte/BSIM_CM_Matlab/BSIM_model_development_v2/DM_Verilog_Hspice/Models_Verilog/BSIMIMG/benchmark_tests/modelcard.nmos')#add path to model card of device under study
    self.fitfolder = ''#sim1.updateparameter('simulationfolder',rootfolder+'/cmdp/user1/project1/hspicesimulations/idvg/')#add path to folder which will contain simulation files
    self.fitfilename = 'anyname'#sim1.updateparameter('simfilename','hspicesimaux')#define simulation file name
    self.fitresultfilename = 'outname.txt'#sim1.updateparameter('simresultfilename','hspicesimauxresult.txt')#define simulation final results file 
    self.nodes = ['Vd', 'Vg', 'Vs', 'Vb'] #sim1.updateparameter('nodes',['Vd', 'Vg', 'Vs', 'Vb'])#include node names
    self.deviceparameter = ['L']#sim1.updateparameter('deviceparameter',['L'])#device parameters defined to sweep in simulation     
    self.biasrange = [[0,1],[0,1],[0,1],[0,1]]#indicate range where fit has to be done
    self.deviceparameterrange = [[50e-9,100e-9]]#indicate range where fit has to be done
    self.vartofitdata = ['Id']#variable that will be use as reference for fitting, from data  
    self.vartofitmodel = ['Ids']#variable that will be use as reference for fitting, from model      
    self.alldatafile = 'alldata.txt'
    self.paramtoinclude = ['Ids','Vd']#this is to include parameters in final data file which contains all the data
    self.paramtofit = ['tins']
  def updateparameter(self,name,value):
    #this funtion update a parameter in the model
    exec "self."+name+' = '+'value'    

  def resetdata(self):
    #this function erase file with all data so it can be created again
    if os.path.isfile(self.fitfolder+self.alldatafile): 
      os.remove(self.fitfolder+self.alldatafile)     
 
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
 
    filealldata.close()   
    
  ############################################################################      
  def fitparameters(self):
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
    
    device.returnvar = self.vartofitmodel 
    
    filealldata.close()
    
    functofit = fungen1(self,device)
    outaux = functofit(biasandparam,device.tins)
    print outaux
    sigma = [1]
    x0 = device.tins
    print x0
    print optimization.curve_fit(functofit, biasandparam, ydata, x0, sigma)


