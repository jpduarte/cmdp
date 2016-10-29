#hspicepython: it is a class to wrap hspice from python lenguage
#author: Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley
import os
import time
import supportfunctions

import numpy as np

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

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
    self.abstol = '1e-6' #Sets the absolute error tolerance for branch currents in DC and transient analysis.
    self.reltol = '1e-6' #Sets the relative error tolerance for voltages from iteration to iteration
    self.absv = '1e-6' #Sets the absolute minimum voltage for DC and transient analysis.
  def updateparameter(self,name,value):
    #this funtion update a parameter in the model
    exec "self."+name+' = '+'value'    
  def runsim(self,flagrun=1):
    #check if simulation path exist, if not, it create the folder
    if not os.path.exists(self.simulationfolder):
      os.makedirs(self.simulationfolder)
    #create simulation file
    simfile = open(self.simulationfolder+self.simfilename+'.sp', 'w')
    simfile.write('*script to generate hspice simulation using cmdp, Juan Duarte \n')
    simfile.write('*Date: '+ time.strftime("%m/%d/%Y")+ ', time: ' + time.strftime("%H:%M:%S") + '\n\n')
    #TODO: do not hard code the following 3 lines, use them as parameters?
    simfile.write('.option abstol='+self.abstol +' reltol=' +self.reltol +' post ingold \n')
    simfile.write('.option ABSV='+self.absv + ' \n')
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
    
    #add name of device like nmos1 based on modelcard
    stringdevice = stringdevice + ' '+self.devicetype
    #add all the parameters to be sweeped
    parameterstringaux = ''
    for parameteraux in self.deviceparameter:    
      parameterstringaux = parameterstringaux + ' '+ parameteraux+' = \''+parameteraux+'_value\''
    #put all together, ex: X1 Vd Vg Vs Vb nmos1 L = 'L_value'
    stringdevice = stringdevice + parameterstringaux+'\n'  
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
    if (flagrun==1):
      os.system('hspice ' + self.simulationfolder +self.simfilename+'.sp -o ' + self.simulationfolder+self.simfilename)
    #print(allvaldc)
    #parse results
    #self.hspicetotex('x','y',allvaldc)
    self.allvaldc = allvaldc
   
    '''#parse results
    outputfiletoread = open(self.simulationfolder+self.simfilename+'.lis', 'r') 
    #state machine
    #00: searching for "x1:"
    #01: searching for x1.d and parsing parameters
    #02: search y, meanhwile parsing data id-vg
    state=0

    #TODO: parse case with many output variables, it include more than one read
    countruns=0
    print "print results?"
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
          i+=1'''
          
          
  def hspicetotex(self,stringstart,stringend):#,allvaldc):
    allvaldc = self.allvaldc 
    outputfiletoread = open(self.simulationfolder+self.simfilename+'.lis', 'r') 
    state=0
    #this take data from .lis file and save information in python dictionaries
    countruns=0
    totalvariables = []
    namesaux = ''
    names_vol_curr = []
    valuestoprint = {}
    flagfirstread = 1 #use to save only once variable that is being sweeped for simulation, can be input voltage or time
    for line in outputfiletoread:
      line = line.replace("x1:", "") #this replaces all string x1: for empty, so only variable name is saved
      if (state==0):
        #stays in this state until stringstart is found
        if (line.find(stringstart)==0):#TODO: make this more robust this work only for hspice 2015, 0 for 2012
          state=1
          
      elif (state==1):
        #this states is to add variables names and first values
        stringinline = str.split(line)
        if len(stringinline)>0:
        
          if (is_number(stringinline[0])):
            if flagfirstread==1:
              namesaux = ['sweepvar']+namesaux
            
            totalvariables= totalvariables + namesaux
              
            if (len(vol_curr_string)<(len(namesaux)+1+flagfirstread*-1)):
              vol_curr_string = vol_curr_string + ['unknown']
              
            if flagfirstread==0:
              vol_curr_string = vol_curr_string[1:]
            names_vol_curr = names_vol_curr + vol_curr_string
            
            state=2

            count=0
            if flagfirstread==1:
              addint = 0
            else:
              addint = 1 
            #this is to create cell for given printed variables
            for variable in namesaux:
              valuestoprint[variable] = []
              count+=1
              
            count=0
            for variable in namesaux:
              valuestoprint[variable].append(stringinline[count+addint])
              count+=1              
          else:
            #this saves previous string so names can be saved
            vol_curr_string = namesaux
            namesaux = stringinline
            

      elif (state==2):
        #this state keep adding values to the variables until end string is found then it goes back to state 0
        if (line.find(stringend)==0) : #TODO: make this more robust this work only for hspice 2015, 0 for 2012
          flagfirstread = 0

          state=0
        else:   
          stringinline = str.split(line)
          count=0
          if flagfirstread==1:
            addint = 0
          else:
            addint = 1        
          for variable in namesaux:
            valuestoprint[variable].append(stringinline[count+addint])
            count+=1    
       

    #write all data to text file   
    hspicefileresult = open(self.simulationfolder+self.simresultfilename, 'w')        
    namesfinal = valuestoprint.keys()
    stringtowrite = ''    
    count=0  
    
    for Vbias in self.nodes:
      stringtowrite = stringtowrite +Vbias.lower()+' '
    for parameteraux in self.deviceparameter:
      stringtowrite = stringtowrite+parameteraux.lower()+' '    
    for variable in namesfinal:
      stringtowrite = stringtowrite + variable+ ' '
      count+=1
      
    hspicefileresult.write(stringtowrite[:-1]+'\n')

    lenarrayaux = len(valuestoprint[namesfinal[0]])
    count=0
    while (count<lenarrayaux):
      stringtowrite = '' 
      column=len(allvaldc)

      j=0
      while (j<column):
        stringtowrite = stringtowrite+str(allvaldc[j][count])+' '
        j+=1      
      
      for namevar in namesfinal:
        stringtowrite = stringtowrite + valuestoprint[namevar][count] + ' '
      hspicefileresult.write(stringtowrite[:-1]+'\n')
      count+=1
    hspicefileresult.close() 
  ##########################################################################################################################  
  def hspicetotex2(self,stringstart,stringend, simfilename,simresultfilename):
    #allvaldc = self.allvaldc 
    outputfiletoread = open(simfilename, 'r') 
    state=0
    #this take data from .lis file and save information in python dictionaries
    countruns=0
    totalvariables = []
    namesaux = ''
    names_vol_curr = []
    valuestoprint = {}
    flagfirstread = 1 #use to save only once variable that is being sweeped for simulation, can be input voltage or time
    for line in outputfiletoread:
      line = line.replace("x1:", "") #this replaces all string x1: for empty, so only variable name is saved
      #print (line)
      if (state==0):
        #stays in this state until stringstart is found
        if (line.find(stringstart)==0):#TODO: make this more robust this work only for hspice 2015, 0 for 2012
          state=1
          
      elif (state==1):
        #this states is to add variables names and first values
        #print (line)
        stringinline = str.split(line)
        if len(stringinline)>0:
        
          if (is_number(stringinline[0])):
            if flagfirstread==1:
              namesaux = ['sweepvar']+namesaux
            
            totalvariables= totalvariables + namesaux
              
            if (len(vol_curr_string)<(len(namesaux)+1+flagfirstread*-1)):
              vol_curr_string = vol_curr_string + ['unknown']
              
            if flagfirstread==0:
              vol_curr_string = vol_curr_string[1:]
            names_vol_curr = names_vol_curr + vol_curr_string
            
            state=2

            count=0
            if flagfirstread==1:
              addint = 0
            else:
              addint = 1 
            #this is to create cell for given printed variables
            for variable in namesaux:
              valuestoprint[variable] = []
              count+=1
              
            count=0
            for variable in namesaux:
              valuestoprint[variable].append(stringinline[count+addint])
              count+=1              
          else:
            #this saves previous string so names can be saved
            vol_curr_string = namesaux
            namesaux = stringinline
            

      elif (state==2):
        #this state keep adding values to the variables until end string is found then it goes back to state 0
        if (line.find(stringend)==0) : #TODO: make this more robust this work only for hspice 2015, 0 for 2012
          flagfirstread = 0

          state=0
        else:   
          stringinline = str.split(line)
          count=0
          if flagfirstread==1:
            addint = 0
          else:
            addint = 1        
          for variable in namesaux:
            valuestoprint[variable].append(stringinline[count+addint])
            count+=1    
       

    #write all data to text file   
    #hspicefileresult = open(self.simulationfolder+self.simresultfilename, 'w')    
    hspicefileresult = open(simresultfilename, 'w')      
    namesfinal = valuestoprint.keys()
    print (namesfinal)
    stringtowrite = ''    
    count=0  
    
    '''for Vbias in self.nodes:
      stringtowrite = stringtowrite +Vbias.lower()+' '
    for parameteraux in self.deviceparameter:
      stringtowrite = stringtowrite+parameteraux.lower()+' '  '''  
    for variable in namesfinal:
      stringtowrite = stringtowrite + variable+ ' '
      count+=1
      
    hspicefileresult.write(stringtowrite[:-1]+'\n')

    lenarrayaux = len(valuestoprint[namesfinal[0]])
    count=0
    while (count<lenarrayaux):
      stringtowrite = '' 
      '''column=len(allvaldc)

      j=0
      while (j<column):
        stringtowrite = stringtowrite+str(allvaldc[j][count])+' '
        j+=1'''     
      
      for namevar in namesfinal:
        stringtowrite = stringtowrite + valuestoprint[namevar][count] + ' '
      hspicefileresult.write(stringtowrite[:-1]+'\n')
      count+=1
    hspicefileresult.close()    
                        
