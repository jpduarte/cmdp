#this returns a function depending on the number of parameters being fitting, max number of parameters is 7
#Juan Duarte, jpduarte@berkeley.edu

import numpy as np

def fungen(self,device):
  #create function of x input, and 1 variable to optimize "a"
  if len(self.paramtofit)==1:
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
  #################################################################################2
  if len(self.paramtofit)==2:
    def functofit(x,a,b):
      #update variable to optimize
      values = [a,b]
      i=0
      for paramtofit in self.paramtofit:
        device.updateparameter(paramtofit,values[i])  
        i+=1 
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
  #################################################################################3
  if len(self.paramtofit)==3:
    def functofit(x,a,b,c):
      #update variable to optimize
      values = [a,b,c]
      i=0
      for paramtofit in self.paramtofit:
        device.updateparameter(paramtofit,values[i])  
        i+=1 
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
  #################################################################################4
  if len(self.paramtofit)==4:
    def functofit(x,a,b,c,d):
      #update variable to optimize
      values = [a,b,c,d]
      i=0
      for paramtofit in self.paramtofit:
        device.updateparameter(paramtofit,values[i])  
        i+=1 
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
  #################################################################################5
  if len(self.paramtofit)==5:
    def functofit(x,a,b,c,d,e):
      #update variable to optimize
      values = [a,b,c,d,e]
      i=0
      for paramtofit in self.paramtofit:
        device.updateparameter(paramtofit,values[i])  
        i+=1 
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
  #################################################################################6
  if len(self.paramtofit)==6:
    def functofit(x,a,b,c,d,e,f):
      #update variable to optimize
      values = [a,b,c,d,e,f]
      i=0
      for paramtofit in self.paramtofit:
        device.updateparameter(paramtofit,values[i])  
        i+=1 
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
  #################################################################################7
  if len(self.paramtofit)==7:
    def functofit(x,a,b,c,d,e,f,g):
      #update variable to optimize
      values = [a,b,c,d,e,f,g]
      i=0
      for paramtofit in self.paramtofit:
        device.updateparameter(paramtofit,values[i])  
        i+=1 
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
  
