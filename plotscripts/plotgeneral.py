#plot general class
#author: Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

import supportfunctions as sf
import pylab
from numpy import loadtxt

class plotgeneral:
  def __init__(self):#, model):
    self.version = 'v1'
    #self.model = model
    
    #defaul parameters
    self.symbol  = 'o'
    self.color = 'k'
    self.markerfacecolor = (1, 1, 1, 1)
    self.lw=1
    self.ylogflag = '0'
    
  def updateparameter(self,name,value):
    #this funtion update a parameter in the model
    if type(value) == type(''):
      exec "self."+name+' = '+'\''+value+'\''
      #print "self."+name+' = '+'\''+value+'\''
    else:
      exec "self."+name+' = '+str(value)    
      #print "self."+name+' = '+str(value) 
      
      
  def plotmodel(self,f,fignumber, xvariable,*args):
    #this function plots the function "f" for all the combinations in the evaluation points *args, then, it plot the evalaution parameter in position xvariable
    biaslist = sf.meshgrid2(*args) #list of ndarray, take all combinations of imputs biases
    totalbiases = len(biaslist[0])#check the length of the final array
    count=0
    fval = []
    while count<totalbiases:
      biasplotaux = []
      for biasx in biaslist:
        biasplotaux.append(biasx[count])
        
      #print biasplotaux
      fval.append( f(*tuple(biasplotaux)))
      count+=1  

    #plot
    pylab.figure(fignumber)
    pylab.plot(biaslist[xvariable-1],fval,self.symbol, lw=self.lw, color=self.color )#markerfacecolor=self.markerfacecolor,
    if self.ylogflag==1:
      ax = pylab.gca()
      ax.set_yscale('log') 
        
    
  def plotfiledata(self,pathandfile,xstring,ystring,fignumber):
    #this function open a file pathandfile and plot the columns with xstring and ystring string header
    target = open( pathandfile, 'r')
    header = str.split(target.readline())
    if len(header)>0:
      xindex =  header.index(xstring)
      yindex =  header.index(ystring)
    
      datalist = loadtxt(pathandfile,skiprows = 1)

      #plot
      pylab.figure(fignumber)      
      pylab.plot( datalist[:,xindex], datalist[:,yindex], self.symbol, lw=self.lw, color=self.color  )
      if self.ylogflag==1:
        ax = pylab.gca()
        ax.set_yscale('log')       
    target.close()

