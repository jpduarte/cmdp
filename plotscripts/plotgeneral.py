#plot general class
#author: Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

import supportfunctions as sf
import matplotlib.pyplot as plt
from numpy import loadtxt


def guess_seq_len(seq):
    guess = 1
    max_len = len(seq) / 2
    for x in range(2, max_len):
        if sum(abs(seq[0:x] - seq[x:2*x]))<1e-10 :
            return x
    return guess

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
    plt.figure(fignumber)
    plt.plot(biaslist[xvariable-1],fval,self.symbol, lw=self.lw, color=self.color )#markerfacecolor=self.markerfacecolor,
    if self.ylogflag==1:
      ax = plt.gca()
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
      plt.figure(fignumber)      
      plt.plot( datalist[:,xindex], datalist[:,yindex], self.symbol, lw=self.lw, color=self.color  )
      if self.ylogflag==1:
        ax = plt.gca()
        ax.set_yscale('log')       
    target.close()
    ax = plt.gca()
    ax.set_xlabel(xstring)
    ax.set_ylabel(ystring)  
    #print guess_seq_len(datalist[:,xindex])  

