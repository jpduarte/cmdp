#plot general class
#author: Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

import supportfunctions as sf
import matplotlib.pyplot as plt
from numpy import loadtxt

#########################################################################
#Derivative Support functions
def mkfdstencil(x,xbar,k):
#this funtion is sue to create finite diference method matrix
  maxorder            = len(x)
  h_matrix            = repmat(np.transpose(x)-xbar,maxorder,1)
  powerfactor_matrix  = np.transpose(repmat(np.arange(0,maxorder),maxorder,1))
  factorialindex      = np.transpose(repmat(factorial(np.arange(0,maxorder)),maxorder,1))
  taylormatrix        = h_matrix ** powerfactor_matrix /factorialindex
  derivativeindex     = np.zeros(maxorder)
  derivativeindex[k]  = 1 
  u = np.linalg.solve(taylormatrix,derivativeindex)
  return u

def K_generator(x):
#this return matrix of Poisson Equation in Silicon Fin, with Neuman BC at both ends
  N=len(x);
  K = lil_matrix((N, N))
  order = 1
  K[0,:6]=mkfdstencil(x[0:6],x[0],order)
  K[1,:6]=mkfdstencil(x[0:6],x[1],order)
  K[2,:6]=mkfdstencil(x[0:6],x[2],order)

  i=3
  for xbar in x[3:-4]:
    #print i
    K[i,i-3:i+3]=mkfdstencil(x[i-3:i+3],xbar,order)
    i+=1
  #print i
  K[i,-7:-1]=mkfdstencil(x[-7:-1],x[-3],order)  
  i+=1
  K[i,-7:-1]=mkfdstencil(x[-7:-1],x[-2],order)
  i+=1  
  K[i,-7:-1]=mkfdstencil(x[-7:-1],x[-1],order)  
  return K.tocsr()
#########################################################################

def guess_seq_len(seq):
    guess = 1
    max_len = len(seq) / 2
    for x in range(2, max_len):
        if sum(abs(seq[0:x] - seq[x:2*x]))<1e-10 :
            return x
    return guess

#########################################################################
#plotgeneral class definition
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
    self.derivativeorder = 0
  def updateparameter(self,name,value):
    #this funtion update a parameter in the model
    if type(value) == type(''):
      exec "self."+name+' = '+'\''+value+'\''
      #print "self."+name+' = '+'\''+value+'\''
    else:
      exec "self."+name+' = '+str(value)    
      #print "self."+name+' = '+str(value) 
    
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

