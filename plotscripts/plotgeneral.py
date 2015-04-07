#plot general class
#author: Juan Pablo Duarte, jpduarte@berkeley.edu
#BSIM Group, UC Berkeley

import supportfunctions as sf
import matplotlib.pyplot as plt
from numpy import loadtxt
import numpy as np

from numpy.matlib import repmat
from scipy.misc import factorial
from scipy import sparse
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm

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

def K_generator(x,order):
#this return matrix to find the derivative, x is the variable to be derived and order is the derivative order
  N=len(x);
  K = lil_matrix((N, N))
  K[0,:6]=mkfdstencil(x[0:6],x[0],order)
  K[1,:6]=mkfdstencil(x[0:6],x[1],order)
  K[2,:6]=mkfdstencil(x[0:6],x[2],order)

  i=3
  for xbar in x[3:-3]:
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
#arrange X and Y matrices
def rearrangearray(arrayXa,elementpercylce,numberelement):
#this function reshpae array to be printing 
  arrayXb = arrayXa.reshape((elementpercylce, len(arrayXa)/elementpercylce))
  arrayXc = np.transpose(arrayXb)
  arrayXd = arrayXc.reshape((len(arrayXa)/numberelement,numberelement))
  arrayXe = np.transpose(arrayXd)
  return arrayXe

def findelementpercylce(arrayaux):
#this function find the number of times a element is keep constinously
  firstelement = arrayaux[0]
  lengthaux = len(arrayaux)
  flag=1
  i=1

  elementpercylce = 0
  while (i<(lengthaux+1) & flag):
    elementpercylce = i-1
    print elementpercylce
    if (abs(arrayaux(i)-firstelement)>0): #%TODO: check abs condition
      flag=0
    i=i+1
  
  if (flag):
    elementpercylce = elementpercylce+1
  return elementpercylce

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

      xarray = datalist[:,xindex]
      yarray = datalist[:,yindex]
      
      #this is to identify index how to re-shape matrix for right plotting
      numberelement = 0
      numberelementaux = len(np.unique(xarray))  
      numberelementmaxpossible = len(xarray)
      if( (np.mod(len(xarray),numberelementaux)==0) & ((numberelementmaxpossible-numberelementaux)>0) ):
        numberelement = numberelementaux;
        elementpercylce = findelementpercylce(xarray)*numberelement
      
      if (numberelement==0):
        numberelement = numberelementmaxpossible
        elementpercylce = numberelement
        
      #reshape matrix to plot lines
      xarray = rearrangearray(xarray,elementpercylce,numberelement)
      yarray = rearrangearray(yarray,elementpercylce,numberelement)
      
      #plot
      plt.figure(fignumber)
      
      #plot variable or its derivatives: TODO: it plot derivate with respect to x-axis, update derivative with respect to any variable
      if (self.derivativeorder<1):      
        plt.plot( xarray, yarray, self.symbol, lw=self.lw, color=self.color  )
      else :
        K = K_generator(xarray[:,0],self.derivativeorder) 
        plt.plot( xarray, K*yarray, self.symbol, lw=self.lw, color=self.color  )
      
      #log scale check  
      if self.ylogflag==1:
        ax = plt.gca()
        ax.set_yscale('log') 
    
    #x and y axis label          
    target.close()
    ax = plt.gca()
    ax.set_xlabel(xstring)
    if (self.derivativeorder<1):
      ax.set_ylabel(ystring)  
    else:
      ax.set_ylabel('d^'+str(self.derivativeorder)+' '+ystring+'/d'+xstring+'^'+str(self.derivativeorder))


