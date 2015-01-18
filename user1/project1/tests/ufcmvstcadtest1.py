#first test for class UFCM
#juan pablo

import sys
sys.path.insert(0, '/users/jpduarte/research/cmdp/compactmodels')
sys.path.insert(0, '/users/jpduarte/research/cmdp/plotscripts')

import UFCM
import plotgeneral
import pylab
import numpy as np

X1 = UFCM.UFCM('x1')
print X1.version
print X1.draincurrent(*(1,1,0,0))

###############plot using plotgeneral
P1 = plotgeneral.plotgeneral()

vg=np.linspace(0,0.86,50)
vd=np.linspace(0.05,0.86,2)
vs=np.linspace(0.0,0.0,1)
vb=np.linspace(0.0,0.0,1)
P1.plotmodel(X1.draincurrent,1,2,vd,vg,vs,vb)
#############plot using plotfiledata
pathandfile = '../data/v7/finfet_7LgWt0.0076Lg1.0Hfin0.042tox0.0008WbR0.0038WbL0.0038vd0.86DATAREADY'
P1.updateparameter('symbol','-')
P1.plotfiledata(pathandfile,'gateOuterVoltage','drainTotalCurrent',1)
pathandfile = '../data/v7/finfet_7LgWt0.0076Lg1.0Hfin0.042tox0.0008WbR0.0038WbL0.0038vd0.05DATAREADY'
P1.updateparameter('color','b')
P1.plotfiledata(pathandfile,'gateOuterVoltage','drainTotalCurrent',1)
#####################################

P1.updateparameter('ylogflag',1)
###############plot using plotgeneral
P1.updateparameter('symbol','o')
vg=np.linspace(0,0.86,50)
vd=np.linspace(0.05,0.86,2)
vs=np.linspace(0.0,0.0,1)
vb=np.linspace(0.0,0.0,1)
P1.plotmodel(X1.draincurrent,2,2,vd,vg,vs,vb)
#############plot using plotfiledata
pathandfile = '../data/v7/finfet_7LgWt0.0076Lg1.0Hfin0.042tox0.0008WbR0.0038WbL0.0038vd0.86DATAREADY'
P1.updateparameter('symbol','-')
P1.plotfiledata(pathandfile,'gateOuterVoltage','drainTotalCurrent',2)
pathandfile = '../data/v7/finfet_7LgWt0.0076Lg1.0Hfin0.042tox0.0008WbR0.0038WbL0.0038vd0.05DATAREADY'
P1.updateparameter('color','b')
P1.plotfiledata(pathandfile,'gateOuterVoltage','drainTotalCurrent',2)
#####################################
print X1.ni

pylab.show() 
