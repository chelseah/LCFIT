#!/usr/bin/env python
# this file is to solve long trend filting by least square method using python's own routine.
import math
import numpy as np
import scipy as sp
import copy
from scipy import linalg
from scipy.optimize import leastsq
from scipy import interpolate
from scipy import signal
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt

def matrixgen(time,n,T):
	'''generate least square fitting matrix with n cos filter,t is the total time span, formulism refer to Huang and Bacos (2012) eq [1]'''
	tn=len(time)
	a=sp.zeros((tn,n))
	a=np.rollaxis(np.sin(np.array(np.r_['c',0:n]*time[np.arange(tn)])*math.pi/T),1,0)
	
	return a

def gentran(time,period,epoch,q):
    ftime=sp.zeros(len(time))
    ftime=(time-epoch-0.5*period)/period-((time-epoch-0.5*period)/period).astype(int)
    ind=ftime<0
    ftime[ind]+=1
    #print min(ftime),max(ftime)
    intran=(ftime > (0.5-q/2.0))*(ftime < (0.5+q/2.0))
    #print ftime[intran])
    return intran


def lssubrecon(otime,oflux,intran,n=30,noplot=True):
    length=len(oflux)
    ctime = otime[-intran]
    cflux = oflux[-intran]-np.mean(oflux[-intran]) 
    k = (cflux[-1]-cflux[0])/(ctime[-1]-ctime[0])
    b = cflux[0]-k*ctime[0]
    cflux -=(k*ctime+b)
    E0=min(ctime)	
    T=max(ctime)-min(ctime)
    #print E0,len(cflux),len(ctime)
    A=matrixgen(ctime-E0,n,T)
    #print A.shape,n
    c,resid,rank,sigma = linalg.lstsq(A,cflux)
    #print resid	
    rflux=sp.zeros(len(otime))	
    e=np.rollaxis(np.sin(np.array(np.r_['c',0:n]*(otime[np.arange(length)]-E0))*math.pi/T),1,0)
    rflux=np.dot(e,c)
    eflux = k*otime+b
    #print rflux.shape,eflux.shape 
    rflux+=eflux
    if(not noplot):
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(otime,rflux,'.',otime,eflux,'+',ctime,cflux,'x')
        plt.show()
    	
    dflux=oflux-rflux+np.mean(oflux)
    
    return dflux

	
def flatfunc(c,ctime,T):
    n = len(c)-1
    b=np.rollaxis(np.sin(np.array(np.r_['c',0:n]*(ctime))*math.pi/T),1,0)
    rflux=np.dot(b,c[0:n])
    rflux-=np.mean(rflux)
    return rflux

	
def tranfunc(c,ctime,intran,T):
    n = len(c)-1
    b=np.rollaxis(np.sin(np.array(np.r_['c',0:n]*(ctime))*math.pi/T),1,0)
    rflux=np.dot(b,c[0:n])
    try:
        rflux[intran]+=c[n]
    except TypeError:
        print c
        raise
    rflux-=np.mean(rflux)
    return rflux

def lsfitrecon(otime,oflux,intran,n=30,noplot=True,dipguess=0.0001):

    length=len(oflux)
    
    ctime=np.array(otime)
    E0=min(ctime)	
    T=max(ctime)-min(ctime)
    cflux=np.array(oflux)-np.mean(oflux)
    e = lambda c,time,index,flux,T:(tranfunc(c,time,index,T)-flux)
    c0=list(np.zeros(n))
    c0.append(dipguess)
    c,success = leastsq(e,c0,args=(ctime-E0,intran,cflux,T),maxfev=10000)
    rflux=sp.zeros(len(ctime))	
    b=np.rollaxis(np.sin(np.array(np.r_['c',0:n]*(ctime[np.arange(length)]-E0))*math.pi/T),1,0)
    rflux=np.dot(b,c[0:n])
    tran = sp.zeros(len(ctime))
    tran[intran]+=c[n]
    tran-=np.mean(tran)
    #print c[n],len(tran[intran])
    if(not noplot):
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(ctime,cflux,'x',ctime,rflux,'.',ctime,tran,'+')
        plt.show()
        plt.close(fig)
    	
    rL=len(rflux)
    nn=len(otime)
    dflux=oflux-rflux+np.mean(oflux)
    
    return dflux

