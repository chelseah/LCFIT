#!/usr/bin/env python
import numpy as np
import scipy as sp
from math import sqrt,cos
import oblateness as Obl
import oblatenessfast as OblF
from occultquad import occultquad
#from Eastman_occultquad import occultquad
import matplotlib
from matplotlib import pyplot as plt
import time

def Scaling():
    #rmean = 0.1
    rmeanarr = (1+np.arange(20))/20.*0.2
    period = 100.
    u1 = 0
    u2 = 0
    #f = 0.1
    #farr = (1+np.arange(5))/5.*0.2
    farr = (1+np.arange(10))/10.*0.2
    #farr = np.array([0.1])
    #alphaarr =(np.arange(50)/50.*180.-90)
    alphaarr = np.array([45])
    sma = (period/365)**(2./3.)*1.5e13/7.e10
    #inc = 89.709/180.*np.pi
    b0arr = -0.8+np.arange(101)/100.*1.6
    #b0arr = np.array([0.5])
    #b0 = 0.8
    #b0 = sma*cos(inc)
    #print sma,b0
    Npoint = 1000
    percent = 0.025
    dflux = np.zeros(Npoint)
    dfluxF = np.zeros(Npoint)
    for l in xrange(len(rmeanarr)):
        for k in xrange(len(b0arr)):
            for i in xrange(len(alphaarr)):
                for j in xrange(len(farr)):
                    rmean = rmeanarr[l]
                    f = farr[j]
                    alpha=alphaarr[i]/180.*np.pi
                    inc = np.arccos(b0arr[k]/sma)
                    b0 = b0arr[k]
                    req = rmean/sqrt(1-f)
                    rpol = sqrt(1-f)*rmean
                    oblf = OblF.Oblateness(req,rpol,alpha,sma,inc,u1,u2)
                    if(b0>(1+req)):
                        print 'no transit happen'
                        return
                    dphi = 2*percent/(Npoint-1)
                    totalFlux = np.pi*(1.0-u1/3-u2/6)
                    phi=-1*percent+dphi*np.arange(Npoint)
                    #call
                    oblf.relativeFlux(phi,dfluxF)
                    #print phi,dflux
                    z=sma*np.sqrt(np.sin(phi*2*np.pi)*np.sin(phi*2*np.pi)+np.cos(phi*2*np.pi)*np.cos(phi*2*np.pi)*cos(inc)*cos(inc));
                    #z = abs(np.sin(phi*np.pi))*sma#/abs(np.sin(inc))
                    #print '3',time.time(),time.clock()-start
                    start = time.clock()
                    circularflux = occultquad(z,u1,u2,rpol)[0]
                    index = np.cos(phi*2*np.pi)<0
                    circularflux[index]=1.0
                    circularfluxmean = occultquad(z,u1,u2,rmean)[0]
                    circularfluxmean[index]=1.0
                    residual = (-circularfluxmean+circularflux-dfluxF/totalFlux)/1.e-6
                    #plt.xlim([-0.01,0.01])
                    ##plt.plot(phi,circularfluxmean)
                    #plt.plot(phi,residual)
                    #plt.show()
                    print rmeanarr[l],alphaarr[i],farr[j],b0arr[k],max(residual) 
    return



def main():
    #rmean = 0.154679
    rmean = 0.08453
    f = 0.1
    f2 = 0.1
    #alpha =45./180.*np.pi
    alpha =-70./180.*np.pi
    alpha2 =70./180.*np.pi
    #sma = 8.924
    sma = 49.584
    #period = 2.218573 
    period = 110.3216229
    #inc = 85.749/180.*np.pi
    #inc = 89.209/180.*np.pi
    inc = 89.51/180.*np.pi
    #inc = 90./180.*np.pi
    u1 = 0.242
    u2 = 0.289
    Npoint = 1000
    percent = 0.025
    #percent = 1.0
    dflux = np.zeros(Npoint)
    dfluxF = np.zeros(Npoint)
    req = rmean/sqrt(1-f)
    rpol = sqrt(1-f)*rmean
    req2 = rmean/sqrt(1-f2)
    rpol2 = sqrt(1-f2)*rmean
    #initial
    start = time.clock()
    #print '0',time.time(),time.clock()
    obl = Obl.Oblateness(req,rpol,alpha,sma,inc,u1,u2)
    oblf = OblF.Oblateness(req2,rpol2,alpha2,sma,inc,u1,u2)
    print '0',time.time(),time.clock()-start
    b0 = sma*cos(inc)
    print b0
    if(b0>(1+req)):
        print 'no transit happen'
        return
    dphi = 2*percent/(Npoint-1)
    totalFlux = np.pi*(1.0-u1/3-u2/6)
    phi=-1*percent+dphi*np.arange(Npoint)
    #call
    start = time.clock()
    #print '1',time.time(),time.clock()-start
    obl.relativeFlux(phi,dflux)
    print '1old',time.time(),time.clock()-start
    start = time.clock()
    oblf.relativeFlux(phi,dfluxF)
    print '1',time.time(),time.clock()-start
    #print phi,dflux
    z=sma*np.sqrt(np.sin(phi*2*np.pi)*np.sin(phi*2*np.pi)+np.cos(phi*2*np.pi)*np.cos(phi*2*np.pi)*cos(inc)*cos(inc));
    #z = abs(np.sin(phi*np.pi))*sma#/abs(np.sin(inc))
    #print '3',time.time(),time.clock()-start
    start = time.clock()
    circularflux = occultquad(z,u1,u2,rpol)[0]
    circularflux2 = occultquad(z,u1,u2,rpol2)[0]
    print '3',time.time(),time.clock()-start
    index = np.cos(phi*2*np.pi)<0
    circularflux[index]=1.0
    #print z.shape,circularflux.shape
    #plt.plot(phi,circularflux)
    #plt.plot(phi,z)
    #plt.plot(phi,dflux/totalFlux)
    #plt.xlim([-0.006,0.006])
    start = time.clock()
    #print '5',time.time(),time.clock()-start
    circularfluxmean = occultquad(z,u1,u2,rmean)[0]
    print '4',time.time(),time.clock()-start
    #plt.plot(phi,circularflux-dflux/totalFlux)
    plt.xlim([-0.01,0.01])
    plt.plot(phi,(-circularfluxmean+circularflux-dflux/totalFlux)/1.e-6)
    plt.plot(phi,(-circularfluxmean+circularflux2-dfluxF/totalFlux)/1.e-6,'+')
    plt.show()
    return

if __name__=='__main__':
    main()
    #Scaling()
