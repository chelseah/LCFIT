#!/usr/bin/env python
import numpy as np
import scipy as sp
from math import sqrt,cos
import oblateness as Obl
from occultquad import occultquad
#from Eastman_occultquad import occultquad
import matplotlib
from matplotlib import pyplot as plt
def main():
    #rmean = 0.154679
    rmean = 0.08453
    f = 0.1
    alpha = 45./180.*np.pi
    #sma = 8.924
    sma = 49.584
    #period = 2.218573 
    period = 110.3216229
    #inc = 85.749/180.*np.pi
    inc = 89.209/180.*np.pi
    u1 = 0.242
    u2 = 0.289
    Npoint = 500
    percent = 0.025
    #percent = 1.0
    dflux = np.zeros(Npoint)
    req = rmean/sqrt(1-f)
    rpol = sqrt(1-f)*rmean
    #initial
    obl = Obl.Oblateness(req,rpol,alpha,sma,inc,u1,u2)
    b0 = sma*cos(inc)
    if(b0>(1+req)):
        print 'no transit happen'
        return
    dphi = 2*percent/(Npoint-1)
    totalFlux = np.pi*(1.0-u1/3-u2/6)
    phi=-1*percent+dphi*np.arange(Npoint)
    #call
    obl.relativeFlux(phi,dflux)
    #print phi,dflux
    z=sma*np.sqrt(np.sin(phi*2*np.pi)*np.sin(phi*2*np.pi)+np.cos(phi*2*np.pi)*np.cos(phi*2*np.pi)*cos(inc)*cos(inc));
    #z = abs(np.sin(phi*np.pi))*sma#/abs(np.sin(inc))
    circularflux = occultquad(z,u1,u2,rpol)[0]
    index = np.cos(phi*2*np.pi)<0
    circularflux[index]=1.0
    #print z.shape,circularflux.shape
    #plt.plot(phi,circularflux)
    #plt.plot(phi,z)
    #plt.plot(phi,dflux/totalFlux)
    #plt.xlim([-0.006,0.006])
    circularfluxmean = occultquad(z,u1,u2,rmean)[0]
    #plt.plot(phi,circularflux-dflux/totalFlux)
    #plt.xlim([-0.01,0.01])
    plt.plot(phi,(-circularfluxmean+circularflux-dflux/totalFlux)/1.e-6)
    plt.show()
    return

if __name__=='__main__':
    main()
