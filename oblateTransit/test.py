#!/usr/bin/env python
import numpy as np
import scipy as sp
from math import sqrt,cos
import oblateness as Obl
from occultquad import occultquad
import matplotlib
from matplotlib import pyplot as plt
def main():
    rmean = 0.1
    f = 0.1
    alpha = 45/180.*np.pi
    sma = 8.924
    period = 2.218573 
    inc = 89.7/180.*np.pi
    u1 = 0.076
    u2 = 0.034
    Npoint = 500
    percent = 0.025
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
    z = abs(np.sin(phi*np.pi*2.))*sma
    circularflux = occultquad(z,u1,u2,rpol)[0]
    #print z.shape,circularflux.shape
    #plt.plot(phi,circularflux)
    #plt.plot(phi,circularflux-dflux/totalFlux)
    circularfluxmean = occultquad(z,u1,u2,rmean)[0]
    plt.plot(phi,circularfluxmean-circularflux+dflux/totalFlux)
    plt.show()
    return

if __name__=='__main__':
    main()
