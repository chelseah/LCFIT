#!/usr/bin/env python
import quadpolar
import quadsphere
import ZeipelModel
from Zeipel import Zeipel
#from Zeipel_all import Zeipel
from occultquad import occultquad
import numpy as np
import scipy as sp
from scipy.integrate import dblquad
import matplotlib 
from matplotlib import pyplot as plt
from const import *
import math
from math import cos,sin,tan,sqrt


def testgd(b,phi,theta):
    #phi and theta should be in radiant
    
    #parameters match Barnes 2009
    rp = 0.1
    a = 0.05
    P = 3.04
    beta = 0.19
    fratio = 0.1947
    Ms = 1.8
    Req = 2.029
    Rpole = Req*(1-fratio)
    Tpole = 8450
    u1 = 0.32; u2 = 0.32
    Protot = 8.64 
   
    #derived parameters and compute circular flux
    inc = np.arccos(Rpole*b*rsun/a/AU)
    tarr = -10000+np.arange(1000)/1000.*20000.
    t0 = 0.
    phase = (tarr-t0)/P/day
    sma = a*AU*np.sin(inc)/Req/rsun
    #print sma
    z=sma*np.sqrt(np.sin(phase*2*np.pi)*np.sin(phase*2*np.pi)+np.cos(phase*2*np.pi)*np.cos(phase*2*np.pi)*cos(inc)*cos(inc));
    circularflux = occultquad(z,u1,u2,rp/Req)[0]

    #derived parameters for gravity darkening bit 
    G_MRsun_hour = 4*math.pi**2*(1.5e3/7)**3./(365.*24.)**2
    Omega = 2*math.pi/Protot

    ggraveq0 = G_MRsun_hour*Ms/Req**2.
    groteq0 = Omega**2.*Req
    gpole = ggraveq0*Req**2./Rpole**2.
    #inite Zeipel 
    #Warning, use beta to init Zeipel, instead of y; F = g**(4*beta)
    gdmodel = Zeipel(fratio,phi,Req,ggraveq0,groteq0,beta)
    F = np.zeros(len(phase))
    Rpole = Req*(1-fratio)

    #call Zeiple
    #Warning, a is the semimajor axis in units of AU/rsun, instead of rstar;
    #input b*Rpole(in rsun,see above) instead of b
    gdmodel.Cal_F(phase,F,theta,a*AU/rsun,b*Rpole)
    lum = gdmodel.Cal_Lum()
    print lum/gpole**(4.*beta)/(4./3*np.pi*Req**3.)
    #correct the circular model
    model = (circularflux-max(circularflux))*F+1

    return model

def main():
    #barr = np.array([-0.9,-0.6,-0.3,0,0.3,0.6,0.9])
    barr = np.array([-0.6,-0.3,0,0.3,0.6])
    phi = 0.
    #phi = 30./180.*math.pi
    theta = 0.
    #theta = 30./180.*math.pi
    #barr = np.array([0,0.3,0.6,0.9])
    #barr = np.array([0.0])
    for b in barr:
        model = testgd(b,phi,theta)
        plt.plot(model)
    plt.ylim([0.995,1.0])
    plt.show()


if __name__ == '__main__':
    main()
