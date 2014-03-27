#!/usr/bin/env python
import quadpolar
import quadsphere
import LaraModel
from Lara import Lara
#from Zeipel_all import Zeipel
from occultquad import occultquad
import numpy as np
import scipy as sp
from scipy.integrate import dblquad
import matplotlib 
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import cm
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
    #Req = 1.0
    Rpole = Req*(1-fratio)
    Tpole = 8450
    u1 = 0.32; u2 = 0.32
    Protot = 8.2 
   
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
    w = (groteq0/ggraveq0)**0.5
    fratio = 1.-(1./(0.5*w**2.+1))

    #inite a gravity darkening model using Zeipel 
    #Warning, use beta to init Zeipel, instead of y; F = g**(4*beta)
    gdmodel = Lara(fratio,phi,Req,ggraveq0,groteq0)
    F = np.zeros(len(phase))
    Rpole = Req*(1-fratio)

    #use Zeiple to compute the flux correction for a transit
    #Warning, a is the semimajor axis in units of AU/rsun, instead of rstar;
    #input b*Rpole(in rsun,see above) instead of b
    gdmodel.Cal_F(phase,F,theta,a*AU/rsun,b*Rpole)
    
    F0 = np.array([0.])
    #use Zeiple to compute the total Flux, integrated along line of sight
    gdmodel.Cal_F0(F0)
    #the ratio of total Flux compare to use the Teff at pole 
    F0ratio = F0/(np.pi*Req**2.)
    print 'total Flux ratio along line of sight compare to piReq^2Tpole^4',F0ratio
    #use Zeiple to compute the luminosity, integrated over the 3-d luminosity
    lum = gdmodel.Cal_Lum()
    #the ratio of total lum compare to use the Teff at pole 
    lumratio = lum/(4.*np.pi*Req**2.)  
    print 'total luminosity compare to 4piReq^2Tpole^4',lumratio
    #print lum/gpole**(4.*beta)/(4./3*np.pi*Req**3.)
    #correct the circular model
    model = (circularflux-max(circularflux))*F+1

    return model

def test_theta():
    w = 0.7
    fr = 1.-(1./(0.5*w**2.+1))
    print fr
    Rp = 1.*(1-fr)
    e = math.sqrt(1.-(1-fr)**2.)
    S = 2*math.pi+math.pi*Rp**2./e*math.log((1.+e)/(1.-e))
    theta = np.linspace(0,np.pi/2.,40)
    inversR = np.sin(theta)**2./Rp**2.+np.cos(theta)**2.
    gdmodel = Lara(fr,0.,1.,1.,w**2.)
    Ts = np.zeros(len(theta))
    absR = 1.
    print np.sqrt(1./inversR)
    gdmodel.CalTtheta(theta,Ts,absR)
    lum = gdmodel.Cal_Lum()
    #print lum,S,Ts[0],Ts[-1]
    plt.plot(theta/np.pi*180.,Ts/(lum/4./math.pi)**0.25)
    #plt.plot(theta/np.pi*180.,Ts)
    
    w = 0.9
    fr = 1.-(1./(0.5*w**2.+1))
    print fr
    Rp = 1.*(1-fr)
    e = math.sqrt(1.-(1-fr)**2.)
    S = 2*math.pi+math.pi*Rp**2./e*math.log((1.+e)/(1.-e))
    theta = np.linspace(0,np.pi/2.,40)
    inversR = np.sin(theta)**2./Rp**2.+np.cos(theta)**2.
    gdmodel = Lara(fr,0.,1.,1.,w**2.)
    Ts = np.zeros(len(theta))
    absR = 1.
    print np.sqrt(1./inversR)
    gdmodel.CalTtheta(theta,Ts,absR)
    lum = gdmodel.Cal_Lum()
    #print lum,S,Ts[0],Ts[-1]
    plt.plot(theta/np.pi*180.,Ts/(lum/4./math.pi)**0.25)
    #
    plt.ylim([0.8,1.3]) 
    plt.show()
    return

def test_theta2():
    w = 0.7
    fr = 1.-(1./(0.5*w**2.+1))
    #print fr
    Req = 1.0
    Rpole = 1.*(1-fr)
    x = np.linspace(-Req,Req,1000)
    y = np.linspace(-Rpole,Rpole,1000)
    X, Y = np.meshgrid(x, y)

    gdmodel = Lara(fr,0.,1.,1.,w**2.)
    Txy = np.zeros([1000,1000])
    for i in xrange(1000):
        gdmodel.Cal_Txy(X[i,:],Y[i,:],Txy[i,:])

    e = math.sqrt(1.-(1-fr)**2.)
    S = 2*math.pi+math.pi*Rpole**2./e*math.log((1.+e)/(1.-e))
    #theta = np.linspace(0,np.pi/2.,40)
    #inversR = np.sin(theta)**2./Rp**2.+np.cos(theta)**2.
    #Ts = np.zeros(len(theta))
    #absR = 1.
    #print np.sqrt(1./inversR)
    #gdmodel.CalTtheta(theta,Ts,absR)
    #lum = gdmodel.Cal_Lum()
    #print lum,S
    #plt.plot(theta/np.pi*180.,Ts/(lum/S)**0.25)
    #plt.show()
    return




def test_fratio():
    rp = 0.1
    a = 0.05
    P = 3.04
    beta = 0.19
    fratio = 0.1947
    Ms = 1.8
    Req = 1.0
    #Req = 2.029
    Rpole = Req*(1-fratio)
    Tpole = 8450
    Protot = 8.2 
    #derived parameters for gravity darkening bit 
    G_MRsun = 4*math.pi**2*(1.5e3/7)**3./(365.*24.*3600)**2
    Omega = 2*math.pi/(Protot*3600.)

    ggraveq0 = G_MRsun*Ms/Req**2.
    groteq0 = Omega**2.*Req
    gpole = ggraveq0*Req**2./Rpole**2.

    x = np.linspace(-Req,Req,1000)
    y = np.linspace(-Rpole,Rpole,1000)
    X, Y = np.meshgrid(x, y)
    fratios = np.linspace(0,0.33,20)
    #fratios=[0.]
    #phis = [0.0]
    #fig = plt.figure(figsize=[8,8])
    count = 1
    print groteq0/ggraveq0
    phi = 0.
    Tminarr=np.zeros(20)
    Tminthe=np.zeros(20)
    Tmaxarr=np.zeros(20)
    j=0
    for fr in fratios:
        #phi-=np.pi/2.
        w = (2.*fr/(1.-fr))**0.5
        gdmodel = Lara(fr,phi,Req,ggraveq0,ggraveq0*w**2.)

        Txy = np.zeros([1000,1000])
        for i in xrange(1000):
            gdmodel.Cal_Txy(X[i,:],Y[i,:],Txy[i,:])
        exterior = Txy==0 
        T_mask = np.ma.masked_array(Txy,mask=exterior)
        Tminarr[j] = np.min(T_mask)
        Tmaxarr[j] = np.max(T_mask)
        Tminthe[j] = np.sqrt(2./(2.+w*w))*(1-w**2.)**(1./12.)*np.exp(-4./3.*w*w/(2+w*w)**3.)
        print fr,w,Tminarr[j],Tmaxarr[j]
        j+=1
    plt.plot(fratios,Tminarr)
    plt.plot(fratios,Tminthe)
    plt.show()
    return
def test_Teff():
    #phi = -2.*np.pi/3.
    rp = 0.1
    a = 0.05
    P = 3.04
    beta = 0.19
    fratio = 0.1947
    Ms = 1.8
    #Req = 2.029
    Req = 2.0
    Rpole = Req*(1-fratio)
    Tpole = 8450
    Protot = 8.2 
    #derived parameters for gravity darkening bit 
    G_MRsun = 4*math.pi**2*(1.5e3/7)**3./(365.*24.*3600)**2
    Omega = 2*math.pi/(Protot*3600.)

    ggraveq0 = G_MRsun*Ms/Req**2.
    groteq0 = Omega**2.*Req
    gpole = ggraveq0*Req**2./Rpole**2.

    x = np.linspace(-Req,Req,100)
    y = np.linspace(-Rpole,Rpole,100)
    X, Y = np.meshgrid(x, y)
    phis = [0.0,30.0,60.0,90.0,120.0,150.0,180.0]
    #phis = [0.0]
    fig = plt.figure(figsize=[35,3.5])
    #fig = plt.figure(figsize=[8,8])
    count = 1
    print groteq0/ggraveq0
    for phi in phis:
        #w = (1./(1.-fratio)**2.-1.)**0.5
        phi/=(180./np.pi)
        #gdmodel = Lara(fratio,phi,1.0,ggraveq0,ggraveq0*w)
        #phi-=np.pi/2.
        w = math.sqrt(groteq0/ggraveq0)
        fratio = 1.-(1./(0.5*w**2.+1))
        print fratio,w
        gdmodel = Lara(fratio,phi,Req,ggraveq0,groteq0)

        Txy = np.zeros([100,100])
        for i in xrange(100):
            gdmodel.Cal_Txy(X[i,:],Y[i,:],Txy[i,:])
        #print Txy
        subindex = 170+count
        #subindex = 110+count
        count+=1
        ax = fig.add_subplot(subindex)
        exterior = np.sqrt(((X/Req)**2) + ((Y/Rpole)**2)) > 1.0
        T_mask = np.ma.masked_array(Txy,mask=exterior)
        print np.max(T_mask),np.min(T_mask)
        ellipse = Ellipse([0,0],2*Req,2*Rpole,fc='none',ec='k') 
        #fig = plt.figure()
        #colors = (Txy-min(Txy))*1.0/(max(Txy)-min(Txy))
        #arr = [cmap[i-1] for i in index]
        #CS = ax.scatter(X,Y,c=colors,cmap=cm.copper)
        #CS.set_clim([0.1,1.0])
        CS = ax.contourf(X,Y,Tpole*T_mask,cmap=cm.copper)
        #CS.set_clim([7600,8500])
        cbar = plt.colorbar(CS)
        ax.set_xlim([-3,3])
        ax.set_ylim([-3,3])
        ax.add_patch(ellipse)
        ax.set_title(phi)
    plt.show()
    #for i in xrange(len(X)):
    #    print X[i],Y[i],Txy[i],colors[i]
def main():
    #barr = np.array([-0.9,-0.6,-0.3,0,0.3,0.6,0.9])
    barr = np.array([-0.6,-0.3,0,0.3,0.6])
    phi = 0.
    #phi = 90./180.*math.pi
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
    #test_Teff()
    #test_fratio()
    #test_theta()
    #test_theta2()
