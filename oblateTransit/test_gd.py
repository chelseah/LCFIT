#!/usr/bin/env python
import Zeipel
from occultquad import occultquad
import numpy as np
import scipy as sp
from scipy.integrate import dblquad
import matplotlib 
from matplotlib import pyplot as plt
from const import *
import math
from math import cos,sin,tan,sqrt

def int_area():
    #def f(x,y):
    #    return 1.
    b = 0.0
    rp = 0.1
    a = 0.05
    P = 3.04
    beta = 0.19
    fratio = 0.1947
    Ms = 1.8
    Req = 2.029
    Tpole = 8450
    theta = 0.0
    phi = 0.0
    inc = np.arccos(Req*b*rsun/a/AU)
    u1 = 0.32; u2 = 0.32
    Protot = 8.64 

    G_MRsun_hour = 4*math.pi**2*(1.5e3/7)**3./(365.*24.)**2
    Omega = 2*math.pi/Protot
    ggraveq0 = G_MRsun_hour*Ms/Req**2.
    groteq0 = Omega**2.*Req
    ggraveq = 1.0
    groteq = ggraveq/ggraveq0*groteq0
    Rpole = Req*(1-fratio)
    gpole=ggraveq*Req**2/Rpole**2.
    beta = 1.0

    def f(x,y):
        #x,y is r, theta
        x0 = x*cos(y)
        y0 = x*sin(y)
        #print x0,y0
        g = Zeipel.cal_geff(x0,y0,fratio,phi,1.,ggraveq,groteq)
        return g**beta
    area = dblquad(lambda x,y: x*f(x,y),0,2.*np.pi,lambda x: 0,lambda x: Req)
    print area[0]/np.pi
    return
def testgd(b):
    rp = 0.1
    a = 0.05
    P = 3.04
    beta = 0.19
    fratio = 0.1947
    Ms = 1.8
    Req = 2.029
    Rpole = Req*(1-fratio)
    Tpole = 8450
    #theta = 90.0
    theta = 0.0
    #phi = 0.0
    phi = 30./180.*math.pi
    #phi = 90./180.*math.pi
    #inc = 90.0/180.*math.pi
    #inc = np.arccos(Req*b*rsun/a/AU)
    inc = np.arccos(Rpole*b*rsun/a/AU)
    #print inc/math.pi*180.
    tarr = -10000+np.arange(500)/500.*20000.
    t0 = 0.
    phase = (tarr-t0)/P/day
    u1 = 0.32; u2 = 0.32
    Protot = 8.64 
    sma = a*AU*np.sin(inc)/Req/rsun
    #print sma
    z=sma*np.sqrt(np.sin(phase*2*np.pi)*np.sin(phase*2*np.pi)+np.cos(phase*2*np.pi)*np.cos(phase*2*np.pi)*cos(inc)*cos(inc));
    circularflux = occultquad(z,u1,u2,rp/Req)[0]
    #plt.plot(phi,circularflux)
    #plt.show()

    theta = theta*math.pi/180.
    
    #Protot = Protot #Protot is in hour
    G_MRsun_hour = 4*math.pi**2*(1.5e3/7)**3./(365.*24.)**2
    Omega = 2*math.pi/Protot

    ggraveq0 = G_MRsun_hour*Ms/Req**2.
    groteq0 = Omega**2.*Req
    #ggraveq = 1.0
    #groteq = ggraveq/ggraveq0*groteq0
    Rpole = Req*(1-fratio)
    #gpole=ggraveq*Req**2/Rpole**2.
    #x = tarr*math.pi*a*AU*cos(theta)/P/day+b*sin(theta)
    #y = -1*x*tan(theta)+rsun*Req*b*cos(theta)+rsun*Req*b*sin(theta)*tan(theta)
    x = phase*2.*math.pi*a*AU*cos(theta)-Rpole*rsun*b*sin(theta)
    #y = -1*x*tan(theta)+rsun*Rpole*b*cos(theta)+rsun*Rpole*b*sin(theta)*tan(theta)
    y = -phase*2.*math.pi+a*AU*sin(theta)+rsun*Rpole*b*cos(theta)
    #x = x/(rsun*Req)
    #y = y/(rsun*Req)
    x = x/(rsun)
    y = y/(rsun)
    #print x,y
    #print phi,fratio,ggraveq,groteq
    g = Zeipel.cal_geff(x,y,fratio,phi,Req,ggraveq0,groteq0)
    #g = Zeipel.cal_geff(x,y,fratio,phi,1.,1,groteq/ggraveq)
    #g0 = Zeipel.cal_geff(0,0,fratio,phi,1.,ggraveq,groteq)
   
    def func(x,y):
        #x,y is r, theta
        x0 = x*cos(y)
        y0 = x*sin(y)
        #print x0,y0
        g = Zeipel.cal_geff(x0,y0,fratio,phi,Req,ggraveq0,groteq0)
        #g = Zeipel.cal_geff(x0,y0,fratio,phi,1.0,1.0,groteq/ggraveq)
        return g**(4.*beta)
    #F0 = dblquad(lambda x,y: x*f(x,y),0,2.*np.pi,lambda x: 0,lambda x: 1.0)
    F0 = dblquad(lambda x,y: x*func(x,y),0,2.*np.pi,lambda x: 0,lambda x: Req)
    feff = 1- ((1-fratio)**2.*np.cos(phi)**2.+np.sin(phi)**2.)**0.5

    F = (g)**(4.*beta)/(F0[0]/np.pi/Req**2.*(1-feff))
    #F = (g/gpole)**(4.*beta)
    
    #model = (model-max(model))*(1-fratio_i)*F+1
    
    #model = (circularflux-max(circularflux))*(1-fratio)*F+1
    model = (circularflux-max(circularflux))*F+1
    #print F
    #model = circularflux
    #model = (circularflux-max(circularflux))*F+1
    #plt.plot(phi,circularflux)
    #plt.show()

    return model

def main():
    barr = np.array([-0.6,-0.3,0,0.3,0.6])
    #barr = np.array([0,0.3,0.6])
    #barr = np.array([0.0])
    for b in barr:
        model = testgd(b)
        plt.plot(model)
    plt.ylim([0.995,1.0])
    plt.show()

def testint():
    int_area()

if __name__ == '__main__':
    main()
    #testint()
