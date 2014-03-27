#!/usr/bin/env python
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
from dataio import readcolumn 
from HATlc import tran
from const import *
import math
from occultquad import occultquad
from lcfromdata import fluxtomag
import time
def main():
    planets = ["Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"]
    infile = "solarsystem.txt"
    period = []; readcolumn(period,1,infile)
    inc = []; readcolumn(inc,2,infile)
    size = []; readcolumn(size,3,infile)
    a = []; readcolumn(a,4,infile)
    phase = -0.05+0.1*np.arange(10000)/10000.
    offsetarr = np.linspace(-0.5,8,400)
    #offsetarr = np.array([0.])
    #plt.ion()
    #plt.show()
    time.sleep(1.0)
    count=0
    for offset in offsetarr:
        fig = plt.figure(figsize=[8,8])
        plt.clf()
        fig.suptitle("offset from the eliptical = %.2f" % offset, color='red')
        for i in xrange(len(planets)):
            transit = tran()
            transit.P = period[i]*365.
            transit.inc = (90.-offset+inc[i])/180.*math.pi
            transit.u1 = 0.4352 
            transit.u2= 0.2965
            transit.sma = (a[i] / rsun*AU)
            z = transit.calZ(phase)
            #print z
            #break
            rmean = size[i]*6.378e8/rsun
            circularflux = occultquad(z,transit.u1,transit.u2,rmean)[0]
            #print transit.inc,rmean
            yformatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        
            ax = fig.add_subplot(4,2,i+1)
            ax.text(-48,1+rmean**2.*0.1,planets[i])
            #if i%2==0:
            #    ax.set_title(planets[i],loc='left')
            #else:
            #    ax.set_title(planets[i],loc='right')

            ax.set_xlim([-60,60])
            ax.set_ylim([1-rmean**2.*1.4,1+rmean**2.*0.3])
            ax.plot(phase*transit.P*24.,circularflux)
            ax.yaxis.set_major_formatter(yformatter)
            #plt.plot(z,circularflux)
        #plt.tight_layout()
        plt.savefig('img-%.3d.png' % count)
        #plt.savefig('solar.png')
        count+=1
        plt.close()
        #break
        #plt.tight_layout()
        #plt.draw()
        #time.sleep(0.001)
    return

def cal_mag(dis):
    d0 = AU
    m0 = -26.74
    dis *= 3.0857e18
    mag = m0-2.5*np.log10(d0**2./dis**2.)
    return mag

def cal_sigma(mag):
    c=3.46*10**(0.4*(12-mag)+8)
    sigma = np.sqrt(c+7.e6*max(1,mag/14.)**4.)/c
    return sigma

def distance():
    planets = ["Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"]
    infile = "solarsystem.txt"
    period = []; readcolumn(period,1,infile)
    inc = []; readcolumn(inc,2,infile)
    size = []; readcolumn(size,3,infile)
    a = []; readcolumn(a,4,infile)
    phase = -0.01+0.02*np.arange(10000)/10000.
    disarr = np.logspace(np.log10(1.),np.log10(1.e4),100)
    #plt.ion()
    #plt.show()
    time.sleep(1.0)
    count=0
    for dis in disarr:
        magsun = cal_mag(dis)
        sigma = cal_sigma(magsun)
        fig = plt.figure(figsize=[8,8])
        #plt.clf()
        #fig.suptitle("offset from the eliptical = %.2f" % offset, color='red')
        fig.suptitle("distance from sun = %.2f pc, Kepmag=%.2f, $\sigma$ = %.2f ppm" % (dis, magsun, sigma/1.e-6),color='red')
        for i in xrange(len(planets)):
            transit = tran()
            transit.P = period[i]*365.
            transit.inc = (90.)/180.*math.pi
            transit.u1 = 0.4352 
            transit.u2= 0.2965
            transit.sma = (a[i] / rsun*AU)
            z = transit.calZ(phase)
            #print z
            #break
            rmean = size[i]*6.378e8/rsun
            circularflux = occultquad(z,transit.u1,transit.u2,rmean)[0]
            circularmag = np.zeros(len(phase))+np.random.normal(loc=0,scale=sigma,size=phase.shape)+fluxtomag(circularflux)
            #circularmag = fluxtomag(circularflux)
            #print transit.inc,rmean
            yformatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        
            ax = fig.add_subplot(4,2,i+1)
            ax.text(-48,1-rmean**2.*0.3,planets[i])
            #if i%2==0:
            #    ax.set_title(planets[i],loc='left')
            #else:
            #    ax.set_title(planets[i],loc='right')

            ax.set_xlim([-60,60])
            ax.set_ylim([1+rmean**2.*1.4,1-rmean**2.*0.5])
            ax.plot(phase*transit.P*24.,1+circularmag)
            ax.yaxis.set_major_formatter(yformatter)
            #plt.plot(z,circularflux)
        #plt.tight_layout()
        plt.savefig('imgmag-%.3d.png' % count)
        plt.savefig('solar.png')
        #plt.show()
        count+=1
        plt.close()
        #break

    return

if __name__=='__main__':
    main()
    #distance()
