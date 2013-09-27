#!/usr/bin/env python
import math
import numpy as np
import scipy as sp
import os
from dataio import *
from HATlc import lightcurve as lc
from HATlc import tran
import random
import oblateness as Obl
from occultquad import occultquad

def binlc(filltime,time,mag):
    cadence = 30./60./24.
    nt = len(filltime)
    #nt = (int)((time[-1]-time[0])/cadence)
    #filltime = cadence*np.arange(nt)+time[0]
    fillmag = np.zeros(nt)
    for i in xrange(nt):
        samples=mag[np.abs(time-filltime[i])<cadence/2.]
        if(len(samples)>3):
            fillmag[i] = np.mean(samples)
    return fillmag

def gentran(time,period,epoch,q):
    ftime=sp.zeros(len(time))
    ftime=(time-epoch-0.5*period)/period-((time-epoch-0.5*period)/period).astype(int)
    ind=ftime<0
    ftime[ind]+=1
    #print min(ftime),max(ftime)
    intran=(ftime > (0.5-q/2.0))*(ftime < (0.5+q/2.0))
    #print ftime[intran])
    return intran


def trangen(time,mag,transit,lcflag=False):
    '''functions to generate fake lightcurves'''
    rmean = math.sqrt(transit.dip)
    f = 0.1
    inc = 89.209*math.pi/180.
    alpha = 0/180.*np.pi
    sma = 1./transit.q/np.pi
    #inc = math.acos(b/sma)
    u1 = 0.242
    u2 = 0.289
    req = rmean/math.sqrt(1-f)
    rpol = math.sqrt(1-f)*rmean
    fkmag=sp.zeros(len(time))+mag
    #fkmag=medianmag+np.random.randn(len(time))*stdmag
    ntran=int((max(time)-transit.epoch)/transit.P)-int((min(time)-transit.epoch)/transit.P)+1

    print ntran,transit.epoch,transit.q,max(time)
    tranwidth=transit.q*transit.P/2.
    obl = Obl.Oblateness(req,rpol,alpha,sma,inc,u1,u2)
    totalFlux = np.pi*(1.0-u1/3-u2/6)
    intran = gentran(time,transit.P,transit.epoch,transit.q)
    medianmag=np.median(mag[-intran])
    print rpol,medianmag
    stdmag=np.std(mag[-intran])
    phase = (time-transit.epoch)/transit.P-((time-transit.epoch)/transit.P).astype(int)
    if(lcflag):
        cadence = 1./60./24.
        nt = (int)((max(time)-min(time))/cadence)
        fktime = min(time)+np.arange(nt)*cadence
        fkintran = gentran(fktime,transit.P,transit.epoch,transit.q) 
        fkphase = (fktime-transit.epoch)/transit.P-((fktime-transit.epoch)/transit.P).astype(int)
        fkdflux = np.zeros(len(fktime[fkintran]))
        obl.relativeFlux(fkphase[fkintran],fkdflux)
        dflux = binlc(time[intran],fktime[intran],fkdflux)
    else:
        dflux = np.zeros(len(time[intran]))
        obl.relativeFlux(phase[intran],dflux)
    z=sma*np.sqrt(np.sin(phase*2*np.pi)*np.sin(phase*2*np.pi)+np.cos(phase*2*np.pi)*np.cos(phase*2*np.pi)*np.cos(inc)*np.cos(inc))
    circularflux = occultquad(z,u1,u2,rpol)[0]
    index = np.cos(phase*2*np.pi)<0
    circularflux[index]=1.0
    fkmag[intran] = medianmag+np.random.randn(len(fkmag[intran]))*stdmag
    fkmag[intran]+=1-(circularflux[intran]-dflux/totalFlux)
    #fkmag[intran]+=1-(circularflux[intran])
    #fkmag[intran]+=dflux/totalFlux

    return fkmag

	
def main():
    lcfile=lc()
    transit=tran()
    if(1):	
        path='./'
        os.chdir(path)
        #infile='HAT-365-0001481.epdlc'
        #infile='K972_short_notran.txt'
        #infile='data/kplr006603043-2011024051157_slc.tab'
        #infile='data/kplr006603043-2011145075126_slc.tab'
        infile='data/kplr006603043.ltf'
        coljd=1
        colmag=7
        #colmag=12
        transit.P=110.3216229;  
        rpstar=0.08453;
        transit.dip=rpstar**2.
        transit.q=transit.calq(rstar=2.5,logg=4.07)
        time=[];readcolumn(time,coljd,infile);time=np.array(time)
        #transit.epoch=min(time)+random.random();
        transit.epoch=1030.36382
        mag=[];readcolumn(mag,colmag,infile);mag=np.array(mag)
        fkmag=trangen(time,mag,transit,lcflag=True)
        outfile=os.path.splitext(infile)[0]+'.fktrap'
        fout=open(outfile,mode='w')
        for i in range(len(time)):
            fout.write('%13.6f %13.6f %13.6f\n' % (time[i],mag[i],fkmag[i]))
        fout.close()
    return


if __name__=='__main__':
	main()
