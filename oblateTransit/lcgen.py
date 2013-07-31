#!/usr/bin/env python
import math
import numpy as np
import scipy as sp
import os
from dataio import *
from util import *
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


def trangen(time,mag,transit):
    '''functions to generate fake lightcurves'''
    rmean = math.sqrt(transit.dip)
    f = 0.098
    alpha = 0/180.*np.pi
    sma = 1./transit.q/np.pi
    inc = math.arccos(b/sma)
    u1 = 0.076
    u2 = 0.034
    req = rmean/sqrt(1-f)
    rpol = sqrt(1-f)*rmean

    medianmag=np.median(mag)
    stdmag=np.std(mag)
    fkmag=sp.zeros(len(time))+mag
    #fkmag=medianmag+np.random.randn(len(time))*stdmag
    ntran=int((max(time)-transit.epoch)/transit.P)-int((min(time)-transit.epoch)/transit.P)+1

    print ntran,transit.epoch,transit.q,max(time)
    tranwidth=transit.q*transit.P/2.
    obl = Obl.Oblateness(req,rpol,alpha,sma,inc,u1,u2)
    totalFlux = np.pi*(1.0-u1/3-u2/6)
    for i in range(ntran):
        midtran=transit.epoch+i*transit.P+int((min(time)-transit.epoch)/transit.P)*transit.P 
        phase = (time-midtran)/transit.P
        intran = gentran(time,transit.P,transit.epoch,transit.q)
        dflux = np.zeros(len(time[intran]))
        obl.relativeFlux(phi,dflux)
        z=sma*np.sqrt(np.sin(phi*2*np.pi)*np.sin(phi*2*np.pi)+np.cos(phi*2*np.pi)*np.cos(phi*2*np.pi)*cos(inc)*cos(inc))
        circularflux = occultquad(z,u1,u2,rpol)[0]
        index = np.cos(phi*2*np.pi)<0
        circularflux[index]=1.0
        fkmag[intran]+=1-(circularflux-dflux/totalflux)

    return fkmag

	
def main():
	lcfile=lc()
	transit=tran()
	if(1):	
		path='./'
		os.chdir(path)
		#infile='HAT-365-0001481.epdlc'
		#infile='K972_short_notran.txt'
		infile='K972_llc_notran.txt'
		coljd=1
		#colmag=3
		#colmag=12
		colmag=2
		transit.P=13.118908; rpstar=0.019;transit.dip=rpstar**2.
  		transit.q=transit.calq(rstar=2.92,logg=3.82)
		time=[];readcolumn(time,coljd,infile);time=np.array(time)
		#transit.epoch=min(time)+random.random();
		transit.epoch=1617.2
		mag=[];readcolumn(mag,colmag,infile);mag=np.array(mag)
		fkmag=trangen(time,mag,transit)
		outfile=os.path.splitext(infile)[0]+'.fktrap'
		fout=open(outfile,mode='w')
		for i in range(len(time)):
			fout.write('%13.6f %13.6f %13.6f\n' % (time[i],mag[i],fkmag[i]))
		fout.close()
	return


if __name__=='__main__':
	main()
