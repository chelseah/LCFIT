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
import matplotlib
from matplotlib import pyplot as plt
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
    f = transit.f
    inc = transit.inc*math.pi/180.
    alpha = transit.alpha*np.pi/180.
    u1 = transit.u1
    u2 = transit.u2
    #u1 = 0.076; u2 = 0.034; 
    sma = 1./transit.q/np.pi
    print sma
    #inc = math.acos(b/sma)
    #rmean = 0.1
    #f = 0.098

    req = rmean/math.sqrt(1-f)
    rpol = math.sqrt(1-f)*rmean
    #fkmag=sp.zeros(len(time))+mag
    ntran=int((max(time)-transit.epoch)/transit.P)-int((min(time)-transit.epoch)/transit.P)+1
    
    #print ntran,transit.epoch,transit.q,max(time),u1,u2,sma,alpha*180./np.pi,inc*180./np.pi
    #print rmean,req,rpol
       
    tranwidth=transit.q*transit.P/2.
    obl = Obl.Oblateness(req,rpol,alpha,sma,inc,u1,u2)
    totalFlux = np.pi*(1.0-u1/3-u2/6)
    intran = gentran(time,transit.P,transit.epoch,transit.q)
    medianmag=np.median(mag[-intran])
    stdmag=np.std(mag[-intran])
    stdmag = 1.e-4
    #print rpol,medianmag
    phase = (time-transit.epoch)/transit.P-((time-transit.epoch)/transit.P).astype(int)
    ind = phase>0.5
    phase[ind]-=1.0

    fkmag=medianmag+np.random.randn(len(time))*stdmag
    if(lcflag):
        circularfluxmeanlc=medianmag+np.random.randn(len(time))*stdmag
        cadence = 1./60./24.
        nt = (int)((max(time)-min(time))/cadence)
        fktime = min(time)+np.arange(nt)*cadence
        fkmagsc=medianmag+np.random.randn(len(fktime))*stdmag
        fkintran = gentran(fktime,transit.P,transit.epoch,transit.q) 
        fkphase = (fktime-transit.epoch)/transit.P-((fktime-transit.epoch)/transit.P).astype(int)
        ind = fkphase>0.5
        fkphase[ind]-=1.0
        fkdflux = np.zeros(len(fktime[fkintran]))
        obl.relativeFlux(fkphase[fkintran],fkdflux)
        z=sma*np.sqrt(np.sin(fkphase*2*np.pi)*np.sin(fkphase*2*np.pi)+np.cos(fkphase*2*np.pi)*np.cos(fkphase*2*np.pi)*np.cos(inc)*np.cos(inc))
        circularfluxmean = occultquad(z,u1,u2,rmean)[0]
        circularflux = occultquad(z,u1,u2,rpol)[0]
        index = np.cos(fkphase*2*np.pi)<0
        circularflux[index]=1.0
        circularfluxmean[index]=1.0
        fkmagsc[fkintran] = medianmag
        fkmagsc[fkintran]+=1-(circularflux[fkintran]-fkdflux/totalFlux)
        circularfluxmean=1-circularfluxmean+medianmag
        fkmag[intran] = binlc(time[intran],fktime[fkintran],fkmagsc[fkintran])
        residual = (fkmagsc-circularfluxmean)/1.e-6
        circularfluxmeanlc[intran] = binlc(time[intran],fktime[fkintran],circularfluxmean[fkintran]) + np.random.randn(len(time[intran]))*stdmag
        circularfluxmean = circularfluxmeanlc 
        #plt.plot(fkphase[fkintran],(-circularfluxmean[fkintran]+circularflux[fkintran]-fkdflux/totalFlux)/1.e-6,'r') 
        #plt.xlim([-0.01,0.01])
        #plt.show()
        #plt.plot(fktime,1-circularflux+medianmag,'r',time,mag) 
        #plt.show()
        #dflux = binlc(time[intran],fktime[fkintran],fkdflux)
    else:
        dflux = np.zeros(len(time[intran]))
        obl.relativeFlux(phase[intran],dflux)
            
        z=sma*np.sqrt(np.sin(phase*2*np.pi)*np.sin(phase*2*np.pi)+np.cos(phase*2*np.pi)*np.cos(phase*2*np.pi)*np.cos(inc)*np.cos(inc))
    
        circularflux = occultquad(z,u1,u2,rpol)[0]
    
        index = np.cos(phase*2*np.pi)<0
        circularflux[index]=1.0
    
        circularfluxmean = occultquad(z,u1,u2,rmean)[0]
        circularfluxmean[index]=1.0


        fkmag[intran] = medianmag+np.random.randn(len(fkmag[intran]))*stdmag
        fkmag[intran]+=1-(circularflux[intran]-dflux/totalFlux)
        circularfluxmean=1-circularfluxmean+medianmag+np.random.randn(len(fkmag))*stdmag

    plt.plot(phase[intran],(-circularfluxmean[intran]+fkmag[intran])/1.e-6,'o',mec='b',mfc='None',ms=1.5,mew=1)
    plt.plot(fkphase[fkintran],residual[fkintran],'r')
    plt.show()
    
    return [fkmag,circularfluxmean]

	
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
        #replace with the parsing inside transit class
        transit.P=110.3216229;  
        rpstar=0.08453;
        transit.dip=rpstar**2.
        transit.q=transit.calq(rstar=2.5,logg=4.07)
        print transit.q
        transit.f = 0.1
        transit.inc = 89.209
        transit.alpha = 45.
        transit.u1 = 0.242
        transit.u2 = 0.289
        #
        time=[];readcolumn(time,coljd,infile);time=np.array(time)
        #transit.epoch=min(time)+random.random();
        transit.epoch=1030.36382
        mag=[];readcolumn(mag,colmag,infile);mag=np.array(mag)
        #fkmag,circularflux=trangen(time,mag,transit)
        fkmag,circularflux=trangen(time,mag,transit,lcflag=True)
        outfile=os.path.splitext(infile)[0]+'.fktrap'
        fout=open(outfile,mode='w')
        fout.write('#Generated by lcgen with the following paramters:\n')
        fout.write('#Period = %f, rpstar=%f, Tdur= %f, inc=%f, alpha=%f, f=%f\n' % (transit.P,rpstar,transit.q*transit.P*24.,transit.inc,transit.alpha,transit.f))
        for i in range(len(time)):
            #print i,time[i],mag[i],fkmag[i],circularflux[i]
            fout.write('%13.6f %13.6f %13.6f %13.6f\n' % (time[i],mag[i],fkmag[i],circularflux[i]))
        fout.close()
    return


if __name__=='__main__':
	main()
