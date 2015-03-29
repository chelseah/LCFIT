#!/usr/bin/env python
import math
import numpy as np
import scipy as sp
import os
from dataio import *
from HATlc import lightcurve as lc
from HATlc import tran
import random
import oblatenessfast as Obl
from occultquad import occultquad
#from Eastman_occultquad import occultquad
import matplotlib
from matplotlib import pyplot as plt
import cmd_parse as cmdp
import cfg_parse as cfgp

def fluxtomag(flux):
    flux0 = 1.0
    mag = -2.5*np.log10((flux/flux0))
    #mag = -2.3141*np.log10((flux/flux0))
    return mag

def binlc(filltime,time,mag):
    cadence = 30./60./24.
    nt = len(filltime)
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
    intran=(ftime > (0.5-1.2*q/2.0))*(ftime < (0.5+1.2*q/2.0))
    return intran


def trangen(time,mag,transit,lcflag=0):
    '''functions to generate fake lightcurves'''
    rmean = math.sqrt(transit.dip)
    print rmean,transit.q*transit.P*24.
    f = transit.f
    #inc = transit.inc*math.pi/180.
    #alpha = transit.alpha*np.pi/180.
    transit.calLD()
    u1 = transit.u1
    u2 = transit.u2
    alpha = transit.alpha
    sma = transit.sma
    inc = transit.inc
    req = rmean/math.sqrt(1-f)
    rpol = math.sqrt(1-f)*rmean
    print sma,transit.P,transit.u1,transit.u2,transit.epoch,transit.q,transit.qg,rpol
    ntran=int((max(time)-transit.epoch)/transit.P)-int((min(time)-transit.epoch)/transit.P)+1
    
       
    tranwidth=transit.q*transit.P/2.
    #print req,rpol,alpha,transit.sma,transit.inc,transit.u1,transit.u2
    #obl = Obl.Oblateness(req,rpol,alpha,transit.sma,transit.inc,transit.u1,transit.u2)
    totalFlux = np.pi*(1.0-u1/3-u2/6)
    intran = gentran(time,transit.P,transit.epoch,transit.q)
    medianmag=np.median(mag[-intran])
    transit.stdmag = 0.0012*np.sqrt(2.)
    if transit.stdmag==-1:
        stdmag=np.std(mag[-intran])
    else:
        stdmag=transit.stdmag
    #print stdmag/np.sqrt(2.)
    stdmag/=np.sqrt(2.)
    print stdmag
    #print rpol,medianmag
    phase = (time-transit.epoch)/transit.P-((time-transit.epoch)/transit.P).astype(int)

    ind = phase>0.5
    phase[ind]-=1.0
    fkmag=medianmag+np.random.randn(len(time))*stdmag
    
    if(lcflag==0):
        #initial short cadence data
        fktime = time
        fkmagsc=medianmag+np.random.randn(len(fktime))*stdmag
        #print 'doing sc stuff' 
        #initial short cadence pahse
        fkintran = gentran(fktime,transit.P,transit.epoch,transit.q) 
        fkphase = abs((fktime-transit.epoch)/transit.P)-(abs((fktime-transit.epoch)/transit.P)).astype(int)
        ind = fkphase>0.5
        fkphase[ind]-=1.0
        z=sma*np.sqrt(np.sin(fkphase*2*np.pi)*np.sin(fkphase*2*np.pi)+np.cos(fkphase*2*np.pi)*np.cos(fkphase*2*np.pi)*np.cos(inc)*np.cos(inc))
        
        #compute corrections
        #fkdflux = np.zeros(len(fktime))
        fkdflux = np.zeros(len(fktime[fkintran]))
        obl.relativeFlux(fkphase[fkintran],fkdflux)
        #obl.relativeFlux(fkphase,fkdflux)
        #print (fkdflux==0).all()
        #return
        #print fkphase[fkintran]
        #fkdflux/=totalFlux
        circularfluxmean = occultquad(z,u1,u2,rmean)[0]
        #circularfluxmean = 1-(1-circularfluxmean)/(totalFlux/np.pi)
        circularflux = occultquad(z,u1,u2,rpol)[0]
        circularflux[fkintran]-=fkdflux/totalFlux
        #circularflux = 1-(1-circularflux)/(totalFlux/np.pi)
        index = np.cos(fkphase*2*np.pi)<0
        circularflux[index]=1.0
        circularfluxmean[index]=1.0
        
        #fkmagsc[fkintran] = medianmag
        #fkmagsc[fkintran]+=1-(circularflux[fkintran]-fkdflux)
        #print fkdflux
        circularmag = fluxtomag(circularflux)
        circularmagmean = fluxtomag(circularfluxmean)
        #plt.plot(phase[fkintran],circularmag[fkintran])
        #plt.show()
        fkmag = fkmagsc+circularmag
        circularfluxmean=circularmagmean+fkmagsc
        #circularmagmean=circularmagmean+fkmagsc
        #ax.plot(fkphase[fkintran],residual[fkintran],'r') 
        #ax.plot(fkphase[fkintran],fkmagsc[fkintran],'r') 
        #ax.set_xlim([-0.01,0.01])
        #plt.show()

    else:
        #dflux = np.zeros(len(time[intran]))
        #obl.relativeFlux(phase[intran],dflux)
        #dflux/=totalFlux
        z=sma*np.sqrt(np.sin(phase*2*np.pi)*np.sin(phase*2*np.pi)+np.cos(phase*2*np.pi)*np.cos(phase*2*np.pi)*np.cos(inc)*np.cos(inc))
    
        circularflux = occultquad(z,u1,u2,rpol)[0]
    
        index = np.cos(phase*2*np.pi)<0
        circularflux[index]=1.0
    
        circularfluxmean,check = occultquad(z,u1,u2,rmean)
        circularfluxmean[index]=1.0
        
        #circularflux = 1-(1-circularflux)/(totalFlux/np.pi)
        #circularfluxmean = 1-(1-circularfluxmean)/(totalFlux/np.pi)
        #ax.plot(phase[intran],circularfluxmean[intran],'r')
        #ax.plot(phase[intran],(circularfluxmean[intran]-circularflux[intran]+dflux)/1.e-6,'r')
        #ax.set_xlim([-0.01,0.01])
        #plt.show()
        fkmag = mag
        print len(circularflux[intran])
        #fkmag[intran] = medianmag+np.random.randn(len(fkmag[intran]))*stdmag
        #fkmag[intran]+=1-(circularflux[intran]-dflux)
        fkmag[intran]+=1-(circularflux[intran])
        circularfluxmean=1-circularfluxmean+medianmag
        #+np.random.randn(len(fkmag))*stdmag
    #ax.plot(phase[intran],(-circularfluxmean[intran]+fkmag[intran])/1.e-6,'o',mec='b',mfc='None',ms=1.5,mew=1)
    #ax.set_xlabel('phase')
    #ax.set_ylabel('residual(ppm)')
    #plt.show()
    #fig=plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(phase[intran]*transit.P*24.,circularfluxmean[intran],'o',mec='b',mfc='None',ms=1.5,mew=1)
    #ax.plot(phase[intran]*transit.P*24.,fkmag[intran],'.',mec='r')
    #ax.plot(phase[intran],-circularfluxmean[intran]+mag[intran],'o',mec='b',mfc='None',ms=1.5,mew=1)
    #ax.set_xlabel('phase')
    #ax.set_ylabel('relativemag(ppm)')
    #yformatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    #ax.yaxis.set_major_formatter(yformatter)
    #plt.show()
    #fig=plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(phase[intran],-circularfluxmean[intran]+mag[intran]+medianmag+1-(circularflux[intran]-dflux),'o',mec='b',mfc='None',ms=1.5,mew=1)
    #ax.set_xlabel('phase')
    #ax.set_ylabel('relativemag(ppm)')
    #plt.show()
     
    return [phase,fkmag,circularfluxmean]

	
def main():
    lcfile=lc()
    transit=tran()
    options = cmdp.ltf_parse()
    infileflag = options.infileflag
    outfileflag = options.outfileflag
    inpathflag = options.inpathflag
    noplot=options.noplot
    cfgfile = options.cfg
    uflag = options.uflag
    print cfgfile
    lcflag = int(cfgp.File_parse(cfgfile,'lcflag'))
    infile = cfgp.File_parse(cfgfile,'infile')
    inlist = cfgp.File_parse(cfgfile,'inlist')
    inpath = cfgp.File_parse(cfgfile,'inpath')
    outfile = cfgp.File_parse(cfgfile,'outfile')
    coljd= int(cfgp.File_parse(cfgfile,'coljd'))
    colmag= int(cfgp.File_parse(cfgfile,'colmag'))
    transit.readpara(cfgfile)
    transit.inc /= (180./np.pi)
    transit.alpha /= (180./np.pi)
    names = []; readcolumn(names,1,inlist,datformat='str')
    if not inpath=='':
        os.chdir(inpath)
        
    lcfile = lc()
    if lcflag==1:
        lcfile.cadence = 'long'
    else:
        lcfile.cadence = 'short'
    print names
    for x in names:
                                             
        time=[];readcolumn(time,coljd,x);time=np.array(time)
        mag=[];readcolumn(mag,colmag,x);mag=np.array(mag)
        phase,fkmag,circularflux=trangen(time,mag,transit,lcflag=lcflag)
        outfile=os.path.splitext(x)[0]+'.fktrap'
        print outfile
        fout=open(outfile,mode='w')
        fout.write('#Generated by lcgen with the following paramters:\n')
        fout.write('#Period = %f, rpstar=%f, Tdur= %f, inc=%f, alpha=%f, f=%f\n' % (transit.P,math.sqrt(transit.dip),transit.q*transit.P*24.,transit.inc,transit.alpha,transit.f))
        for i in range(len(time)):
            fout.write('%13.6f %13.6f %13.6f %13.6f %13.6f %13.6f\n' % (time[i],phase[i],mag[i],fkmag[i],fkmag[i]-np.median(fkmag),circularflux[i]))
        fout.close()
    return


if __name__=='__main__':
	main()
