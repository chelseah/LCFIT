#!/usr/bin/env python
from refold import *
from dataio import *
import numpy as np
import scipy as sp
import random
from lcgen import gentran,binlc
import oblatenessfast as Obl
from occultquad import occultquad
import matplotlib
from matplotlib import pyplot as plt
from math import sqrt,cos,sin
import sys

def fluxtomag(flux):
    flux0 = 1.0
    mag = -2.5*np.log10((flux/flux0))
    #mag = -2.3141*np.log10((flux/flux0))
    return mag
def purmute(mag,wn):
    length=len(mag)
    mag=np.array(mag)
    newmag=np.array_split(mag,int(length/wn))
    for x in newmag:
        random.shuffle(x)
    #for i in newmag:
    return np.concatenate(newmag)

def windowtran(time,mag,epoch,period,q,wd=0.6):
    ftime=phasefold(time,epoch,period)
    newtime=time[abs(ftime-0.5)>wd*q]
    newmag=mag[abs(ftime-0.5)>wd*q]
    newftime=ftime[abs(ftime-0.5)>wd*q]
    
    return [newtime,newmag,newftime]

def windowphase(time,mag,epoch,period,q,ph,wd=4):
    ftime=phasefold(time,epoch,period)
    newtime=time[abs(ftime-ph)<wd*q]
    newmag=mag[abs(ftime-ph)<wd*q]
    newftime=ftime[abs(ftime-ph)<wd*q]
    
    return [newtime,newmag,newftime]

def simple():
    inpath='data/'
    #cadence = 'short'
    cadence = 'long'
    os.chdir(inpath)
    #infile='kplr006603043-2011024051157_slc.tab'
    #infile='kplr006603043-2013011073258_llc.tab'
    #infile='kplr006603043-2012277125453_llc.tab'
    #infile='kplr006603043-2012179063303_llc.pflt'
    #infile='kplr006603043-2012088054726_llc.pflt'
    #infile='kplr006603043-2012004120508_llc.pflt'
    #infile = 'kplr006603043-2011271113734_llc.pflt'
    #infile = 'kplr006603043-2011177032512_llc.pflt'
    #infile = 'kplr006603043-2011073133259_llc.pflt'
    #infile = 'kplr006603043-2009166043257_llc.pflt'   
    #infile = 'kplr006603043-2009259160929_llc.pflt'
    #infile = 'kplr006603043-2009350155506_llc.pflt'
    #infile = 'kplr006603043-2010078095331_llc.pflt'
    #infile = 'kplr006603043-2010174085026_llc.pflt'
    #infile = 'kplr006603043-2010265121752_llc.pflt'
    #infile = 'kplr006603043-2010355172524_llc.pflt'
    #infile='kplr011391018.ltf.0'
    #infile='kplr011391018.tab_recon.7'
    #infile = 'kplr007906882-2011073133259_slc.tab'
    infile = 'kplr007906882-2013098041711_llc.pflt'
    outfile=''
    coljd=1
    colmag=2
    #KIC1139108
    #epoch=981.091128085308 
    #period= 30.3604520162
    #rmean = 0.134185445174
    #rsum = 0.0274010731669
    #inc = 89.5034018837/180.*np.pi
    #sma = (1+rmean)/rsum
    #qvar = 1./sma/np.pi
    #q1 = 0.439524664543
    #q2 = 0.476624913031
    #b = np.cos(inc)
    #u1 = 2.*np.sqrt(q1)*q2 
    #u2 = np.sqrt(q1)*(1-2*q2)
    #KIC7906882
    #0.248946023677 0.388995838747 52.5135743015 0.673667358282 0.0108005175588 0.120193698155 89.5954361728
    epoch=833+171.67373 
    period= 52.513565
    rmean = 0.120193698155
    rsum = 0.0108005175588
    inc = 89.5954361728/180.*np.pi
    sma = (1+rmean)/rsum
    qvar = 1./sma/np.pi
    b = np.cos(inc)*sma
    #b = 0.675
    #qvar = 3.2639/period/24.
    #sma = 1./(qvar*np.pi)
    #inc = np.arccos(b/sma)
    q1 = 0.248946023677
    q2 = 0.388995838747

    u1 = 2.*np.sqrt(q1)*q2 
    u2 = np.sqrt(q1)*(1-2*q2)
    #u1 = 0.347
    #u2 = 0.258
    print b, sma,u1,u2,inc/np.pi*180.
    #u1 = 0.242 
    #u2 = 0.289
    #f = 0.1
    #alpha =45./180.*np.pi
    #f = 0.058
    #alpha =11.22/180.*np.pi
    f = 0.0
    alpha =0./180.*np.pi
    req = rmean/sqrt(1-f)
    rpol = sqrt(1-f)*rmean
    totalFlux = np.pi*(1.0-u1/3-u2/6)
    if cadence =='short':
        wn = 24.*60.
        minlen=qvar*period*24.*60.*4.
    else:
        wn = 24.*2.
        minlen=qvar*period*24.*2.*4.

    time,mag=readcol(infile,coljd,colmag)
    time=np.array(time)
    mag=np.array(mag)
    #rm transit
    newtime,newmag,newftime=windowtran(time,mag,epoch,period,qvar)
    fragtime=[]
    print minlen
    #return
    #generate random segment of lc 
    ideal = False
    if(ideal):
        fragtime = time
        fragmag = np.zeros(len(mag))+np.median(newmag)
        fragftime = phasefold(time,epoch,period)
        ph=0.5
        #fragtime,fragmag,fragftime=windowphase(newtime,newmag,epoch,period,qvar,ph)
    else:
        while len(fragtime)<minlen:
            ph = random.random()
            print ph,qvar
            while (abs(ph-0.5)<2.*qvar):
                ph = random.random()
            fragtime,fragmag,fragftime=windowphase(newtime,newmag,epoch,period,qvar,ph)
    newepoch = epoch+(ph-0.5)*period
    #print newepoch
    if cadence =='short':
        fkintran = gentran(fragtime,period,newepoch,qvar)
        fkphase = abs((fragtime-newepoch)/period)-(abs((fragtime-newepoch)/period)).astype(int)
        ind = fkphase>0.5
        fkphase[ind]-=1.0
        #print fkphase
        z=sma*np.sqrt(np.sin(fkphase*2*np.pi)*np.sin(fkphase*2*np.pi)+np.cos(fkphase*2*np.pi)*np.cos(fkphase*2*np.pi)*np.cos(inc)*np.cos(inc))
        
        #compute corrections
        fkdflux = np.zeros(len(fragtime[fkintran]))
        obl = Obl.Oblateness(req,rpol,alpha,sma,inc,u1,u2)
        obl.relativeFlux(fkphase[fkintran],fkdflux)
        
        circularfluxmean = occultquad(z,u1,u2,rmean)[0]
        #plt.plot(fkphase,z)
        #plt.show()
        circularflux = occultquad(z,u1,u2,rpol)[0]
        index = np.cos(fkphase*2*np.pi)<0
        circularflux[index]=1.0
        circularflux[fkintran]-=fkdflux/totalFlux
        circularfluxmean[index]=1.0
        circularmag = fluxtomag(circularflux)
        circularmagmean = fluxtomag(circularfluxmean)
        fragmagcircular=fragmag+circularmagmean
        fragmag+=circularmag
        residual = (fragmag-fragmagcircular)/1.e-6
    elif cadence=='long':
        npoints = int(round((max(fragtime)-min(fragtime))*24.*60.))
        fragtimesc = np.arange(npoints)*(1./24./60.)+min(fragtime) 
        fkintransc = gentran(fragtimesc,period,newepoch,qvar) 
        fkphasesc = abs((fragtimesc-newepoch)/period)-(abs((fragtimesc-newepoch)/period)).astype(int)
        ind = fkphasesc>0.5
        fkphasesc[ind]-=1.0
        z=sma*np.sqrt(np.sin(fkphasesc*2*np.pi)*np.sin(fkphasesc*2*np.pi)+np.cos(fkphasesc*2*np.pi)*np.cos(fkphasesc*2*np.pi)*np.cos(inc)*np.cos(inc))
        print max(fragtime), min(fragtime) 
        #compute corrections
        fkdfluxsc = np.zeros(len(fragtimesc[fkintransc]))
        obl = Obl.Oblateness(req,rpol,alpha,sma,inc,u1,u2)
        obl.relativeFlux(fkphasesc[fkintransc],fkdfluxsc)
        circularfluxmeansc = occultquad(z,u1,u2,rmean)[0]
        circularfluxsc = occultquad(z,u1,u2,rpol)[0]
        index = np.cos(fkphasesc*2*np.pi)<0
        circularfluxsc[index]=1.0
        circularfluxmeansc[index]=1.0
        circularfluxsc[fkintransc]-=fkdfluxsc/totalFlux
        circularmagsc = fluxtomag(circularfluxsc)
        circularflux=binlc(fragtime,fragtimesc,circularmagsc)
        circularmagmean = fluxtomag(circularfluxmeansc)
        circularfluxmean=binlc(fragtime,fragtimesc,circularmagmean)
        fragmagcircular=fragmag+circularfluxmean
        fragmag+=circularflux
        residual = (fragmag-fragmagcircular)/1.e-6
    fragtime -= (newepoch-epoch)
        
    if (outfile==''):
        outfile=os.path.splitext(infile)[0]+'.rmtran' 
    print outfile
    ftfan=open(outfile,mode='w')
    for i in range(len(fragtime)):
        ftfan.write('%-14.7f %-14.7f %14.7f %14.7f %14.7f\n' % (fragtime[i],fragmag[i],residual[i],fragftime[i],fragmagcircular[i]))
    ftfan.close()				
		
    return


def testsimple():
    inpath='data'
    cadence = 'short'
    if cadence =='short':
        wn = 24.*60.
    else:
        wn = 24.*2.
    os.chdir(inpath)
    infile='kplr006603043-2011024051157_slc.tab'
    outfile=''
    coljd=1
    colmag=2
    epoch=1030.3645
    period= 110.3216229
    qvar= 13.32/period/24.
    #epoch+=period*qvar*0.5
    time,mag=readcol(infile,coljd,colmag)
    time=np.array(time)
    mag=np.array(mag)
    newtime,newmag,newftime=windowtran(time,mag,epoch,period,qvar)
    fkmag = purmute(newmag,wn)
    if (outfile==''):
        outfile=os.path.splitext(infile)[0]+'.rmtran' 
    print outfile
    ftfan=open(outfile,mode='w')
    for i in range(len(newtime)):
        ftfan.write('%-14.7f %-14.7f %14.7f\n' % (newtime[i],fkmag[i],newftime[i]))
    ftfan.close()				
		
    return

if __name__=='__main__':
    sys.argv[0]='lcfromdata'
    simple()
