#!/usr/bin/env python
from refold import *
from dataio import *
import numpy as np
import scipy as sp
import random
from lcgen import gentran
from occultquad import occultquad
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
    inpath='data'
    cadence = 'short'
    os.chdir(inpath)
    infile='kplr006603043-2011024051157_slc.tab'
    outfile=''
    coljd=1
    colmag=2
    epoch=2455030.3645
    period= 110.3216229
    rmean = 0.08453
    qvar= 13.32/period/24.
    inc = 90./180.*np.pi
    sma = 51.1812
    u1 = 0.242 
    u2 = 0.289
    #epoch+=period*qvar*0.5
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
    #generate random segment of lc 
    while len(fragtime)<minlen:
        ph = random.random()
        while (abs(ph-0.5)<2.*qvar):
            ph = random.random()
        print ph
        fragtime,fragmag,fragftime=windowphase(newtime,newmag,epoch,period,qvar,ph)
   
    #circularfluxmeanlc=medianmag+np.random.randn(len(time))*stdmag
    #cadence = 1./60./24.
    #nt = (int)((max(time)-min(time))/cadence)
    #fktime = min(time)+np.arange(nt)*cadence
    #fkmagsc=medianmag+np.random.randn(len(fktime))*stdmag
    #print 'doing sc stuff' 
    #initial short cadence pahse
    newepoch = epoch+(ph-0.5)*period
    fkintran = gentran(fragtime,period,newepoch,qvar) 
    fkphase = abs((fragtime-newepoch)/period)-(abs((fragtime-newepoch)/period)).astype(int)
    ind = fkphase>0.5
    fkphase[ind]-=1.0
    z=sma*np.sqrt(np.sin(fkphase*2*np.pi)*np.sin(fkphase*2*np.pi)+np.cos(fkphase*2*np.pi)*np.cos(fkphase*2*np.pi)*np.cos(inc)*np.cos(inc))
    
    #compute corrections
    #fkdflux = np.zeros(len(fktime[fkintran]))
    #obl.relativeFlux(fkphase[fkintran],fkdflux)
    ##fkdflux/=totalFlux
    circularfluxmean = occultquad(z,u1,u2,rmean)[0]
    #circularfluxmean = 1-(1-circularfluxmean)/(totalFlux/np.pi)
    #circularflux = occultquad(z,u1,u2,rpol)[0]
    #circularflux = 1-(1-circularflux)/(totalFlux/np.pi)
    index = np.cos(fkphase*2*np.pi)<0
    #circularflux[index]=1.0
    circularfluxmean[index]=1.0
    
    #fkmagsc[fkintran]+=1-(circularflux[fkintran]-fkdflux)
    fragmag+=1-circularfluxmean
    #fkmag[intran] = binlc(time[intran],fktime[fkintran],fkmagsc[fkintran])
    #residual = (fkmagsc-circularfluxmean)/1.e-6
    #circularfluxmeanlc[intran] = binlc(time[intran],fktime[fkintran],circularfluxmean[fkintran]) + np.random.randn(len(time[intran]))*stdmag
    #circularfluxmean = circularfluxmeanlc 

        
    if (outfile==''):
        outfile=os.path.splitext(infile)[0]+'.rmtran' 
    print outfile
    ftfan=open(outfile,mode='w')
    for i in range(len(fragtime)):
        ftfan.write('%-14.7f %-14.7f %14.7f\n' % (fragtime[i],fragmag[i],fragftime[i]))
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
    epoch=2455030.3645
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
	simple()
