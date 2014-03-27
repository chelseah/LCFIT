#!/usr/bin/env python
from refold import *
from dataio import *
import numpy as np
import scipy as sp
import random
from dataio import *
from HATlc import lightcurve as lc
from HATlc import tran
from lcgen import gentran,binlc
import oblatenessfast as Obl
from occultquad import occultquad
import matplotlib
from matplotlib import pyplot as plt
from math import sqrt,cos,sin
import sys
import cmd_parse as cmdp
import cfg_parse as cfgp

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

def windowphase(time,mag,epoch,period,q,ph,wd=8):
    ftime=phasefold(time,epoch,period)
    newtime=time[abs(ftime-ph)<wd*q]
    newmag=mag[abs(ftime-ph)<wd*q]
    newftime=ftime[abs(ftime-ph)<wd*q]
    
    return [newtime,newmag,newftime]

def simpletemp():
    infile = 'HAT-266-0008524.tfalc'
    time = []; readcolumn(time,2,infile); time = np.array(time)
    mag = []; readcolumn(mag,18,infile); mag = np.array(mag)
    period = 9.199166; epoch = 54774.2606359; q = 0.0243
    epoch+=period*q/2.
    rp = math.sqrt(0.0122)
    inc = math.pi/2.
    phase = abs((time-epoch)/period)-(abs((time-epoch)/period)).astype(int)
    ind = phase>0.5
    phase[ind]-=1.0
    sma = 1./(q*math.pi)
    #phase = phasefold(time,epoch,period)
    z=sma*np.sqrt(np.sin(phase*2*math.pi)*np.sin(phase*2*math.pi)+np.cos(phase*2*math.pi)*np.cos(phase*2*math.pi)*np.cos(inc)*np.cos(inc))
    circularflux = occultquad(z,0.0,0.0,rp)[0]
    index = np.cos(phase*2*np.pi)<0
    circularflux[index]=1.0
    circularmag = fluxtomag(circularflux)

    for i in xrange(len(time)):
        print time[i],phase[i],circularmag[i]+np.median(mag),mag[i]
    return 

def simple():
    inpath='data/'
    os.chdir(inpath)
    lcfile = lc()
    transit =tran()
    lcfile.cadence = 'short'
    #lcfile.cadence = 'long'
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
    lcfile.name='kplr011391018-2010009094841_slc.tab'
    #infile='kplr011391018.tab_recon.7'
    #lcfile.name = 'kplr007906882-2011073133259_slc.tab'
    #lcfile.name = 'HAT-266-0008524.tfalc'
    #lcfile.name = 'kplr007906882-2013098041711_llc.pflt'
    outfile=''
    lcfile.coljd=1
    lcfile.colmag=2
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
    transit.epoch= 981.0908934987607
    transit.P= 30.3604338541 
    transit.rmean = 0.130681438524 
    transit.rsum = 0.0274694247696
    transit.inc = 89.5179576467/180.*np.pi
    transit.sma = (1+transit.rmean)/transit.rsum
    transit.q = 1./transit.sma/np.pi
    transit.b = np.cos(transit.inc)*transit.sma
    #b = 0.675
    #qvar = 3.2639/period/24.
    #sma = 1./(qvar*np.pi)
    #inc = np.arccos(b/sma)
    transit.q1 = 0.541139148419 
    transit.q2 = 0.402402665

    transit.calLD()
    #u1 = 0.347
    #u2 = 0.258
    print transit.b, transit.sma,transit.u1,transit.u2,transit.inc/np.pi*180.
    #u1 = 0.242 
    #u2 = 0.289
    #f = 0.1
    #alpha =45./180.*np.pi
    #f = 0.058
    #alpha =11.22/180.*np.pi
    transit.f = 0.089
    #transit.alpha =45./180.*np.pi
    #transit.alpha =-46.379/180.*np.pi
    transit.alpha =-34.9191359645/180.*np.pi
    #GenLC(lcfile,transit)
    #return
    fragtime,fragmag,residual,fragftime,fragmagcircular = GenLC(lcfile,transit)
    if (outfile==''):
        outfile=os.path.splitext(lcfile.name)[0]+'.rmtran' 
    print outfile
    ftfan=open(outfile,mode='w')
    for i in range(len(fragtime)):
        ftfan.write('%-14.7f %-14.7f %14.7f %14.7f %14.7f\n' % (fragtime[i],fragmag[i],residual[i],fragftime[i],fragmagcircular[i]))
    ftfan.close()				
    return	


def GenLC(lcfile,transit):
    req = transit.rmean/sqrt(1-transit.f)
    rpol = sqrt(1-transit.f)*transit.rmean
    totalFlux = np.pi*(1.0-transit.u1/3-transit.u2/6)
    if lcfile.cadence =='short':
        wn = 24.*60.
        minlen=transit.q*transit.P*24.*60.*4.
    else:
        wn = 24.*2.
        minlen=transit.q*transit.P*24.*2.*4.

    time,mag=readcol(lcfile.name,lcfile.coljd,lcfile.colmag)
    time=np.array(time)
    mag=np.array(mag)
    #rm transit
    newtime,newmag,newftime=windowtran(time,mag,transit.epoch,transit.P,transit.q)
    fragtime=[]
    print minlen,transit.q,transit.P,lcfile.cadence
    #return
    #generate random segment of lc 
    ideal = False
    if(ideal):
        fragtime = time
        fragmag = np.zeros(len(mag))+np.median(newmag)
        fragftime = phasefold(time,transit.epoch,transit.P)
        ph=0.5
        #fragtime,fragmag,fragftime=windowphase(newtime,newmag,epoch,period,qvar,ph)
    else:
        while len(fragtime)<minlen:
            ph = random.random()
            print ph,transit.q
            while (abs(ph-0.5)<2.*transit.q):
                ph = random.random()
            fragtime,fragmag,fragftime=windowphase(newtime,newmag,transit.epoch,transit.P,transit.q,ph)
    newepoch = transit.epoch+(ph-0.5)*transit.P
    #print newepoch
    if lcfile.cadence =='short':
        fkintran = gentran(fragtime,transit.P,newepoch,transit.q)
        fkphase = abs((fragtime-newepoch)/transit.P)-(abs((fragtime-newepoch)/transit.P)).astype(int)
        ind = fkphase>0.5
        fkphase[ind]-=1.0
        #print fkphase
        z = transit.calZ(fkphase) 
        #compute corrections
        fkdflux = np.zeros(len(fragtime[fkintran]))
        #print req,rpol,transit.alpha,transit.sma,transit.inc,transit.u1,transit.u2
        obl = Obl.Oblateness(req,rpol,transit.alpha,transit.sma,transit.inc,transit.u1,transit.u2)
        obl.relativeFlux(fkphase[fkintran],fkdflux)
        #print fkphase[fkintran] 
        circularfluxmean = occultquad(z,transit.u1,transit.u2,transit.rmean)[0]
        #plt.plot(fkphase,z)
        #plt.show()
        circularflux = occultquad(z,transit.u1,transit.u2,rpol)[0]
        index = np.cos(fkphase*2*np.pi)<0
        circularflux[index]=1.0
        circularflux[fkintran]-=fkdflux/totalFlux
        circularfluxmean[index]=1.0
        circularmag = fluxtomag(circularflux)
        circularmagmean = fluxtomag(circularfluxmean)
        fragmagcircular=fragmag+circularmagmean
        fragmag+=circularmag
        residual = (fragmag-fragmagcircular)/1.e-6
    elif lcfile.cadence=='long':
        npoints = int(round((max(fragtime)-min(fragtime))*24.*60.))
        fragtimesc = np.arange(npoints)*(1./24./60.)+min(fragtime) 
        fkintransc = gentran(fragtimesc,transit.P,newepoch,transit.q) 
        fkphasesc = abs((fragtimesc-newepoch)/transit.P)-(abs((fragtimesc-newepoch)/transit.P)).astype(int)
        ind = fkphasesc>0.5
        fkphasesc[ind]-=1.0
        z = transit.calZ(fkphasesc)
        print max(fragtime), min(fragtime) 
        #compute corrections
        fkdfluxsc = np.zeros(len(fragtimesc[fkintransc]))
        obl = Obl.Oblateness(req,rpol,transit.alpha,transit.sma,transit.inc,transit.u1,transit.u2)
        obl.relativeFlux(fkphasesc[fkintransc],fkdfluxsc)
        circularfluxmeansc = occultquad(z,transit.u1,transit.u2,transit.rmean)[0]
        circularfluxsc = occultquad(z,transit.u1,transit.u2,rpol)[0]
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
    fragtime -= (newepoch-transit.epoch)
    return [fragtime,fragmag,residual,fragftime,fragmagcircular] 
    

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
    inlist = cfgp.File_parse(cfgfile,'inlist')
    inpath = cfgp.File_parse(cfgfile,'inpath')
    outfile = cfgp.File_parse(cfgfile,'outfile')
    coljd= int(cfgp.File_parse(cfgfile,'coljd'))
    colmag= int(cfgp.File_parse(cfgfile,'colmag'))
    transit =tran()
    transit.readpara(cfgfile)
    transit.calLD()
    transit.inc /= (180./np.pi)
    transit.alpha /= (180./np.pi)
    names = []; readcolumn(names,1,inlist,datformat='str')
    print transit.b, transit.sma,transit.u1,transit.u2,transit.inc/np.pi*180.,transit.alpha/np.pi*180.
    if not inpath=='':
        os.chdir(inpath)
    
    lcfile = lc()
    if lcflag==1:
        lcfile.cadence = 'long'
    else:
        lcfile.cadence = 'short'
    for x in names:
        lcfile.name = x
        lcfile.coljd = coljd
        lcfile.colmag = colmag
        GenLC(lcfile,transit)
        fragtime,fragmag,residual,fragftime,fragmagcircular = GenLC(lcfile,transit)
        #break
        #if (outfile==''):
        outfile=os.path.splitext(lcfile.name)[0]+'.rmtran' 
        print outfile
        ftfan=open(outfile,mode='w')
        for i in range(len(fragtime)):
            ftfan.write('%-14.7f %-14.7f %14.7f %14.7f %14.7f\n' % (fragtime[i],fragmag[i],residual[i],fragftime[i],fragmagcircular[i]))
        ftfan.close()				
    return	



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
    main()
    #simple()
    #simpletemp()
