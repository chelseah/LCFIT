#!/usr/bin/env python

from dataio import *
from lslft_recon import *

def main():
    '''Given the position of transit to window (given by tran,tran[0]=ts, tran[1]=te) and do a better fit for the base line\n
	  Only use when reprocess the selected lightcurves, apply for a short baseline with one deep known transit\n
	  It is a useful thing to do before the kepfit_simple.sh
	  Warning: the output file does not containing the full information of the inputfile, only 2 columns with BJD and detrended flux. The replaceflag is default off to prevent unwantted replace of files. '''

    infile = 'test.txt'
    outfile = 'test_recon.txt'
    coljd = 1
    colmag = 4
    time=[]
    mag=[]
    #n=60
    #tmin=10.0
    tmin=0.2
    qn=1
    P = 125.6306458; epoch = 1086.33958; Tdur = 9.0432; q = Tdur/P/24.
    readcolumn(time,coljd,infile)
    readcolumn(mag,colmag,infile)
    time=np.array(time)
    mag=np.array(mag)
    print len(time),len(mag),max(time),min(time)
    intran = gentran(time,P,epoch,q)
    #dip = np.mean(mag[intran])-np.mean(mag[-intran])
    #dip = 0.00775
    dip = 0.07871**2. 
    if(dip<0):
        dip=0.01
    #print mag[intran]
    print np.mean(mag),np.mean(mag[intran]),np.mean(mag[-intran]),dip
    #time-=2454000
    n=round((max(time)-min(time))/tmin)
    #tran = np.zeros(len(time))
    #intran = (1719.82<time)*(time<1720.08)
    #intran = (2086.46<time)*(time<2086.79)
    #print time[intran]
    #print len(time[intran])
    #dflux=lsfitrecon(time,mag,intran,n=n,qn=qn,wn=7,dipguess=dip)
    dflux=lssubrecon(time,mag,intran,n=n,qn=qn,wn=7)
    dflux=dflux+(min(mag)-np.median(dflux))
    if((not os.path.exists(outfile)) or replace):	
        fout=open(outfile,mode='w')
        for i in range(len(time)):
            #line='%10.6f %10.6f\n' % (time[i]+2454000,dflux[i])
            line='%10.6f %10.6f\n' % (time[i],dflux[i])
            fout.write(line)
        fout.close()
    else:
        print '%s already exists, turn replaceflag on if want to replace it' % (outfile)
    return	
	
if __name__=='__main__':
    main()
