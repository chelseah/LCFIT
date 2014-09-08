#!/usr/bin/env python
import os
import math
import pyfits
import shutil
import numpy as np
from info import *
import cmd_parse as cmdp
import cfg_parse as cfgp
from dataio import *

# this function read Kepler FITS file
def openfit(infile,outfile='',outpath='',field=[],magflag=True,momflag=False,BJD=-1,z0=0,fluxflag=False):
    '''function openfit: extract Kepler FITS file into ASCII files
    Default setting: output .tab file including BJD,SAPMAG,SAPMAGERR, PDCMAG,PDCMAGERR
    Default format: BJD-2454000, Magnitude with z0=25.2415
    Optional output: momentum centriod files .mom including BJD,MOMX,MOMXERR,MOMY,MOMYERR\n
    Example for default input\n
    openfit(infile,outfile='',field=[],magflag=True,momflag=False,BJD=-1,z0=0)\n    Warning: please turn off magflag and momflag if want to specify output fields\n'''	
    if (not outpath==''):
        if(not os.path.isdir(outpath)):
            print 'OSError: No such file or directory: %s for the output path' % (outpath)
            exit(2)

#	print len(field),field
    if(infile.endswith('.fits')):
        print infile
        hdulist=pyfits.open(infile)
        tbdata=hdulist[1].data
        if(BJD==-1):
            BJD=2454000
        if(z0==0):
            z0=25.2415
        BJDREFI=hdulist[1].header['BJDREFI']
        BJDREFF=hdulist[1].header['BJDREFF']
        time=tbdata.field('TIME')+BJDREFI+BJDREFF-BJD
        base=os.path.basename(infile)
        if (outfile==''):
            outfile=os.path.splitext(base)[0]+'.tab'
            outfile2=os.path.splitext(base)[0]+'.mom'
        outfile=outpath+outfile
        outfile2=outpath+outfile2
        
        if (magflag):
            flux=tbdata.field('SAP_FLUX')
            fluxerr=tbdata.field('SAP_FLUX_ERR')
            cflux=tbdata.field('PDCSAP_FLUX')
            cfluxe=tbdata.field('PDCSAP_FLUX_ERR')
            fout=open(outfile,mode='w')
            print 'Extrant mag to .tab file: %s\n' % (outfile)
            head='# BJD-2454000 Mag Magerr KcMag KcMagerr\n'
            fout.write(head)
            if(fluxflag==True):
                mag=flux
                mage = fluxerr
                cmag=cflux
                cmage=cfluxe
            else:
                mag=z0-2.5*np.log10(np.array(flux))
                mage=abs(2.5/np.log(10)*fluxerr/flux)
                cmag=z0-2.5*np.log10(np.array(cflux))
                cmage=abs(2.5/np.log(10)*cfluxe/cflux)
            for i in range(len(time)):
                if ((not math.isnan(mag[i])) and (not math.isnan(cmag[i]))):
                    line='%.13f %.13f %.13e %.13f %.13e\n' % (time[i],mag[i],mage[i],cmag[i],cmage[i])
                    fout.write(line)
            fout.close()
        
        if (momflag):
            momx=tbdata.field('MOM_CENTR1')
            momxe=tbdata.field('MOM_CENTR1_ERR')
            momy=tbdata.field('MOM_CENTR2')
            momye=tbdata.field('MOM_CENTR2_ERR')
            fout2=open(outfile2,mode='w')
            print 'Extrant momentum centriods to .mom file: %s\n' % (outfile)
            head2='# BJD-2454000 Momx Momxerr Momy Momyerr\n'
            fout2.write(head2)
            for i in range(len(time)):
                if ((not math.isnan(momx[i])) and (not math.isnan(momy[i]))):
                    line2='%.13f %.13f %.13e %.13f %.13e\n' % (time[i],momx[i],momxe[i],momy[i],momye[i])
                    fout2.write(line2)
            fout2.close()
        
        if(not len(field)==0):
            varchar=[]
            fout=open(outfile,mode='w')
            print 'Extract Field: BJD %s\n' % (' '.join(field))
            for x in field: 	
                name=str(x)
                var=tbdata.field(name)
                varchar.append(var)
                nf=len(field)
            head='#'+' '.join(field)+'\n'
            fout.write(head)
            for i in range (len(time)):
                fout.write('%.13f' % (time[i]))			
                for j in range(nf):
                    fout.write(' %13.7f' % (varchar[j][i]))
                fout.write('\n')
            
        hdulist.close()
        
    else:
        errmes='INPUTERROR: -- openfit: %s is not a FITS file' % (infile)
        print errmes
        exit(1)
    return


########################################################################

def main():
    '''Module Readfits:
	   For usage and options: use python readfits --help (-h)
	   Make sure a configuration file is in the directory or provided by the comand line.
	   If output path (-o or -O) is not provided, the default is the directory where the input file is;
	   If output file is not provided (in configuration file), the default is the basename of the input file with extension '.tab' (if output default magnitude files or user choice of fileds) or '.mom' (if output default centroids files)
	   If output format is not provided ('field' in configuration file), there are two choices of standard output format, magnitude file and centroids file, depend on the magflag and momflag in the configuration file. Default setting is magflag=True and momflag=False. If both are True, routine output two files for the same input FITS files. For more information, read the openfits section of help. 	
	   If conflicted options are provided: the configuration file always have higher priority than the conmand line option.
	   Other cases of confliction: 
	   if a infile is provided (either through -i or -I) the routine only works on that single file, the inlist options are ignored and the inpath option be treated as where the file is from .
	   if a infile is not provided, and a input list is provided (either through -l or -L), the routine works on the files provided as input list, any option on outfile is ignored.
	   if neither a infile and a input list is provided, but a input path is provided (either through -p or -P), the routine interprete as reading all the FITS files in the provided directory, any option on outfile is ignored.'''

#Parse command line
    options=cmdp.readfits_parse()
    
    inpathflag=options.inpathflag
    inlistflag=options.inlistflag
    infileflag=options.infileflag
    outpathflag=options.outpathflag
    fluxflag=options.nomagflag
    cfgfile=options.cfg
    uflag=options.uflag	
#Parse general configuration
    BJD=float(cfgp.gen_parse(cfgfile,'BJD'))
    z0=float(cfgp.gen_parse(cfgfile,'z0'))	
    #handelfile='/home/chelsea/src/newpipe/missreport'
    #fout=open(handelfile,mode='w')
#Parse configuration
    if(inpathflag):
        inpath=cfgp.readfits_parse(cfgfile,'inpath')
    else:
        inpath=options.inpath
    
    if(inlistflag):
        inlist=cfgp.readfits_parse(cfgfile,'inlist')
    else:
        inlist=options.inlist
    if(infileflag):
        infile=cfgp.readfits_parse(cfgfile,'infile')
    else:
        infile=options.infile
    
    if(outpathflag):
        outpath=cfgp.readfits_parse(cfgfile,'outpath')
    else:
        outpath=options.outpath	
    #always readin outputfile, the default is using what the default setting gives
    outfile=cfgp.readfits_parse(cfgfile,'outfile')

    if(uflag):
        print 'Parsing user supplied parameters from %s\n' % (cfgfile)
        rawfield=cfgp.readfits_parse(cfgfile,'field')
        if(not len(rawfield)==0):
            field=rawfield.split(',')
        else:
            field=[]
        
        momflag=bool(int(cfgp.readfits_parse(cfgfile,'momflag')))
        magflag=bool(int(cfgp.readfits_parse(cfgfile,'magflag')))	
        print momflag,magflag,field
    else:		
        field=[]
        momflag=False
        magflag=True
    
#Set operation directory
    if(not inpath==''):
        os.chdir(inpath)
#The priority of infile is the first
    if (not infile==''):
        print 'read file %s?(y/n)\n' % (infile)
        confirm=raw_input()
        if (confirm=='y'):
            openfit(infile,outfile=outfile,field=field,magflag=magflag,momflag=momflag,BJD=BJD,z0=z0,fluxflag=fluxflag)		
            return
        else:
            print 'Exit...\n'
            exit(0)
        
#The priority of inlist is the second
    if(not inlist==''):
        print 'read all the FITS files from input list %s?(y/n)\n' % (inlist)
        confirm=raw_input()
        if (confirm=='y'):
            if(not inpath==''):
                os.chdir(inpath)
            flist=open(inlist,mode='r')
            line=flist.readline().split()
            while(not len(line)==0):
       	        if (not line[0].startswith('#')):
       	            iname=line[0]
                    if(os.path.exists(iname)):      
                        print iname
                        openfit(iname,outfile='',field=field,magflag=magflag,momflag=momflag,BJD=BJD,z0=z0,fluxflag=fluxflag)		
                line=flist.readline().split()
            flist.close()
            return
        else:
            print 'Exit...\n'
            exit(0)
        
#The priority of inpath is the last
    if(not inpath==''):
        os.chdir(inpath)
        print 'read all the FITS files from the input directory %s?(y/n)\n' % (inpath)
        s0=int(cfgp.readfits_parse(cfgfile,'s0'))
        se=int(cfgp.readfits_parse(cfgfile,'se'))
        qmin=int(cfgp.readfits_parse(cfgfile,'qmin'))
        qmax=int(cfgp.readfits_parse(cfgfile,'qmax'))
        inlist=cfgp.readfits_parse(cfgfile,'inlist')
        binpath=os.path.split(inlist)[0]+'/'
        ext=os.path.splitext(inlist)[1]
        print s0,se,qmin,qmax,binpath
        for sgn in range(s0,se+1):
            outpath=inpath+'Sgn-%d/'% (sgn)
            inlist=binpath+'Sgn-%d%s'% (sgn,ext)
            names=[]
            readcolumn(names,1,inlist)
            for q in range(qmin,qmax+1):
                qpath='Q%d_public/' % (q)
                for starid in names:
                    #if(1):
                    #starid=names[0]
                    x=qpath+'kplr%.9d' % (int(starid))+'-'+find_ext(q,sgn)+'_llc.fits'
        	        #print x
                    if (os.path.exists(x)):
                        openfit(x,outfile='',outpath=outpath,field=field,magflag=magflag,momflag=momflag,BJD=BJD,z0=z0,fluxflag=fluxflag)		
                    #else:
                    #    fout.write('%s %d %d\n' % (x,q,sgn))		
        #fout.close()
if __name__=='__main__':
	main()	
