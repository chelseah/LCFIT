#!/usr/bin/env python
# openfile.py this file is used for data io control
from types import *
import pyfits
import math
import os
import string
import numpy as np

#===============================================================
#=====open different type of files==============================
#===============================================================
def open_fits(curve,infile):
	BJD=2454000
	if(infile.endswith('.fits')):	
                print infile
		hdulist=pyfits.open(infile)
                try:
                       	tbdata=hdulist[1].data
                  	BJDREFI=hdulist[1].header['BJDREFI']
                       	BJDREFF=hdulist[1].header['BJDREFF']
                        curve.info.BJD=BJDREFI+BJDREFF-BJD
                        time=tbdata.field('TIME')
                        flux=tbdata.field('SAP_FLUX')
                        mcent1=tbdata.field('MOM_CENTR1')
                        mcent2=tbdata.field('MOM_CENTR2')
			for i in range(len(time)):
				if not math.isnan(flux[i]):
					curve.data.time.append(time[i]+BJDREFI+BJDREFF-BJD)
					curve.data.flux.append(flux[i]/1e5)
					curve.data.mom1.append(mcent1[i])
					curve.data.mom2.append(mcent2[i])
			curve.info.id=hdulist[1].header['KEPLERID']
			if(os.path.splitext(infile)[0].endswith('llc')):
				curve.info.type='lc'
			else:
				curve.info.type='sc'
			curve.info.sn=hdulist[0].header['SEASON']
			curve.info.sgn=hdulist[0].header['SKYGROUP']		
			curve.info.mod=hdulist[0].header['MODULE']
			curve.info.chan=hdulist[0].header['CHANNEL']	
                        curve.info.ra=hdulist[1].header['RA_OBJ']
                        curve.info.dec=hdulist[1].header['DEC_OBJ']
			if type(hdulist[0].header['KMAG']) is FloatType:
				curve.info.Kmag=hdulist[0].header['KMAG']
			
			curve.info.Kpmag=hdulist[0].header['KEPMAG']
                except IOError:
                        print 'file read failed from %s' % (infile)
                finally:
                        hdulist.close()
		return True
	else:
		return False

