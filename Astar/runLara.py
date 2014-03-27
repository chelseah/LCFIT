#!/usr/bin/env python
import quadpolar
import quad2d
import LaraModelfilter
from Larafilter import Larafilter as Lara
import math
import numpy as np
import scipy as sp
import pyfits as pyf
from dataio import readcolumn
def getspecfile(filtername):
    specfile = 'spec_intens_fullgrid.' + filtername + '.fits'
    return specfile
def parseline(line):
    return [float(line[1]),float(line[2]),float(line[3]),float(line[6]),float(line[7]),float(line[8]),float(line[9])]

def main():
    table = 'rotating_table_z0.014.dat'
    filtertable = 'spec_table.txt'
    #Filters = ['2MASS_J','2MASS_H','2MASS_Ks','Johns_U','Johns_B','Johns_V','Johns_R','Johns_I','SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z','TYCHO_B','TYCHO_V']
    Filters = []; readcolumn(Filters,1,filtertable,datformat='str')
    ZP = []; readcolumn(ZP,7,filtertable)
    phigrid = np.array([0.]) #mu=cos(theta)
    fin = open(table,mode='r')
    for line in fin.readlines():
        if line.startswith('#'):
            continue
        #line = '1.000     1.000     0.157     3.762    1.0000  8.4954e+09    1.0000  0.0000e+00  8.3353e+05    0.0000\n'
        line = '1.000     1.000     0.005     3.762    0.5581  8.4954e+09    1.0000  0.0000e+00  6.9902e+05    0.0000\n'
        Rsun = 6.955e10
        Msun = 2.e33
        Mass, logL, logTe,Obl,Omega,Req,OOc = parseline(line.split())
        Req/=(Rsun/1.e5)
        fratio = 1-Obl
        #print line.rstrip(),
        #G_MRsun = 4*math.pi**2*(1.5e3/6.955)**3.
        ggraveq0 = 6.67e-8*Mass*Msun/(Req*Rsun)**2.
        print Mass, Req, 10.**logTe ,np.log10(ggraveq0)
        groteq0 = Omega**2.*(Req*Rsun)
        #loggsolar = np.log10(6.67e-8*2.e33/(6.955e10)**2.)
        #print np.log10(ggraveq0), logTe, Mass, Req
        for filtername in Filters:
            specfile = getspecfile(filtername)
            hdu = pyf.open(specfile)[0]
            # axes are log(g), log(T), mu (fits header should be right)
            nx, ny, nz = hdu.data.shape
            #print nx,ny,nz,len(hdu.data.flatten())
            logg_val = hdu.header['CRVAL1']; logg_delt=hdu.header['CDELT1']
            logT_val = hdu.header['CRVAL2']-logTe; logT_delt=hdu.header['CDELT2']
            mu_val = hdu.header['CRVAL3']; mu_delt=hdu.header['CDELT3']
            # I might want to scale the units for logg and T first before put in Lara 
            #print "in python", hdu.data[0,0,0], hdu.data[3,4,5], hdu.data[11,12,13], hdu.data[14,1143,0], nx, ny, nz
            
            #print type(hdu.data)
            phi = phigrid[0]
            #ggraveq0 = 1.
            #groteq0 = 1.
            gdmodel = Lara(hdu.data.flatten().astype(np.float64),nx,ny,nz,logg_val,logg_delt,logT_val,logT_delt,mu_val,mu_delt,fratio,phi,Req,ggraveq0,groteq0)
            
            for i in xrange(len(phigrid)):
                if i>0:
                    #F0=2
                    gdmodel.Resetphi(phigrid[i])
                F0=gdmodel.Cal_F0()
                print F0,
                lum=gdmodel.Cal_Lum()
                print lum
                #print F0/(10.*3.0857e18/6.955e10)**2.,
            #break
        print ''
        break
    return
if __name__=='__main__':
    main()
