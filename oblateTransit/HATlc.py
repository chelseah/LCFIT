#!/usr/bin/env python
import math
import numpy as np
from const import *
import cfg_parse as cfgp
class lightcurve():
    def __init__(self):
        self.name=''
        self.cadence='long'
        self.coljd=-1
        self.colmag=-1
        self.outext=''
        self.ap = [True,True,True]

class tran():
    def __init__(self):
        self.P= -1
        self.q= -1
        self.qg = 0.25
        self.dip= -1
        self.epoch= -1
        self.f = 0.0
        self.inc = 90.
        self.alpha = 0.0
        self.u1 = -1
        self.u2 = -1
        self.q1 = -1
        self.q2 = -1
        self.stdmag = -1
        self.rmean = -1
        self.rsum = -1
        self.b = -1
    def __str__(self):
        return '%13.6f %13.6f %6.4e %4.2e \n' % (self.P, self.epoch,self.dip,self.q)
    def readpara(self,cfgfile):
        self.P = float(cfgp.tran_parse(cfgfile,'period'))
        self.u1 = cfgp.tran_parse(cfgfile,'u1')
        self.u2 = cfgp.tran_parse(cfgfile,'u2')
        self.q1 = cfgp.tran_parse(cfgfile,'q1')
        self.q2 = cfgp.tran_parse(cfgfile,'q2')

        if not self.u1 == '':
            self.u1 = float(self.u1)
        if not self.u2 == '':
            self.u2 = float(self.u2)
        if not self.q1 == '':
            self.q1 = float(self.q1)
        if not self.q2 == '':
            self.q2 = float(self.q2)
        self.rmean = float(cfgp.tran_parse(cfgfile,'rpstar'))
        self.rsum = float(cfgp.tran_parse(cfgfile,'rsum'))
        self.inc = float(cfgp.tran_parse(cfgfile,'inc'))
        self.alpha = float(cfgp.tran_parse(cfgfile,'alpha'))
        self.f = float(cfgp.tran_parse(cfgfile,'planetf'))
        self.dip = float(cfgp.tran_parse(cfgfile,'rpstar'))**2.
        self.epoch =float(cfgp.tran_parse(cfgfile,'epoch'))
        self.sma = (1+self.rmean)/self.rsum

        Tdur = cfgp.tran_parse(cfgfile,'Tdur')
        if not Tdur == '':
            self.q = float(Tdur)/24./self.P 
        else: 
            self.q = 1./self.sma/np.pi
            #self.q = self.calq() 
            #print 'Warning: the default q is computed with solar radii and inc=90.'
        qg = cfgp.tran_parse(cfgfile,'qgress')
        if not qg == '': 
            self.qg = float(qg)
        f =cfgp.tran_parse(cfgfile,'planetf')
        if not f == '':
            self.f =float(f)
        inc = cfgp.tran_parse(cfgfile,'inc')
        if not inc=='':
            self.inc = float(inc)
            self.b = np.cos(math.pi/180.*self.inc)*self.sma
        alpha = cfgp.tran_parse(cfgfile,'alpha')
        if not alpha=='':
            self.alpha = float(alpha)

        stdmag = cfgp.tran_parse(cfgfile,'stdmag')
        if not stdmag=='':
            self.stdmag = float(stdmag)
    def calLD(self):
        self.u1 = 2.*math.sqrt(self.q1)*self.q2
        self.u2 = math.sqrt(self.q1)*(1-2*self.q2)
        return
    def calq(self,rstar=1,logg=4.5):
        if(self.dip<=0):
            raise ValueError('A real planet to calculate qexp should have positive depth')
            return
        else:
            rpstar=math.sqrt(self.dip)
            sqrtgm=math.sqrt(10**(logg)*(rstar*rsun)**2.)
            periods=self.P*day
            sqrtsemia=(periods*sqrtgm/(2*math.pi))**(1./3.)
            vcir=2*math.pi*sqrtsemia**2./periods
            dur=2*(rpstar+1)*rstar*rsun/vcir
            qexp=dur/periods
            return qexp
    def calZ(self,phase):
    
        z=self.sma*np.sqrt(np.sin(phase*2*math.pi)*np.sin(phase*2*math.pi)+np.cos(phase*2*math.pi)*np.cos(phase*2*math.pi)*np.cos(self.inc)*np.cos(self.inc))

        return z

