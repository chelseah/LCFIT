#!/usr/bin/env python
import math
from const import *
class lightcurve():
    def __init__(self):
        self.name=''
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
    def __str__(self):
        return '%13.6f %13.6f %6.4e %4.2e \n' % (self.P, self.epoch,self.dip,self.q)
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


