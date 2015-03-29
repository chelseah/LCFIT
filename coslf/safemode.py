#!/usr/bin/env python
import numpy as np
import scipy as sp
from lslft_recon import lspolyordern
from dataio import readcolumn,pendcol
import matplotlib 
from matplotlib import pyplot as plt
def safemode_cor(time,mag):
    order = 6 
    dmag,c = lspolyordern(time-time[0],mag,order)
    return [dmag,c]

def safemode_corr_q6(infile):
    time = []; readcolumn(time,1,infile); time = np.array(time)
    mag = []; readcolumn(mag,6,infile); mag = np.array(mag)
    safemode = [1400,1403.5]
    indexa = time>safemode[0]
    indexb = time<safemode[1]
    intran = [1401.00,1401.59]
    if not(intran==[]):
        indexc = time>intran[0]
        indexd = time<intran[1]
        index = indexa*indexb-indexc*indexd
        index1 = indexa*indexb 
    else: 
        index = indexa*indexb
        index1 = index
    dmag,c=safemode_cor(time[index],mag[index])
    order = 6 
    rflux = np.zeros(len(time[index1]))
    for i in range(order+1):
        rflux +=c[i]*(time[index1]-time[index1][0])**(order-i)
    dflux = mag[index1]-rflux
    print np.median(dflux),mag[index1][-1]
    mag[index1] = dflux+mag[index1][-1]
    pendcol(infile,mag,ext='.smpltf') 
    return

def main():
    inlist = ""
    inpath = ""
    inpath = "/Users/xuhuang/Documents/2014Spring/HEK/KOI351/PixelLC/safemode/"
    infile = "kplr011442793-2010265121752_llc.pflt"
    names = [infile]
    #names = []; readcolumn(names,1,inlsit,datformat='str')
    for i in xrange(len(names)):
        safemode_corr_q6(inpath+names[i])
    
    
    return 


def test():
    inpath = "/Users/xuhuang/Documents/2014Spring/HEK/KOI351/PixelLC/safemode/"
    infile = "kplr011442793-2010265121752_llc.pflt"
    time = []; readcolumn(time,1,inpath+infile); time = np.array(time)
    mag = []; readcolumn(mag,6,inpath+infile); mag = np.array(mag)
    safemode = [1400,1403.2]
    intran = [1401.00,1401.59]
    indexa = time>safemode[0]
    indexb = time<safemode[1]
    indexc = time>intran[0]
    indexd = time<intran[1]
    
    index = indexa*indexb-indexc*indexd
    dmag,c=safemode_cor(time[index],mag[index])
    order = 6 
    rflux = np.zeros(len(time[indexa*indexb]))
    for i in range(order+1):
        rflux +=c[i]*(time[indexa*indexb]-time[index][0])**(order-i)
    plt.plot(time,mag,'.',time[indexa*indexb],mag[indexa*indexb]-rflux+np.median(mag),'+')
    plt.xlim([1400,1403.2])
    plt.show()

if __name__=='__main__':
    #test()
    main()
