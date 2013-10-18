#!/usr/bin/env python
import numpy as np
import scipy as sp
import math
import os
import matplotlib 
from matplotlib import plt
from dataio import readcolumn
from lcgen import gentran
def main()
    inlist = 'B13candidate.txt'
    colid = 2  
    colp = 6
    cole = 9
    colTdur =  15
    colrpstar = 12 
    coljd = 1
    colmag = 4
    names = []; readcolumn(names,colid,inlist,datformat='str')
    period = []; readcolumn(period,colp,inlist); period=np.array(period)
    epoch = []; readcolumn(epoch,cole,inlist);epoch = np.array(epoch)
    rpstar = []; readcolumn(rpstar,colrpstar,inlist); rpstar = np.array(rpstar)
    Tdur = []; readcolumn(Tdur,colTdur,inlist); Tdur = np.array(Tdur)

    fout = open('KOIootvinfo','w')
    inpath = '/home/chelsea/mntproj/fullcand/TFA/'
    for i in len(names):
        infile = inpath+'kplr%.9d.tfa' % (int(names[i]))
        if(not os.path.exsits(infile)):
            print infile
            continue
        oot1 = gentran(time,period[i],epoch[i]-Tdur[i],Tdur[i]/period[i]/24.)
        oot2 = gentran(time,period[i],epoch[i]+Tdur[i],Tdur[i]/period[i]/24.)
        ootv1 = np.std(mag[oot1])
        ootv2 = np.std(mag[oot2])
        ootv = ((ootv1**2.+ootv2**2.)/2.)**0.5
        fout.write("%s %f %f %fi %f\n" % (names[i],period[i],epoch[i],rpstrar[i],ootv))        
    fout.close()
    return

if __name__=='__main__':
    main()
