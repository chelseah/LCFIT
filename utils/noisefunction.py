#!/usr/bin/env python
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
def cal_noise(Kepmag):
    c = 1.28*10**(0.4*(12-Kepmag)+7)
    N = 49032
    #sigma = 1.e6*(c+9.5e5*(14./Kepmag)**5.)**0.5/(c*N**(0.5))
    sigma = 1.e6*(c+9.5e5*(14./Kepmag)**5.)**0.5/c

    return sigma

def main():
    kp = 11+np.arange(100)/100.*3.
    sigma = cal_noise(kp)
    plt.plot(kp,sigma)
    plt.show()
    return

if __name__=='__main__':
    main()
