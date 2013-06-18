import os
import sys
from numpy import *
import matplotlib.pyplot as plt
#from astropy_transit import *
from scipy import optimize
from scipy import interpolate
import emcee
import random
import functions
from ctypes import *
from numpy import *
import string

import warnings
warnings.filterwarnings('ignore')

### Compile shared transit library
os.system("rm jktebop_lib.so")
os.system("gfortran jktebop_lib.f -ffree-form -fpic -shared -o jktebop_lib.so")
jktebop = cdll.LoadLibrary("./jktebop_lib.so")

def transitmodel(V,phase,la,lb):

    ldtype = array([4,4])
    dtype = c_long(1)
    la = c_double(la)
    lb = c_double(lb)
    
    mag = []
    for i in phase:

        f = c_double(0.0)
        i = c_double(i)
        
        jktebop.getmodel_(V.ctypes.data_as(POINTER(c_double)),ldtype.ctypes.data_as(POINTER(c_long)),byref(i),byref(dtype),byref(la),byref(lb),byref(f))

        mag.append(ctypeslib.as_array(f)[0])

    mag = array(mag)

    flux = 10**(mag/-2.5)

    return flux

### Syntax for transit model
### transitmodel(V.ctypes.data_as(POINTER(c_double)),ldtype.ctypes.data_as(POINTER(c_long)),byref(phase),byref(dtype),byref(la),byref(lb),byref(flux))
### dtype = 1, ldtype = [4,4]

### Global constants
au = 1.496*10**11
msun = 1.988435*10**30
rsun = 6.955*10**8
mjup = 1.8988*10**27 
rjup = 6.9173*10**7
day = 60.*60.*24.
gconst = 6.67*10**(-11)

t0_global = floor(eval(functions.read_config_file("T0")))

tested_params = []
chisq_log = []

def lc_chisq(initial_params,free_param_names,fixed_param_names,fixed_param_values,lc,plot_pdf,output_data):
    #print initial_params
    
    global period_i,t0_i,rsum_i,rratio_i,i_0_i,ld1_i,ld2_i
    ### Give dummy values to avoid error
    [period_i,t0_i,rsum_i,rratio_i,i_0_i,ld1_i,ld2_i] = [1,1,1,1,1,1,1]

    #offsets = []
    ld1_coeff = []
    ld2_coeff = []
    ### Set parameter names
    for i in range(len(free_param_names)):
        globals()[free_param_names[i]+"_i"] = initial_params[i]

        if free_param_names[i] == "lc_ld1":
            ld1_coeff.append(initial_params[i])
        if free_param_names[i] == "lc_ld2":
            ld2_coeff.append(initial_params[i])

    for i in range(len(fixed_param_names)):
        globals()[fixed_param_names[i]+"_i"] = fixed_param_values[i]

        if fixed_param_names[i] == "lc_ld1":
            ld1_coeff.append(fixed_param_values[i])
        if fixed_param_names[i] == "lc_ld2":
            ld2_coeff.append(fixed_param_values[i])

    t0_i = t0_i + t0_global

    chisq = 0
    npoints = 0

    ### Input for transit model 
    # ! V(1) = surface brightness ratio    V(15) = third light
    # ! V(2) = sum of fractional radii     V(16) = phase correction
    # ! V(3) = ratio of stellar radii      V(17) = light scaling factor
    # ! V(4) = linear LD for star A        V(18) = integration ring size (deg)
    # ! V(5) = linear LD for star B        V(19) = orbital period (days)
    # ! V(6) = orbital inclination         V(20) = ephemeris timebase (days)
    # ! V(7) = e cos(omega) OR ecentricity V(21) = nonlinear LD for star A
    # ! V(8) = e sin(omega) OR omega       V(22) = nonlinear LD for star B
    # ! V(9) = gravity darkening 1         
    # ! V(10)= gravity darkening 2         
    # ! V(11) = primary reflected light    
    # ! V(12) = secondary reflected light  
    # ! V(13) = stellar mass ratio        
    # ! V(14) = tidal lead/lag angle (deg)

    hjd_i,flux_i,fluxerr_i = lc[:,0],lc[:,1],lc[:,2]

    model_input = array([0,rsum_i,rratio_i,ld1_coeff[0],0,i_0_i,ecosw_i,esinw_i,0,0,0,0,0,0,0,0,0,1,period_i,t0_i,ld2_coeff[0],0])

    #print len(model_input),"model_input_length"
    model = transitmodel(model_input,hjd_i,1.,0.)

    ### Apply offset
    x0 = [median(flux_i)]
    def minfunc(x0):
        flux_ii = flux_i + x0[0]
        chisq_i =  sum(((flux_ii-model)/fluxerr_i)**2)
        return chisq_i

    x0 = optimize.fmin(minfunc,x0,disp=0)
    flux_i = flux_i + x0[0]

    diff = flux_i - model

    chisq = chisq + sum((diff/fluxerr_i)**2)
    npoints = npoints+float(len(diff))

    model_i = model

    if plot_pdf:

        Mt=2*pi*((hjd_i-t0_i)/period_i - floor((hjd_i-t0_i)/period_i))

        plt.errorbar(Mt/(2*pi),flux_i,fluxerr_i,marker="o",markersize=5,linestyle="None",color="k")
        plt.errorbar(Mt/(2*pi)+1,flux_i,fluxerr_i,marker="o",markersize=5,linestyle="None",color="k")

        hjd_m = arange(t0_i,t0_i+period_i,0.001)

        Mt=2*pi*((hjd_m-t0_i)/period_i - floor((hjd_m-t0_i)/period_i))

        model = transitmodel(model_input,hjd_m,1.0,0.0)

        plt.scatter(Mt/(2*pi),model,s=5,facecolor="r",edgecolor="r")
        plt.scatter(Mt/(2*pi)+1,model,s=5,facecolor="r",edgecolor="r")

        plt.xlim(0.99,1.01)

        plt.show()

    if output_data:
        Mt=2*pi*((hjd_i-t0_i)/period_i - floor((hjd_i-t0_i)/period_i))
        Mt_data = Mt/(2*pi)
        return Mt_data,flux_i,fluxerr_i,model_i

    else:
        return chisq


def gaussian(x,x0,c):
    return exp(-(x-x0)**2 / (2*c**2))

def box(x,x0,c):
    if abs(x-x0) > c or x < 0:
        return NaN
    else:
        return 1.0

def calc_master_chisq(initial_params,default_params,free_param_names,fixed_param_names,fixed_param_values,prior_params,prior_mean,prior_std,prior_func,lc,plot_pdf):
    initial_params = array(initial_params)*array(default_params)+array(default_params)
    
    chisq = lc_chisq(initial_params,free_param_names,fixed_param_names,fixed_param_values,lc,plot_pdf,False)

    for i in range(len(prior_params)):
        for j in range(len(free_param_names)):
            if prior_params[i] == free_param_names[j]:
                if prior_func[i] == "b":
                    chisq = chisq * box(initial_params[j],prior_mean[i],prior_std[i])

    #print chisq
    return chisq


def calc_probability(initial_params,default_params,free_param_names,fixed_param_names,fixed_param_values,prior_params,prior_mean,prior_std,prior_func,lc,chisq_base,plot_pdf):

    initial_params = array(initial_params)*array(default_params)+array(default_params)

    chisq = lc_chisq(initial_params,free_param_names,fixed_param_names,fixed_param_values,lc,plot_pdf,False)

    global stellar_params,tested_params,chisq_log

    prob = (chisq_base-chisq)/2

    ### Calculate stellar parameters
    global period_i,t0_i,rsum_i,rratio_i,i_0_i,ld1_i,ld2_i
    ### Give dummy values to avoid error
    [period_i,t0_i,rsum_i,rratio_i,i_0_i,ld1_i,ld2_i] = [1,1,1,1,1,1,1]

    ### Set parameter names
    for i in range(len(free_param_names)):
        globals()[free_param_names[i]+"_i"] = initial_params[i]

    for i in range(len(fixed_param_names)):
        globals()[fixed_param_names[i]+"_i"] = fixed_param_values[i]

    i_0_i = i_0_i * pi / 180.
    
    e_0_i = sqrt(ecosw_i**2+esinw_i**2)
    w_0_i = arccos(ecosw_i/e_0_i)

    a_0_i = (rsum_i)/(1+rratio_i)

    for i in range(len(free_param_names)):
        if prior_func[i] == "g":
            prob = prob * gaussian(initial_params[i],prior_mean[i],prior_std[i])
        if prior_func[i] == "b":
            prob = prob * box(initial_params[i],prior_mean[i],prior_std[i])
            #print prob,initial_params[i],prior_mean[i],prior_std[i]
    
    if functions.isnan(prob):

        tested_params = [list(initial_params)]
        chisq_log = [[chisq]]


        f = open("test","w")
        functions.write_table(tested_params,f)
        f.close()
        os.system("cat test >> mcmc_tested_params")

        f = open("test","w")
        functions.write_table(chisq_log,f)
        f.close()
        os.system("cat test >> mcmc_chisq_log")

    print prob
    return prob

def mcmc_loop(initial_params,default_params,free_param_names,fixed_param_names,fixed_param_values,prior_params,prior_mean,prior_std,prior_func,lc,plot_pdf):

    import random

    global chisq_log,tested_params
    
    ### Run an intial round to get baseline chisq

    initial_params_temp = initial_params * default_params + default_params

    chisq_base = lc_chisq(initial_params_temp,free_param_names,fixed_param_names,fixed_param_values,lc,plot_pdf,False)

    print "Running MCMC to sample the parameter space"

    nwalkers = int(eval(functions.read_config_file("WALKERS")))
    nmcmc = int(eval(functions.read_config_file("MCMC")))
    nburn = int(eval(functions.read_config_file("NBURN")))
    nthreads = 1

    param_tolerance_names = ["ecosw","esinw","period","t0"]
    param_tolerances = [0.1,0.1,0.00001,0.000001]
    default_tolerance = 0.001

    chisq_log,stellar_params,tested_params= [],[],[]

    ndim = len(initial_params)
    #nwalkers = 20
    p0 = []
    for i in range(nwalkers):
        pi = []
        for j in range(len(initial_params)):
            tolerance = default_tolerance
            if j < len(free_param_names):
                for k in range(len(param_tolerance_names)):
                    if param_tolerance_names[k] == free_param_names[j]:
                        tolerance = param_tolerances[k]
                        break
            pi.append(random.gauss(initial_params[j],tolerance))
        p0.append(pi)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, calc_probability, args=[default_params,free_param_names,fixed_param_names,fixed_param_values,prior_params,prior_mean,prior_std,prior_func,lc,chisq_base,plot_pdf],threads=nthreads)

    pos, prob, state = sampler.run_mcmc(p0, nburn)
    sampler.reset()

    master_pos = []
    master_prob = []
    best_prob = -9999999
    best_pos = p0[0]

    for result in sampler.sample(pos, iterations=nmcmc, storechain=False):
        position,probability = result[0],result[1]

        for i in range(len(position)):
            if functions.isnan(probability[i]):
                if probability[i] > best_prob:
                    best_pos = list(position[i])
                    best_prob = probability[i]

                master_prob.append(probability[i])
                master_pos.append(list(position[i]))


    print "iteration finished"
        
    print best_pos

    return array(best_pos)


def inflate_errors(x0,free_param_names,free_param_vals,fixed_param_names,fixed_param_vals,lc):

    ### Calculate stellar parameters
    global period_i,t0_i,rsum_i,rratio_i,i_0_i,ld1_i,ld2_i
    ### Give dummy values to avoid error
    [period_i,t0_i,rsum_i,rratio_i,i_0_i,ld1_i,ld2_i] = [1,1,1,1,1,1,1]

    temp_params = x0*free_param_vals+free_param_vals
    ### Set parameter names

    ld1_coeff = []
    ld2_coeff = []

    ### Set parameter names
    for i in range(len(free_param_names)):
        globals()[free_param_names[i]+"_i"] = temp_params[i]
        if free_param_names[i] == "lc_ld1":
            ld1_coeff.append(temp_params[i])
        if free_param_names[i] == "lc_ld2":
            ld2_coeff.append(temp_params[i])

    for i in range(len(fixed_param_names)):
        globals()[fixed_param_names[i]+"_i"] = fixed_param_vals[i]
        if fixed_param_names[i] == "lc_ld1":
            ld1_coeff.append(fixed_param_vals[i])
        if fixed_param_names[i] == "lc_ld2":
            ld2_coeff.append(fixed_param_vals[i])


    t0_global = floor(eval(functions.read_config_file("T0")))
    t0_i = t0_i + t0_global

    #### Calculate error inflation
    hjd_i,flux_i,fluxerr_i = lc[:,0],lc[:,1],lc[:,2]

    model_input = array([0,rsum_i,rratio_i,ld1_coeff[0],0,i_0_i,ecosw_i,esinw_i,0,0,0,0,0,0,0,0,0,1,period_i,t0_i,ld2_coeff[0],0])
    model = transitmodel(model_input,hjd_i,1,0)

    ### Apply offset
    z0 = [median(flux_i)]
    def minfunc(z0):
        flux_ii = flux_i + z0[0]
        chisq_i =  sum(((flux_ii-model)/fluxerr_i)**2)
        return chisq_i

    z0 = optimize.fmin(minfunc,z0,disp=0)
    flux_i = flux_i + z0[0]

    df = float(len(flux_i))

    def chisq(y,dy,f,df):
        return sum(((f-y)/dy)**2)/df

    def minfunc(factor,y,dy,f,df):
        x2 = chisq(y,dy*factor,f,df)
        x2 = abs(x2-1)
        return x2

    s0 = 1.
    s0 = optimize.fmin(minfunc,s0,args=(flux_i,fluxerr_i,model,df))
    s0 = s0[0]
    print "inflating errors by factor of ",s0
    fluxerr_i = fluxerr_i*s0

    lc[:,2] = fluxerr_i
    return lc
