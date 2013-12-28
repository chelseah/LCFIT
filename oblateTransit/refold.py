# !/usr/bin/env python
import os
import numpy as np
import scipy as sp
import math
from dataio import *

def phasefold(time,epoch,period):
	ftime=sp.zeros(len(time))
	ftime=(time-epoch-0.5*period)/period-((time-epoch-0.5*period)/period).astype(int)
	ind=ftime<0
	ftime[ind]+=1
	return ftime

def hourfold(time,epoch,period):
	ftime=sp.zeros(len(time))
	ftime=(time-epoch-0.5*period)/period-((time-epoch-0.5*period)/period).astype(int)
	ind=ftime<0
	ftime[ind]+=1
	ftime-=0.5
	ftime*=period*24
	return ftime

def rawfold(time,epoch,period):
	ftime=sp.zeros(len(time))
	ftime=time-period*((time-epoch-0.5*period)/period).astype(int)
	return ftime

