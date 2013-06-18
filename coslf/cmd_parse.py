#!/usr/bin/python
#lc_parse.py use to parse all the command lines for functions realted to run lc function. 
import optparse
def ltf_parse():
	p = optparse.OptionParser()
        p.add_option('--infile','-i',default='',help='the name of input file') 
        p.add_option('--outpath','-o',default='',help='the path for the output file(Optional)') 
        p.add_option('--inlist','-l',default='',help='the file contain a list of input files(Optional)') 
        p.add_option('--cfg','-c',default='/home/chelsea/src/newpipe/example.cfg',help='the filename of configuration file that contains parameter settings(Optional)') 
        p.add_option('--inpath','-p',default='',help='the path for gettting the input file(s)(Optional)') 
	p.add_option('--ext','-e',default='_llc.pflt',help='the ext of the files to trend filter and joining')
        p.add_option('--infileflag','-I',default=False,action='store_true',help='the flag for read input file from configration file')
        p.add_option('--outpathflag','-O',default=False,action='store_true',help='the flag for read output path from configration file')
        p.add_option('--inlistflag','-L',default=False,action='store_true',help='the flag for read input file list from configuration file')
        p.add_option('--inpathflag','-P',default=False,action='store_true',help='the flag for read input path from configuration file')
	p.add_option('--extflag','-E',default=False,action='store_true',help='the flag for read ext of the files to trend filter and joining')
        p.add_option('--uflag','-u',default=False,action='store_true',help='the flag for read other parameter settings from the configuration file') 
        p.add_option('--filt','-f',default=True,action='store_false',help='the flag for turn off the filt flat') 
        p.add_option('--noplot','-n',default=True,action='store_false',help='the flag for turn on the plotting flag during trend filtering, only recommended in debuging') 
        options,arguments = p.parse_args()
	return options

def runcfg():
        p = optparse.OptionParser()
	p.add_option('--eflag','-e',default=False,action='store_true',help='Creat example.cfg for all the parameter settings')
	p.add_option('--infile','-i',default='example.cfg',help='the file to write cfg parameters')
	p.add_option('--gapflag','-g',default=False,action='store_true',help='Creat gap_example.cfg for all the anomaly data for every quarter')
	options,arguments=p.parse_args()
	return options
