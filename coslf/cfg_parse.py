#!/usr/bin/python
#cfg_parse.py use to parse all the configure parameters from .cfg file for 
#functions related to run lc function. 

import ConfigParser
import cmd_parse as cmdp

def set_parse(infile):
    p = ConfigParser.RawConfigParser()
###################### Configuration parameter for Longtrendfilter ###########
    p.add_section('Section LTF')
    p.set('Section LTF','s0','1') #the start sgn number
    p.set('Section LTF','se','84') #the end sgn number
    p.set('Section LTF','inpath','')
    p.set('Section LTF','infile','')
    p.set('Section LTF','inlist','')
    p.set('Section LTF','outpath','')
    p.set('Section LTF','outfile','')
    p.set('Section LTF','coljd',1)
    p.set('Section LTF','colmag',2)
    p.set('Section LTF','tmin',1.0)
    p.set('Section LTF','period','')
    p.set('Section LTF','epoch','')
    p.set('Section LTF','tdur','')

#################### The General Setting part #############################
	
    p.add_section('This is a configure file for xxx package')

    with open (infile,'wb') as configfile:
        p.write(configfile)	
        return 

def ltf_parse(infile,name):

    p = ConfigParser.RawConfigParser()
    try:
        p.read(infile)
        var=p.get('Section LTF',name)

    except IOError:
        print '\n'
    except ConfigParser.MissingSectionHeaderError:
        print 'Error: the section LTF is missing, excute set_parse to see example.cfg\n'
    except ConfigParser.NoOptionError:
        print 'Error: the option %s is missing, excute set_parse to see example.cfg\n' % name
    return var 

if __name__=='__main__':
    options=cmdp.runcfg()
    if(options.eflag):	
        infile=options.infile
        set_parse(infile)

