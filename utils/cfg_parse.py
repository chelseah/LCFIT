#!/usr/bin/python
#cfg_parse.py use to parse all the configure parameters from .cfg file for 
#functions related to run lc function. 

import ConfigParser
import cmd_parse as cmdp

def set_parse(infile):
        p = ConfigParser.RawConfigParser()
	
######################## Configuration parameter for Readfits##################
	p.add_section('Section FITS')
	p.set('Section FITS','infile','')
	p.set('Section FITS','inlist','')
	p.set('Section FITS','inpath','')
	p.set('Section FITS','outfile','')
	p.set('Section FITS','outpath','')
	p.set('Section FITS','field',[])
	p.set('Section FITS','momflag',False)
	p.set('Section FITS','magflag',True)
	
#################### The General Setting part #############################
	p.add_section('General setting')
	p.set('General setting','BJD',2454000)
	p.set('General setting','z0',0)
	p.set('General setting','qmax',6)
	p.set('General setting','qmin',1)
	
	p.add_section('This is a configure file for xxx package')

	with open (infile,'wb') as configfile:
		p.write(configfile)	
        return 

def gen_parse(infile,name):	
        p = ConfigParser.RawConfigParser()
	try:
		p.read(infile)
		var=p.get('General setting',name)
	except IOError:
		print '\n'
	except ConfigParser.MissingSectionHeaderError:
		print 'Error: the format of cfg file is not correct, excute set_parse to see example.cfg\n'
	except ConfigParser.NoOptionError:
		print 'Error: the option of cfg file is not correct, excute set_parse to see example.cfg\n'
		
	return var

def readfits_parse(infile,name):

        p = ConfigParser.RawConfigParser()
	try:
		p.read(infile)
		var=p.get('Section FITS',name)
	except IOError:
		print '\n'
	except ConfigParser.MissingSectionHeaderError:
		print 'Error: the format of cfg file is not correct, excute set_parse to see example.cfg\n'
	except ConfigParser.NoOptionError:
		print 'Error: the option of cfg file is not correct, excute set_parse to see example.cfg\n'

	return var

def set_gap():
#set up the cfg file for gapinfo	
        p = ConfigParser.RawConfigParser()
####################### Configuration parameter for SC##################
	p.add_section('Section SC Gaps')
	p.set('Section SC Gaps','q1','')
	p.set('Section SC Gaps','q2','')
	p.set('Section SC Gaps','q3','')
	p.set('Section SC Gaps','q4','')
	p.set('Section SC Gaps','q5','')
	p.set('Section SC Gaps','q6','')

######################## Configuration parameter for LC##################
	p.add_section('Section LC Gaps')
	p.set('Section LC Gaps','q1','985.272,985.474')
	p.set('Section LC Gaps','q2','1002.05,1005.73;1014.5186,1019.7214;1033.1597,1033.3657;1056.4815,1056.8494;1059.7723,1059.77236;1063.30845,1066.4153;1073.0174,1073.5284;1079.1901,1079.19018;1080.8049,1080.82537;1086.60,1086.70;1088.4085,1089.3283;1058.08,1058.26')
	p.set('Section LC Gaps','q3','1092.7222,1092.9;1096.6455,1096.8;1099.9148,1100.0374;1113.0536,1113.8301;1123.0661,1123.9039;1129.0736,1129.2;1132.8538,1133.0;1139.02,1139.05;1153.9617,1156.9400;1179.3607,1179.5')
	p.set('Section LC Gaps','q4','1186.2774,1186.2979;1206.28,1206.8;1213.9445,1213.9650;1217.3,1218.73;1220.1768,1220.1972;1229.3515,1236.0;1266.17,1266.20')
	p.set('Section LC Gaps','q5','1309,1310.8;1319.55559,1319.8;1321.59895,1321.7;1324.97174,1326.1;1336.93,1338.02')
	p.set('Section LC Gaps','q6','1382.45,1382.70;1395.0,1395.05;1400.0,1401.9;1423.99172,1424.1;1432.0,1434.4')
	p.set('Section LC Gaps','q7','1462.6725,1463.76;1487.4993,1487.4993;1492.7916,1494.42;1522.7473,1525.34;1542.7314,1542.7314;1543.4261;1543.9165;1546.6342,1546.6751')
	p.set('Section LC Gaps','q8','1567.8647,1568.83;1572.9935,1572.9935;1582.6995,1582.6995;1593.2228,1593.2228;1593.5702,1598.21;1613.3908,1613.3908;1619.8069,1619.8069;1621.5847,1621.5847')
	p.set('Section LC Gaps','q9','1677.4296,1678.1040;1706.6293,1707.2423;1719.3390,1719.3390;1720.7285;1720.7285')
	p.set('Section LC Gaps','q10','1739.3434,1739.3534;1762.3466,1762.3671;1769.4560,1770.2731;1799.2008,1799.2108;1801.7340,1802.5512;1828.3532,1828.3532')
	p.set('Section LC Gaps','q11','1833.7058,1833.7058;1836.6476,1836.6476;1836.7906,1836.8110;1840.7947,1840.7948;1840.9786,1840.9886;1842.4904,1842.5004;1853.8286,1853.8386;1859.6918,1859.7532;1862.6541,1862.6950;1864.7787,1865.5347;1870.7849,1870.8849;1873.7676,1873.7677;1895.7292,1896.6076;1902.4912,1905.0245;1915.4434,1915.5434;1926.1279,1926.2280;1926.3527,1926.4527;1928.8654,1928.9472;1929.1106,1929.7644')
	p.set('Section LC Gaps','q12','1931.9097,1932.;1935.4642,1935.5643;1936.7300,1936.8308;1942.8185,1942.8389;1945.6989,1945.81;1949.1922,1951.3168;1953.7886,1954.7080;1956.0358,1956.1358;1958.4055,1959.1305;1962.5729,1962.6729;1966.5360,1966.6360;1967.6392,1967.7392;1986.4947,1987.3936;1993.0318,1997.5262;1998.6088704,1998.7702;2004.4310,2004.5310;2009.9263,2010.02;2011.3358,2011.4358')
	p.add_section('This is a configure file for anomaly point information')
	with open ('gap_example.cfg','wb') as configfile:
		p.write(configfile)	
        return 

def LC_parse(infile,name):

    p = ConfigParser.RawConfigParser()
    try:
        p.read(infile)
        rawvar=p.get('Section LC Gaps',name)
        intvar=rawvar.split(';')
        var=[]
        for x in intvar:
        	var.append(map(float,x.split(',')))
			
    except IOError:
        print '\n'
    except ValueError:
        print x.split(',')
        if x=='':
            var = []
        else:
            raise
    except ConfigParser.MissingSectionHeaderError:
        print 'Error: the format of cfg file is not correct, excute set_parse to see example.cfg\n'
    except ConfigParser.NoOptionError:
        print 'Error: the option of cfg file is not correct, excute set_parse to see example.cfg\n'

    return var

if __name__=='__main__':
	options=cmdp.runcfg()
	if(options.eflag):	
		infile=options.infile
		set_parse(infile)
		
	if(options.gapflag):
		set_gap()

