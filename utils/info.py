# !/usr/bin/env python
import numpy as np
import scipy as sp
import cfg_parse as cfgp

def __init__():
	BJDmin=cfgp.gen_parse()

def find_Q(time):
	BJDmin=964.081
	BJDmax=2014.5227
	BJD=2454000
	if (time>BJD):
		time-=BJD
	
	if(964.081<time<997.481):
		Qn=1
	elif(time<964.081<time):
		Qn=0
	elif(1002.0175<time<1090.9649):
		Qn=2
	elif(1092.7222<time<1181.9966):
		Qn=3
	elif(1184.8778<time<1274.7038):
		Qn=4
	elif(1275.9912<time<1367.5096):
		Qn=5
	elif(1371.9473<time<1461.7936):
		Qn=6
	elif(1462.6725<time<1552.0491):
		Qn=7
	elif(1567.8647<time<1634.8460):
		Qn=8
	elif(1641.016958<time<1738.423953):
		Qn=9
	elif(1739.343433<time<1832.765870):
		Qn=10
	elif(1833.7058<time<1930.8267):
		Qn=11
	elif(1931.9097<time<2014.5227):
		Qn=12
	elif(2015.2379<time<2105.5544):
		Qn=13
	elif(2107.1598<time<2204.3214):
		Qn=14
	elif(2206.5083<time<2304.1363):
		Qn=15
	elif(2304.2<time<2425.0):
		Qn=16
	elif(2392.2<time):
		Qn=17
	else:
		print 'The time input exceeds the time range of public data\n'
		print 'tmin(BJD)-%d=%d\n' % (BJD,BJDmin) 	
		print 'tmax(BJD)-%d=%d\n' % (BJD,BJDmax)
		print 'accept input in both BJD itself or BJD-%d\n' % (BJD)
		exit(1)
	return Qn


def find_Qf(infile):
	kid=infile.split('_')[0]
	try:
		qid=int(kid.split('-')[1])
		if (qid==2009131105131):
			Qn=0
		elif (qid==2009166043257):
			Qn=1
		elif (qid==2009259160929):
			Qn=2
		elif (qid==2009350155506):
			Qn=3
		elif (qid==2010078095331):
			Qn=4
		#	name4='/home/chelsea/mntproj/Q4_public/kplr'+x+'-'+'2010009091648_llc.fits'
		elif (qid==2010009091648):
			Qn=4
		elif (qid==2010174085026):
			Qn=5
		elif (qid==2010265121752):
			Qn=6
		elif (qid==2010355172524):
			Qn=7
		elif(qid==2011073133259):
			Qn=8
		elif(qid==2011177032512):
			Qn=9
		elif(qid==2011271113734):
			Qn=10
		elif(qid==2012004120508):
			Qn=11
		elif(qid==2012088054726):
			Qn=12
		elif(qid==2012179063303):
			Qn=13
		elif(qid==2012277125453):
			Qn=14
		elif(qid==2013011073258):
			Qn=15
		elif(qid==2013098041711):
			Qn=16
		elif(qid==2013131215648):
			Qn=17
		else:
			Qn=-1
			print 'FormatError: could not find the quarter number according to the filename\n'	
			print 'The standard file name: kplrxxxxxxxxx-xxxxxxxxxxx_xxx.xxx'
		return Qn
	except ValueError:
		print 'FormatError: could not find the integer part of filename (JD of the quarter) correspoinding to the quarter number\n'	
		print 'The standard file name: kplr%d-%d_xxx.xxx % (starid,quarterJD) \n'	
		print 'Use find_ext(Qn) to construct your filename as \n'
		print '\'kplr\'+starid+\'-\'+find_ext(Qn)+\'_\'+ext\n'
		exit(1)
	except IndexError:
		print 'FormatError: could not split filename into starid and quarternumber\n'
		print 'The standard file name: kplr%d-%d_xxx.xxx % (starid,quarterJD)\n'	
		print 'Use find_ext(Qn) to construct your filename as \n'
		print '\'kplr\'+starid+\'-\'+find_ext(Qn)+\'_\'+ext\n'
		

def find_ext(Qn,sgn=-1):
    qmax=13
    if (Qn==0):
        kext='2009131105131'
    elif (Qn==1):
        kext='2009166043257'
    elif(Qn==2):
        kext='2009259160929'
    elif(Qn==3):
        kext='2009350155506'
    elif(Qn==4):
        if((4<sgn) and (sgn<9)):
            kext='2010009091648'
        else:
            kext='2010078095331'
    elif(Qn==5):
        kext='2010174085026'
    elif(Qn==6):
        kext='2010265121752'
    elif(Qn==7):
        kext='2010355172524'
    elif(Qn==8):
        kext='2011073133259'
    elif(Qn==9):
        kext='2011177032512'
    elif(Qn==10):
        kext='2011271113734'
    elif(Qn==11):
        kext='2012004120508'
    elif(Qn==12):
        kext='2012088054726'
    elif(Qn==13):
        kext='2012179063303'
    elif(Qn==14):
        kext='2012277125453'
    elif(Qn==15):
        kext='2013011073258'
    elif(Qn==16):
        kext='2013098041711'
    elif(Qn==17):
        kext='2013131215648'
    else:
	    print 'Exceed the maximum number of quarter\n'
	    print 'The Default setting of qmax= %d' % (qmax)	
	    exit(1)
    return kext


def gapinfo(qn,usrflag=False,gapfile=''):
#fill this part with Quater information
##* TO BE DONE:	make this function read user supplied gap imformation 

	gap_llc=[]

	if (usrflag):
		qname='q%s' % (qn)
		llcq=cfgp.LC_parse(gapfile,qname)		
	else:	
	#	print 'lcdt',lcdt
		if (qn==0):
			llcq=[]
		if (qn==1):
			llcq=[[985.272,985.474]]
		if (qn==2):	
			llcq=[[1002.05,1005.73],[1014.5186,1019.7214],[1033.1597,1033.3657],[1056.4815,1056.8494],[1059.7723,1059.77236],[1063.30845,1066.4153],[1073.0174,1073.5284],[1079.1901,1080.82537],[1086.60,1086.70],[1088.4085,1089.3283],[1058.08,1058.26]]
		if (qn==3):
			llcq=[[1092.7222,1092.9],[1096.6455,1096.8],[1099.9148,1100.0374],[1113.0536,1113.8301],[1123.0661,1123.9039],[1129.0736,1129.2],[1132.8538,1133.0],[1139.02,1139.05],[1153.9617,1156.9400],[1179.3607,1179.5]]

		if (qn==4):
			llcq=[[1186.2774,1186.2979],[1206.28,1206.8],[1213.9445,1213.9650],[1217.3,1218.73],[1220.1768,1220.1972],[1229.3515,1236.0],[1266.17,1266.20]]
		#llcq=[[1186.50,1186.51],[1200.69,1200.71],[1206.28,1206.8],[1213.14,1213.15],[1217.3,1218.73],[1220.1768,1220.1972],[1226.37,1226.40],[1229.3515,1236.0],[1263.99,1264.01]]
		if (qn==5):
			llcq=[[1309,1310.8],[1319.55559,1319.8],[1321.59895,1321.7],[1324.97174,1326.1],[1336.93,1338.02]]
		if (qn==6):
			llcq=[[1382.45,1382.70],[1395.0,1395.05],[1400.0,1401.9],[1423.99172,1424.1],[1432.0,1434.4]]
		if (qn==7):
			llcq=[[1462.6725,1463.76],[1487.4993,1487.4993],[1492.7916,1494.42],[1522.7473,1525.34],[1542.7314,1542.7314],[1543.4261,1543.9165],[1546.6342,1546.6751]]	
		if (qn==8):
			llcq=[[1567.8647,1568.83],[1572.9935,1572.9935],[1582.6995,1582.6995],[1593.2228,1593.2228],[1593.5702,1598.21],[1613.3908,1613.3908],[1619.8069,1619.8069],[1621.5847,1621.5847]]
		if (qn==9):
			llcq=[[1677.4296,1678.88],[1706.6293,1707.2423],[1719.3390,1719.3390],[1720.7285,1720.7285]]
		if(qn==10):
			llcq=[[1739.3434,1740.3534],[1762.3466,1762.3671],[1769.4560,1771.96],[1799.2008,1799.2108],[1801.7340,1804.66],[1828.3532,1828.3532]]
		if(qn==11):
			llcq=[[1833.7058,1834.92],[1836.6476,1836.6476],[1836.7906,1836.8110],[1840.7947,1840.7948],[1840.9786,1840.9886],[1842.4904,1842.5004],[1853.8286,1853.8386],[1859.6918,1859.7532],[1862.6541,1862.6950],[1864.7787,1866.7],[1870.7849,1870.8849],[1873.7676,1873.7677],[1895.7292,1897.33],[1902.4912,1905.0245],[1915.4434,1915.5434],[1926.1279,1926.2280],[1926.3527,1926.4527],[1928.8654,1928.9472],[1929.1106,1929.7644]]
		if(qn==12):
			llcq=[[1931.9097,1932.],[1935.4642,1935.5643],[1936.7300,1936.8308],[1942.8185,1942.8389],[1945.6989,1945.81],[1949.1922,1951.3168],[1953.7886,1955.21],[1956.0358,1956.1358],[1958.4055,1961.47],[1962.5729,1962.6729],[1966.5360,1966.6360],[1967.6392,1967.7392],[1986.4947,1989.40],[1993.0318,1997.5262],[1998.6088704,1998.7702],[2004.4310,2004.5310],[2009.9263,2010.02],[2011.3358,2011.4358]]
		if(qn==13):
			llcq=[[2015.73,2016.14],[2048.66,2049.90],[2078.02,2079.48]]
		if(qn==14):
			llcq=[[2106.79,2107.79],[2139.00,2141.02],[2170.06,2171.50]]
		if(qn==15):
			llcq=[[2237.40,2238.22],[2246.15,2252.94]]
		if(qn==16):
			llcq=[[2304,2306],[2310,2322.7], [2358,2360.96]]
		if(qn==17):
			llcq=[[2393.83,2394.09],[2414.75,2419.58],[2421.9,2421.95]]
	for x in llcq:
		gap=sp.array(x)
		gap_llc.append(gap)
		
	return gap_llc
					
				
def main():
	return
if __name__=='__main__':
	main()
