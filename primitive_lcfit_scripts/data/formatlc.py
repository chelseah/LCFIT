from numpy import *
import matplotlib.pyplot as plt
import functions

data = loadtxt("kplr006603043-2011145075126_slc.tab")

data_new = []
for i in data:
    if i[0] > 1691.5 and i[0] < 1693.1:
        data_new.append(i)
data = array(data_new)

mag = data[:,3]
flux = 10**((mag-median(mag))/-2.5)

o = open("lc2.dat","w")
output = [data[:,0],flux,data[:,4]]
output = transpose(output)
functions.write_table(output,o)
o.close()

plt.scatter(data[:,0],flux)
plt.show()
