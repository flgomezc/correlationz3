import numpy as np
import random
import timeit

Catalog = np.loadtxt('../data/BolshoiP_bdmV.csv', skiprows=1, delimiter=',',
                        usecols=[3,4,5,6])
index = np.where(Catalog[:,3]>5*10**10)
SmallCatalog = np.array(Catalog[index])
Data = open('../data/SmallCatalog.csv', 'w')
for line in SmallCatalog:
    Data.write(str(line[0])+','+str(line[1])+','+str(line[2])+'\n')
Data.close()


### Create Random Field Catalog

RFCatalog = open('../data/RandomFieldCatalog.cvs', 'w')
random.seed(0)
BoxLength = 250 # Mpc h^-1
M = 1      ### M times the number of data points
Ndata=SmallCatalog[:,0].shape[0]
Nrand=int(M*Ndata)

RandomField = np.zeros([Nrand,3])
for line in RandomField:
    for k in range(3):
        line[k] = random.random()*BoxLength
        RFCatalog.write(str(line[k])+',')
    RFCatalog.write('\n')
RFCatalog.close()
