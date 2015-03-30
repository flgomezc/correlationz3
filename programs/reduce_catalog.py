import numpy as np
import random
import timeit

Catalog = np.loadtxt('../data/BolshoiP_bdmV.csv', skiprows=1, delimiter=',', usecols=[3,4,5,6])

index = np.where(Catalog[:,3]>5*10**10)

SmallCatalog = np.array(Catalog[index])

Data = open('../data/SmallCatalog.csv', 'w')
for line in SmallCatalog:
    Data.write(str(line[0])+','+str(line[1])+','+str(line[2])+'\n')
Data.close()
