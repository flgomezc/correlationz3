import numpy as np
RandomCat = open('../data/RandomFieldCatalog.csv', 'w')
Nrand = 234271
np.random.seed()

for i in range(Nrand):
    RandomCat.write(
                    str(250*np.random.rand())+','+
                    str(250*np.random.rand())+','+
                    str(250*np.random.rand())+'\n'
    )

RandomCat.close()
