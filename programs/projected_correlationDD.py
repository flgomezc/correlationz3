
import numpy as np
import random
import timeit

BoxLength = 250 # Mpc h^-1

Bins = 16
dist = 10**(np.linspace(-1.10,2.10,Bins)) ## Bins
DIST = np.logspace(-1,2,15)               ## Center of Bins
print 'DIST = '
print DIST


# Loading Catalogs
print '\nLoading Small Catalog'
Data = np.loadtxt('../data/SmallCatalog.csv',
                            delimiter=',', usecols=[0,1,2])
Ndata=Data[:,0].shape[0]

print 'Loading Random Field Catalog'
RandomField = np.loadtxt('../data/RandomFieldCatalog.csv',
                            delimiter=',', usecols=[0,1,2])
Nrand = RandomField[:,0].shape[0]


###################### 3D Correlation Functions

def Counting_DD():
    counterDD   = np.zeros(((Bins-1),(Bins-1)))
    
    print '\nCounting DATA-DATA'
    T_start_Countig_DD = timeit.default_timer()
    
    matrix = np.zeros([Ndata,8])    ### x1 y1 z1  x2 y2 z2  R   R_parallel

    for k in range(3):
            matrix[:,k]   = Data[:,k]

    for line in Data:
        aux = 0
        for k in range(3):
            matrix[:,k+3] = line[k]
        for k in range(2):
            aux += (matrix[:,k]-matrix[:,k+3])**2
        matrix[:,6] = np.sqrt(aux)
        
        matrix[:,7] = np.abs(matrix[:,2]-matrix[:,2+3])
        
        for m in range(Bins-1):
            index = np.where((matrix[:,6]>=dist[m])&(matrix[:,6]<dist[m+1]))[0]
            for n in range(Bins-1):
                index2 = np.where((matrix[index,7]>=dist[n])&(matrix[index,7]<dist[n+1]))[0]
                counterDD[m,n] += index2.size
                
    counterDD = counterDD/( Ndata*(Ndata-1)/2.0 )
    
    T_stop_Counting_DD = timeit.default_timer()
    print 'DD counting time: ', T_stop_Counting_DD - T_start_Countig_DD,'s'
    print '   DD='
    print counterDD

    file = open('../results/DD.dat', 'w')
    for i in range(Bins-2):
        file.write(str(counterDD[i])+',')
    file.write(str(counterDD[Bins-2])+'\n')
    file.close()
    return counterDD


##################################### main program

Counting_DD()



