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


def Counting_DD():
    counterDD = np.zeros(Bins-1)
    print '\nCounting DATA-DATA'
    T_start_Countig_DD = timeit.default_timer()
    matrix = np.zeros([Ndata,7])
    for k in range(3):
            matrix[:,k]   = Data[:,k]
    counterDD = np.zeros(Bins-1)
    for line in Data:
        aux = 0
        for k in range(3):
            matrix[:,k+3] = line[k]
            aux += (matrix[:,k]-matrix[:,k+3])**2

        matrix[:,6] = np.sqrt(aux)
        for m in range(counterDD.size):
            index = np.where((matrix[:,6]>=dist[m])&(matrix[:,6]<dist[m+1]))[0]
            counterDD[m] += index.shape
    counterDD = counterDD/( 2*Nrand*(Nrand-1) )
    T_stop_Counting_DD = timeit.default_timer()
    print 'DD counting time: ', T_stop_Counting_DD - T_start_Countig_DD,'s\n'
    print '\n   DD=\n'
    print counterDD
    return counterDD

def Counting_RR():
    counterRR = np.zeros(Bins-1)
    print '\nCounting RANDOM-RANDOM'
    T_start_Countig_RR = timeit.default_timer()
    matrix = np.zeros([Nrand,7])
    for k in range(3):
            matrix[:,k]   = RandomField[:,k]
    #counterRR = np.zeros(Bins-1)
    for line in RandomField:
        aux = 0
        for k in range(3):
            matrix[:,k+3] = line[k]
            aux += (matrix[:,k]-matrix[:,k+3])**2

        matrix[:,6] = np.sqrt(aux)
        for m in range(counterRR.size):
            index = np.where((matrix[:,6]>=dist[m])&(matrix[:,6]<dist[m+1]))[0]
            counterRR[m] += index.shape
    counterRR = counterRR/( 2*Nrand*(Nrand-1) )
    T_stop_Counting_RR = timeit.default_timer()
    print 'RR counting time: ', T_stop_Counting_RR - T_start_Countig_RR,'s\n'
    print '\n   RR=\n'
    print counterRR
    return counterRR

def Counting_DR():
    counterDR = np.zeros(Bins-1)
    print '\nCounting DATA-RANDOM'
    T_start_Countig_DR = timeit.default_timer()
    matrix = np.zeros([Ndata,7])
    for k in range(3):
            matrix[:,k]   = Data[:,k]
    counterDR = np.zeros(Bins-1)
    for line in RandomField:
        aux = 0
        for k in range(3):
            matrix[:,k+3] = line[k]
            aux += (matrix[:,k]-matrix[:,k+3])**2

        matrix[:,6] = np.sqrt(aux)
        for m in range(counterDR.size):
            index = np.where((matrix[:,6]>=dist[m])&(matrix[:,6]<dist[m+1]))[0]
            counterDR[m] += index.shape
    counterDR = counterDR/( 2*Ndata*Nrand)
    T_stop_Counting_DR = timeit.default_timer()
    print 'DR counting time: ', T_stop_Counting_DR - T_start_Countig_DR,'s\n'
    print '\n   DR=\n'
    print counterDR
    return counterDR
