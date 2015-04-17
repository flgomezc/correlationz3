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

    counterDD = counterDD/( Ndata*(Ndata-1) )
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

    counterRR = counterRR/( Nrand*(Nrand-1) )
    T_stop_Counting_RR = timeit.default_timer()
    print 'RR counting time: ', T_stop_Counting_RR - T_start_Countig_RR,'s'
    print '   RR='
    print counterRR

    file = open('../results/RR.dat', 'w')
    for i in range(Bins-2):
        file.write(str(counterRR[i])+',')
    file.write(str(counterRR[Bins-2])+'\n')
    file.close()
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
    counterDR = counterDR/( Ndata*Nrand)
    T_stop_Counting_DR = timeit.default_timer()
    print 'DR counting time: ', T_stop_Counting_DR - T_start_Countig_DR,'s'
    print '   DR='
    print counterDR

    file = open('../results/DR.dat', 'w')
    for i in range(Bins-2):
        file.write(str(counterDR[i])+',')
    file.write(str(counterDR[Bins-2])+'\n')
    file.close()
    return counterDR


################################## Projected Correlation Functions



def Proj_Counting_DD():
    counterDD = np.zeros(Bins-1)
    print '\nProjected XY pair Counting DATA-DATA'
    T_start_Countig_DD = timeit.default_timer()
    matrix = np.zeros([Ndata,5])

    for k in range(2):
            matrix[:,k]   = Data[:,k]
    counterDD = np.zeros(Bins-1)

    for line in Data:
        aux = 0
        for k in range(2):
            matrix[:,k+2] = line[k]
            aux += (matrix[:,k]-matrix[:,k+2])**2

        matrix[:,4] = np.sqrt(aux)
        for m in range(counterDD.size):
            index = np.where((matrix[:,4]>=dist[m])&(matrix[:,4]<dist[m+1]))[0]
            counterDD[m] += index.shape

    counterDD = counterDD/( Ndata*(Ndata-1) )
    T_stop_Counting_DD = timeit.default_timer()
    print 'p_DD counting time: ', T_stop_Counting_DD - T_start_Countig_DD,'s'
    print '   p_DD='
    print counterDD

    file = open('../results/p_DD.dat', 'w')
    for i in range(Bins-2):
        file.write(str(counterDD[i])+',')
    file.write(str(counterDD[Bins-2])+'\n')
    file.close()
    return counterDD

def Proj_Counting_RR():
    counterRR = np.zeros(Bins-1)
    print '\nProjected XY pair Counting RANDOM-RANDOM'
    T_start_Countig_RR = timeit.default_timer()
    matrix = np.zeros([Nrand,5])
    for k in range(2):
            matrix[:,k]   = RandomField[:,k]
    #counterRR = np.zeros(Bins-1)
    for line in RandomField:
        aux = 0
        for k in range(2):
            matrix[:,k+2] = line[k]
            aux += (matrix[:,k]-matrix[:,k+2])**2

        matrix[:,4] = np.sqrt(aux)
        for m in range(counterRR.size):
            index = np.where((matrix[:,4]>=dist[m])&(matrix[:,4]<dist[m+1]))[0]
            counterRR[m] += index.shape

    counterRR = counterRR/( Nrand*(Nrand-1) )
    T_stop_Counting_RR = timeit.default_timer()
    print 'p_RR counting time: ', T_stop_Counting_RR - T_start_Countig_RR,'s'
    print '   p_RR='
    print counterRR

    file = open('../results/p_RR.dat', 'w')
    for i in range(Bins-2):
        file.write(str(counterRR[i])+',')
    file.write(str(counterRR[Bins-2])+'\n')
    file.close()
    return counterRR

def Proj_Counting_DR():
    counterDR = np.zeros(Bins-1)
    print '\nProjected XY pair Counting DATA-RANDOM'
    T_start_Countig_DR = timeit.default_timer()
    matrix = np.zeros([Ndata,5])
    for k in range(2):
            matrix[:,k]   = Data[:,k]
    counterDR = np.zeros(Bins-1)
    for line in RandomField:
        aux = 0
        for k in range(2):
            matrix[:,k+2] = line[k]
            aux += (matrix[:,k]-matrix[:,k+2])**2

        matrix[:,4] = np.sqrt(aux)
        for m in range(counterDR.size):
            index = np.where((matrix[:,4]>=dist[m])&(matrix[:,4]<dist[m+1]))[0]
            counterDR[m] += index.shape
    counterDR = counterDR/( Ndata*Nrand)
    T_stop_Counting_DR = timeit.default_timer()
    print 'p_DR counting time: ', T_stop_Counting_DR - T_start_Countig_DR,'s'
    print '   p_DR='
    print counterDR

    file = open('../results/p_DR.dat', 'w')
    for i in range(Bins-2):
        file.write(str(counterDR[i])+',')
    file.write(str(counterDR[Bins-2])+'\n')
    file.close()
    return counterDR
