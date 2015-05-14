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
Random = np.loadtxt('../data/RandomFieldCatalog03.csv',
                            delimiter=',', usecols=[0,1,2])
Nrandom = Random[:,0].shape[0]


###################### 3D Correlation Functions

def Counting_DD():
    """
    Counting_DD(none)

    Makes the projected pair counting in DATA-DATA
    Returns a square matrix where each element is a bin of perpendicular
    and paralell distance given by DIST.
    """
    counter   = np.zeros(((Bins-1),(Bins-1)))

    print '\nCounting DATA-DATA'
    T_start_Countig = timeit.default_timer()

    matrix = np.zeros([Ndata,8])    ### x1 y1 z1  x2 y2 z2  R_perp   R_parallel

    for k in range(3):
            matrix[:,k]   = Data[:,k]

    for line in Data:
        aux = 0
        for k in range(3):
            matrix[:,k+3] = line[k]
        for k in range(2):
            aux += (matrix[:,k]-matrix[:,k+3])**2
        matrix[:,6] = np.sqrt(aux)       #### R_perpendicular

        matrix[:,7] = np.abs(matrix[:,2]-matrix[:,2+3])    #### R_parallel

        for m in range(Bins-1):
            index = np.where((matrix[:,6]>=dist[m])&(matrix[:,6]<dist[m+1]))[0]   ### R_perp
            for n in range(Bins-1):
                index2 = np.where((matrix[index,7]>=dist[n])&(matrix[index,7]<dist[n+1]))[0] ### R_paralell
                counter[m,n] += index2.size

    counter = counter/( Ndata*(Ndata-1)/2.0 )

    T_stop_Counting = timeit.default_timer()
    print 'DD counting time: ', T_stop_Counting - T_start_Countig,'s'
    print '   done\n'
    #print counter

    file = open('../results/DD.dat', 'w')  #a append, w write, r read, r+ read+write
    file.write('#'+str(timeit.time.ctime())+'\n')
    file.write("# DD time execution:"+str(T_stop_Counting - T_start_Countig)+'s\n')
    file.write("# Lines: R_perpendicular, Columns: R_paralell \n")
    for m in range(Bins-1):
        for n in range(Bins-2):
            file.write(str(counter[m,n])+',')
        file.write(str(counter[m,Bins-2])+'\n')
    file.close()

    return counter


def Counting_RR():
    """
    Counting_RR(none)

    Makes the projected pair counting in RANDOM-RANDOM
    Returns a square matrix where each element is a bin of perpendicular
    and paralell distance given by DIST.
    """
    counter   = np.zeros(((Bins-1),(Bins-1)))

    print '\nCounting RANDOM-RANDOM'
    T_start_Countig = timeit.default_timer()

    matrix = np.zeros([Nrandom,8])    ### x1 y1 z1  x2 y2 z2  R_perp   R_parallel

    for k in range(3):
            matrix[:,k]   = Random[:,k]

    for line in Random:
        aux = 0
        for k in range(3):
            matrix[:,k+3] = line[k]
        for k in range(2):
            aux += (matrix[:,k]-matrix[:,k+3])**2
        matrix[:,6] = np.sqrt(aux)       #### R_perpendicular

        matrix[:,7] = np.abs(matrix[:,2]-matrix[:,2+3])    #### R_parallel

        for m in range(Bins-1):
            index = np.where((matrix[:,6]>=dist[m])&(matrix[:,6]<dist[m+1]))[0]   ### R_perp
            for n in range(Bins-1):
                index2 = np.where((matrix[index,7]>=dist[n])&(matrix[index,7]<dist[n+1]))[0] ### R_paralell
                counter[m,n] += index2.size

    counter = counter/( Nrandom*(Nrandom-1)/2.0 )

    T_stop_Counting = timeit.default_timer()
    print 'RR counting time: ', T_stop_Counting - T_start_Countig,'s'
    print '   done\n'
    #print counter

    file = open('../results/RR03.dat', 'w')  #a append, w write, r read, r+ read+write
    file.write('#'+str(timeit.time.ctime())+'\n')
    file.write("# RR time execution:"+str(T_stop_Counting - T_start_Countig)+'s\n')
    file.write("# Lines: R_perpendicular, Columns: R_paralell \n")
    for m in range(Bins-1):
        for n in range(Bins-2):
            file.write(str(counter[m,n])+',')
        file.write(str(counter[m,Bins-2])+'\n')
    file.close()

    return counter


def Counting_DR():
    """
    Counting_DR(none)

    Makes the projected pair counting in DATA-RANDOM
    Returns a square matrix where each element is a bin of perpendicular
    and paralell distance given by DIST.
    """
    counter   = np.zeros(((Bins-1),(Bins-1)))

    print '\nCounting DATA-RANDOM'
    T_start_Countig = timeit.default_timer()

    matrix = np.zeros([Ndata,8])    ### x1 y1 z1  x2 y2 z2  R_perp   R_parallel

    for k in range(3):
            matrix[:,k]   = Data[:,k]

    for line in Random:
        aux = 0
        for k in range(3):
            matrix[:,k+3] = line[k]
        for k in range(2):
            aux += (matrix[:,k]-matrix[:,k+3])**2
        matrix[:,6] = np.sqrt(aux)       #### R_perpendicular

        matrix[:,7] = np.abs(matrix[:,2]-matrix[:,2+3])    #### R_parallel

        for m in range(Bins-1):
            index = np.where((matrix[:,6]>=dist[m])&(matrix[:,6]<dist[m+1]))[0]   ### R_perp
            for n in range(Bins-1):
                index2 = np.where((matrix[index,7]>=dist[n])&(matrix[index,7]<dist[n+1]))[0] ### R_paralell
                counter[m,n] += index2.size

    counter = counter/( Ndata*Nrandom )

    T_stop_Counting = timeit.default_timer()
    print 'DR counting time: ', T_stop_Counting - T_start_Countig,'s'
    print '   done\n'
    #print counter

    file = open('../results/DR03.dat', 'w')  #a append, w write, r read, r+ read+write
    file.write('#'+str(timeit.time.ctime())+'\n')
    file.write("# DR time execution:"+str(T_stop_Counting - T_start_Countig)+'s\n')
    file.write("# Lines: R_perpendicular, Columns: R_paralell \n")
    for m in range(Bins-1):
        for n in range(Bins-2):
            file.write(str(counter[m,n])+',')
        file.write(str(counter[m,Bins-2])+'\n')
    file.close()

    return counter

##################################### main program

#A = Counting_DD()
#B = Counting_DR()
C = Counting_RR()
