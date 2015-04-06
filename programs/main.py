import numpy as np
import random
import timeit
random.seed(0)

T_start_global = timeit.default_timer()

RunningLog = open('running.log', 'w')

T_start_loading = timeit.default_timer()
Data = np.loadtxt('../data/SmallCatalog.csv', delimiter=',', usecols=[0,1,2])
T_stop_loading = timeit.default_timer()

RunningLog.write('Correlation Function At z=3 Log\n\n')
RunningLog.write('Loading Data: '+str(T_stop_loading - T_start_loading)+'s\n')

BoxLength = 250 # Mpc h^-1


Ndata=Data[:,0].shape[0]
RunningLog.write('N data   ='+str(Ndata)+'\n')

Bins = 16
dist = 10**(np.linspace(-1.10,2.10,Bins)) ## Bins
DIST = np.logspace(-1,2,15)                   ## Center of Bins
counterDD = np.zeros(Bins-1)
counterDR = np.zeros(Bins-1)
counterRR = np.zeros(Bins-1)

RunningLog.write("Distance Bins =\n"+str(dist)+'\n\n')



T_start_RandomField = timeit.default_timer()
#### Create Random Field
RandomField = np.loadtxt('../data/RandomFieldCatalog.csv', delimiter=',', usecols=[0,1,2])
Nrand = RandomField[:,0].shape[0]
RunningLog.write('N random ='+str(Nrand)+'\n\n')
print 'RandomField', RandomField.shape
T_stop_RandomField = timeit.default_timer()
RunningLog.write('Loading RandomField: '+str(T_stop_RandomField - T_start_RandomField )+'s\n')
RunningLog.write('RandomField.shape =' +str(RandomField.shape)+'\n')


T_start_Countig_RR = timeit.default_timer()
matrix = np.zeros([Nrand,7])
for k in range(3):
        matrix[:,k]   = RandomField[:,k]
counterRR2 = np.zeros(Bins-1)
for line in RandomField:
    aux = 0
    for k in range(3):
        matrix[:,k+3] = line[k]
        aux += (matrix[:,k]-matrix[:,k+3])**2

    matrix[:,6] = np.sqrt(aux)
    for m in range(counterRR2.size):
        index = np.where((matrix[:,6]>=dist[m])&(matrix[:,6]<dist[m+1]))[0]
        counterRR2[m] += index.shape
counterRR2 = counterRR2/( 2*Nrand*(Nrand-1) )
T_stop_Counting_RR = timeit.default_timer()

RunningLog.write('Random-Random counting\nRR = \n'+str(counterRR2)+'\n\n')
RunningLog.write('Random-Random counting time = '+
                  str(T_stop_Counting_RR - T_start_Countig_RR)+'s\n')
print 'RR counting time: ', T_stop_Counting_RR - T_start_Countig_RR,'s\n'



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

RunningLog.write('Data-Data counting\nDD = \n'+str(counterDD)+'\n\n')
RunningLog.write('Data-Data counting time = '+
                  str(T_stop_Counting_DD - T_start_Countig_DD)+'s\n')
print 'DD counting time: ', T_stop_Counting_DD - T_start_Countig_DD,'s\n'


T_start_Countig_DR = timeit.default_timer()
matrix = np.zeros([Nrand,7])
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

RunningLog.write('Data-Random counting\nDR = \n'+str(counterDR)+'\n\n')
RunningLog.write('Data-Random counting time = '+
                  str(T_stop_Counting_DR - T_start_Countig_DR)+'s\n\n\n')
print 'DR counting time: ', T_stop_Counting_DR - T_start_Countig_DR,'s\n'



CorrFunc = (counterDD - 2*counterDR + counterRR)/counterRR
RunningLog.write('Correlation Function\nDR = \n'+str(CorrFunc)+'\n\n')

T_stop_global = timeit.default_timer()
RunningLog.write('Global running time = '+
                  str(T_stop_global - T_start_global)+'s\n\n\n')
print 'Global running time: ', T_stop_global - T_start_global,'s\n'

RunningLog.close()
