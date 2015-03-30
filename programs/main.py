import numpy as np
import random
import timeit

start = timeit.default_timer()

random.seed(0)

M = 5      ### M times the number of data points

pos = np.loadtxt('../data/cells/0.cvs', delimiter=',', usecols=[0,1,2])
N=pos[:,0].shape[0]

Ntotal = 12.651e6
Nbox = int(Ntotal/(25**3))
Nrand=M*Nbox
print 'N data   =', N
print 'N random =', Nrand


dist = 10**(np.linspace(-2.125,2.125,18)) ## Bins
DIST = np.logspace(-2,2,17)                   ## Center of Bins
counterDD = np.zeros(17)
counterDR = np.zeros(17)
counterRR = np.zeros(17)

#### Create Random Field
RandomField = np.zeros([Nrand,3])
for i in range(Nrand):
    for k in range(3):
        RandomField[i,k] = random.random()*10 # Boxlength

#print RandomField.shape

#### Counting Pairs in Random Field
for i in range(Nrand-1):
    for j in range(i+1,Nrand):
        r = 0.0
        for k in range(3):
            r+=(RandomField[i,k]-RandomField[j,k])**2
        r = r**0.5
        m = 0
        for m in range(counterDD.size):
            if (dist[m]<r) and (r<dist[m+1]):
                counterRR[m] +=1
counterRR = counterRR/( Nrand*(Nrand-1) ) #### For 1 box, already divided by 2

print counterRR, '\n'



### Counting Pairs in Data-Data

for i in range(N-1):
    for j in range(i+1,N):
        r = 0.0
        for k in range(3):
            r+=(pos[i,k]-pos[j,k])**2
        r = r**0.5
        m = 0
        for m in range(counterDD.size):
            if (dist[m]<r) and (r<dist[m+1]):
                counterDD[m] +=1
counterDD = counterDD*1.0/(N*(N-1)) #### For 1 box, already divided by 2
print counterDD


### Counting pairs Data Random
for i in range(N-1):
    for j in range(0,Nrand):
        r = 0.0
        for k in range(3):
            r+=(pos[i,k]-RandomField[j,k])**2
        r = r**0.5
        m = 0
        for m in range(counterDR.size):
            if (dist[m]<r) and (r<dist[m+1]):
                counterDR[m] +=1
counterDR = counterDR*1.0/(N*Nrand)
print counterDR

CorrFunc = (counterDD - 2*counterDR + counterRR)/ counterRR

print 'Correlation Function:\n'
print CorrFunc

stop = timeit.default_timer()
print stop-start, 'seconds\n'
