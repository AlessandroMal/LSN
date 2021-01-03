#import math
import os
import numpy as np

def VMC(mu,sigma,delta):
    inp=[]
    inp.append(str(mu)+'\n')
    inp.append(str(sigma)+'\n')
    inp.append(str(delta)+'\n')
    inp.append(str(M)+'\n')
    inp.append(str(nblk)+'\n')
    inp.append(str(printpsi)+'\n')

    out=open('../data/input.dat', 'w')
    out.writelines(inp)
    out.close()

    os.system('./ex1.exe')
    L,E,err=np.loadtxt('../data/H_av.dat', unpack='True')
    return E[-1]

M=100000
nblk=50
printpsi=1
mu=2
sigma=2
def delta(s): #salto per transizioni di x in psi(x)
    return 2*s
d=1 #salto per transizioni di mu e sigma

t=np.linspace(1,20, 4)
tlow=np.linspace(0.1,1, 4, endpoint='False')
beta=1/np.flip( np.concatenate((tlow,t)) ) #abbasso la temperatura man mano
nstepVMC=3

for i in range(len(beta)):
    acc=0
    print('Annealing',i+1,'/',len(beta),'...')
    print('temperature T=',1/beta[i])

    for j in range(nstepVMC):
        E=VMC(mu,sigma,delta(sigma))

        muTrial=mu + np.random.uniform(-d,d)
        sigmaTrial=sigma + np.random.uniform(-d,d)
        ETrial=VMC(muTrial,sigmaTrial,delta(sigmaTrial))

        if np.random.uniform() < np.exp( -beta[i] * (ETrial-E) ):
            mu=muTrial
            sigma=sigmaTrial
            acc+=1

        if i==0:
            Emin=E
            mumin=mu
            sigmamin=sigma
            os.system('cp ../data/psi.out ../data/psiopt.out')
            os.system('cp ../data/H_av.dat ../data/Hmin.dat')
        else:
            if E<Emin:
                Emin=E
                mumin=mu
                sigmamin=sigma
                os.system('cp ../data/psi.out ../data/psiopt.out')
                os.system('cp ../data/H_av.dat ../data/Hmin.dat')
    print('Acceptance rate for annealing: ',acc/nstepVMC)

print('GS energy:       ',Emin)
print('Optimized mu:    ',mumin)
print('Optimized sigma: ',sigmamin)

out=open('../data/input.dat', 'a')
out.write('\t mu\n \t sigma\n \t delta\n \t M\n \t nblk\n \t printpsi')
out.close()
