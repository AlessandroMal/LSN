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

    out=open('../data/inputVMC.dat', 'w')
    out.writelines(inp)
    out.close()

    os.system('./ex1.exe')
    L,E,err=np.loadtxt('../data/H_av.dat', unpack='True')
    del L
    del err
    return E[-1]

M=100000
nblk=50
printpsi=1
mu=2
sigma=2
def delta(s): #salto per transizioni di x in psi(x)
    return 1.75*s
d=0.75 #salto per transizioni di mu e sigma

t=np.linspace(1,10, 8)
tlow=np.linspace(0.1,1, 8, endpoint=False)
beta=1/np.flip( np.concatenate((tlow,t)) ) #abbasso la temperatura man mano
nstepVMC=100

for i in range(len(beta)):
    acc=0
    print('Annealing',i+1,'/',len(beta),end='; ')
    print('temperature T=',1/beta[i],'>>>')
    E=VMC(mu,sigma,delta(sigma))

    for j in range(nstepVMC):
        muTrial=abs(mu + np.random.uniform(-d,d))
        sigmaTrial=abs(sigma + np.random.uniform(-d,d))
        ETrial=VMC(muTrial,sigmaTrial,delta(sigmaTrial))

        if np.random.uniform() < np.exp( -beta[i] * (ETrial-E) ):
            mu=muTrial
            sigma=sigmaTrial
            E=ETrial
            acc+=1

        if i==0:
            Emin=E
            mumin=mu
            sigmamin=sigma
            os.system('cp ../data/psi.out ../data/psiopt.out')
            os.system('cp ../data/H_av.dat ../data/Hmin.dat')
            param=open('../data/psiopt.out', 'a')
            param.write(str(mumin)+' '+str(sigmamin))
            param.close()
        else:
            if E<Emin:
                Emin=E
                mumin=mu
                sigmamin=sigma
                os.system('cp ../data/psi.out ../data/psiopt.out')
                os.system('cp ../data/H_av.dat ../data/Hmin.dat')
                param=open('../data/psiopt.out', 'a')
                param.write(str(mumin)+' '+str(sigmamin))
                param.close()
    print('mu_opt=',mumin,'  sigma_opt=',sigmamin,'  E_min=',Emin)
    print('Acceptance rate for annealing: ',acc/nstepVMC)
    print('-------------------------------------------')

print('GS energy:       ',Emin)
print('Optimized mu:    ',mumin)
print('Optimized sigma: ',sigmamin)

out=open('../data/inputVMC.dat', 'a')
out.write('\t mu\n \t sigma\n \t delta\n \t M\n \t nblk\n \t printpsi')
out.close()
