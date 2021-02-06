#import math
import os
import numpy as np

#variational MC per l'identificazione dei migliori parametri
def VMC(mu,sigma,delta):
    inp=[]
    inp.append(str(mu)+'\n') #scrivo su file i parametri che il sorgente C++ dovrà leggere
    inp.append(str(sigma)+'\n')
    inp.append(str(delta)+'\n')
    inp.append(str(M)+'\n')
    inp.append(str(nblk)+'\n')
    inp.append(str(printpsi)+'\n')

    out=open('../data/inputVMC.dat', 'w')
    out.writelines(inp) #scrivo su file
    out.close()

    os.system('./ex1.exe') #eseguo il metropolis
    L,E,err=np.loadtxt('../data/H_av.dat', unpack='True') #leggo i risultati
    del L
    del err
    return E[-1]

M=100000   #step totale
nblk=50    #blocchi
printpsi=1 #se stampo o no la funzione psi
mu=2       #parametri iniziali
sigma=2
def delta(s): #salto per transizioni di x in psi(x)
    return 1.75*s
d=0.75 #salto per transizioni di mu e sigma

t=np.linspace(1,10, 8) #valori di temperatura
tlow=np.linspace(0.1,1, 8, endpoint=False)
beta=1/np.flip( np.concatenate((tlow,t)) ) #abbasso la temperatura man mano
nstepVMC=100

for i in range(len(beta)): #ciclo sulle T
    acc=0
    print('Annealing',i+1,'/',len(beta),end='; ')
    print('temperature T=',1/beta[i],'>>>')
    E=VMC(mu,sigma,delta(sigma))

    for j in range(nstepVMC): #variational MC per ogni T
        muTrial=abs(mu + np.random.uniform(-d,d))
        sigmaTrial=abs(sigma + np.random.uniform(-d,d))
        ETrial=VMC(muTrial,sigmaTrial,delta(sigmaTrial))

        if np.random.uniform() < np.exp( -beta[i] * (ETrial-E) ): #acceptance
            mu=muTrial
            sigma=sigmaTrial
            E=ETrial
            acc+=1

        if i==0: #saving if optimal 
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
