import numpy as np
import functions1 as f

valid=np.array([])
predict=np.array([])
weights=np.array([])
losses=np.array([])

left=-2 #estremi asse x
right=2
sigma=0 #rumore
Ntr=100 #numero dati di training
Nval=10 #numero dati di validazione
b=32    #batch size
e=40    #epoche
sd=5    #seed

#-------------------------------------------------------------------------
'''
par=np.array([10,30,60,100]) #vsEpochs
for i in range(len(par)):
    xv,yv,xp,yp,wg,ls=f.lineNet(left,right, sigma, Ntr, Nval, b, par[i], sd)

    if len(predict)==0: predict =np.append(predict, xp)
    predict =np.append(predict, yp)
    weights =np.append(weights, wg)

np.savetxt('../data/weights.dat', weights.reshape(len(par), 2))
np.savetxt('../data/valid.dat', (xv,yv))
np.savetxt('../data/pred.dat', predict.reshape(len(par)+1, Ntr))
np.savetxt('../data/loss.dat', ls)
'''
#-------------------------------------------------------------------------

par=np.array([40,80,160,320]) #vsTrain
Npr=100
for i in range(len(par)):
    xv,yv,xp,yp,wg,ls=f.lineNet(left,right, sigma, par[i], Nval, b, e, sd, Npred=Npr)

    predict =np.append(predict, xp)
    predict =np.append(predict, yp)
    losses =np.append(losses, ls)
    weights =np.append(weights, wg)

np.savetxt('../data/weights.dat', weights.reshape(len(par), 2))
np.savetxt('../data/valid.dat', (xv,yv))
np.savetxt('../data/pred.dat', predict.reshape(2*len(par), Ntr))
np.savetxt('../data/loss.dat', losses.reshape(2*len(par), e))

#-------------------------------------------------------------------------
'''
par=np.array([0.2,0.4,0.6,0.8]) #vsNoise
for i in range(len(par)):
    xv,yv,xp,yp,wg,ls=f.lineNet(left,right, par[i], Ntr, Nval, b, e, sd)#, Npred=Npr)

    
    if len(predict)==0: 
        predict =np.append(predict, xp)
        valid=np.append(valid, xv)
    valid=np.append(valid, yv)
    predict =np.append(predict, yp)
    losses =np.append(losses, ls)
    weights =np.append(weights, wg)

np.savetxt('../data/weights.dat', weights.reshape(len(par), 2))
np.savetxt('../data/valid.dat', valid.reshape(len(par)+1, Nval))
np.savetxt('../data/pred.dat', predict.reshape(len(par)+1, Ntr))
np.savetxt('../data/loss.dat', losses.reshape(2*len(par), e))
'''
#-------------------------------------------------------------------------
print('data printed on files in ../data/')
