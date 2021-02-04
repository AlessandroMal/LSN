import numpy as np
import functions3 as f

Ntr=1000
Nval=100
Npre=4000


xv,yv,xp,yp,w,s=f.Network(-1.5,1.5, 0.25, Ntr, Nval, 32, 150, 5, [150,50,30,1], 'relu', 'sgd', 'mse', Npred=Npre)

val=np.concatenate((xv[:,0],xv[:,1],yv))
pre=np.concatenate((xp[:,0],xp[:,1],yp[:,0]))

#print(np.shape(xv), np.shape(yv), np.shape(xp), np.shape(yp))
#np.savetxt('../data/weights.dat', w)
np.savetxt('../data/valid.dat', val.reshape(3,Nval))
np.savetxt('../data/pred.dat', pre.reshape(3,Npre))
np.savetxt('../data/score.dat', s)
print('data printed on files in ../data/')
