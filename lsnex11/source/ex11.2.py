import numpy as np
import functions2 as f

#valid=np.array([])
predict=np.array([])
#weights=np.array([])
losses=np.array([])
Nval=20
Ntr =1000
sd=5

#obj=[[1,1,1],[1,4,1],[4,4,1],[10,4,1],[30,20,1],[4,4,4,1]] #network layers
#obj=['relu','elu','selu','tanh','softmax','sigmoid'] #activation functions
#obj=['sgd','Adam','Adamax','Adagrad','Adadelta','RMSprop'] #optimizers
obj=['mse','msle','mae','kl_divergence','log_cosh','hinge'] #loss functions
metr=['mse','msle','mae','kullback_leibler_divergence','logcosh','hinge'] #metrics

for i in range(len(obj)):
#    xv,yv,xp,yp,w,l=f.Network(-1,1, 0.2, Ntr, Nval, 32, 50, sd, obj[i], 'relu', 'sgd', 'mse')
#    xv,yv,xp,yp,w,l=f.Network(-1,1, 0.2, Ntr, Nval, 32, 50, sd, [25,25,25,1], obj[i], 'sgd', 'mse')
#    xv,yv,xp,yp,w,l=f.Network(-1,1, 0.2, Ntr, Nval, 32, 50, sd, [25,25,25,1], 'relu', obj[i], 'mse')
    xv,yv,xp,yp,w,l=f.Network(-1,1, 0.2, Ntr, Nval, 32, 50, sd, [25,25,25,1], 'selu', 'Adam', obj[i], met=metr[i])

#    valid=np.append(valid, v)
    if len(predict)==0: predict =np.append(predict, xp)
    predict =np.append(predict, yp)
    losses =np.append(losses, l)
#    weights =np.append(weights, w)

#print(np.shape(valid))
#np.savetxt('../data/weights.dat', weights)
np.savetxt('../data/valid.dat', np.array([xv,yv]))
np.savetxt('../data/pred.dat', predict.reshape(len(obj)+1, Ntr))
np.savetxt('../data/losses.dat', losses.reshape(len(obj), 2))
print('data printed on files in ../data/')
