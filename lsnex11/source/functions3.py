import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers import Dense, Activation
#from tensorflow.keras.utils import get_custom_objects
#from tensorflow.keras import backend as K

def f_targ(x,y): return np.sin(x**2 + y**2)

def Network(left,right,sigma,Ntrain,Nvalid,batch,ep, seed, N_neu, f_act, opt, los, Npred=-1):
    if Npred<0: Npred=Ntrain
    np.random.seed(seed)

    
    x_train=np.random.uniform(left,right,2*Ntrain)
    x_valid=np.random.uniform(left,right,2*Nvalid)
    x_train=x_train.reshape(Ntrain,2)
    x_valid=x_valid.reshape(Nvalid,2)

    #x_valid.sort()
    #y_target= f_targ(x_valid)

    y_train=np.random.normal(f_targ(x_train[:,0],x_train[:,1]), sigma)
    y_valid=np.random.normal(f_targ(x_valid[:,0],x_valid[:,1]), sigma)

    model= tf.keras.Sequential()
    for i in range(len(N_neu)):
        if i==0: model.add(Dense(units=N_neu[i], input_shape=(2,), activation=f_act))
        else:
            if i==len(N_neu)-1: model.add(Dense(units=N_neu[i]))
            else: model.add(Dense(units=N_neu[i], activation=f_act))

    model.compile(optimizer=opt, loss=los, metrics=[los])
    model.summary()

    history=model.fit(x=x_train, y=y_train, batch_size=batch, epochs=ep, shuffle=True, validation_data=(x_valid, y_valid))

    x_pred=np.random.uniform(left,right,2*Npred)
    x_pred=x_pred.reshape(Npred,2)
    score=model.evaluate(x_pred, f_targ(x_pred[:,0],x_pred[:,1]), batch_size=32, verbose=1)

    return x_valid, y_valid, x_pred, model.predict(x_pred), model.get_weights(), score
