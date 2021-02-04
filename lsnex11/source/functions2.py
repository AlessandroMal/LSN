import numpy as np
import tensorflow as tf
from tensorflow import keras
#from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
#from tensorflow.keras.utils import get_custom_objects
#from tensorflow.keras import backend as K

def f_targ(x): return 4 - 3*x - 2*x**2 + 3*x**3

def Network(left,right,sigma,Ntrain,Nvalid,batch,ep, seed, N_neu, f_act, opt, los, met=''):
    if met=='': met=los
    np.random.seed(seed)
    x_train=np.random.uniform(left,right,Ntrain)
    x_valid=np.random.uniform(left,right,Nvalid)
    x_valid.sort()
    #y_target= f_targ(x_valid)

    y_train=np.random.normal(f_targ(x_train), sigma)
    y_valid=np.random.normal(f_targ(x_valid), sigma)

    model= tf.keras.Sequential()
    for i in range(len(N_neu)):
        if i==0: model.add(Dense(units=N_neu[i], input_shape=(1,), activation=f_act))
        else:
            if i==len(N_neu)-1: model.add(Dense(units=N_neu[i]))
            else: model.add(Dense(units=N_neu[i], activation=f_act))

    model.compile(optimizer=opt, loss=los, metrics=[met])
    model.summary()

    x_pred=np.random.uniform(left-(right-left)/4,right+(right-left)/4,Ntrain)
    history=model.fit(x=x_train, y=y_train, batch_size=batch, epochs=ep, shuffle=True, validation_data=(x_valid, y_valid))

    return x_valid, y_valid, x_pred, model.predict(x_pred), model.get_weights(), np.array([history.history['loss'][-1], history.history['val_loss'][-1]])
