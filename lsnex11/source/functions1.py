import numpy as np
import tensorflow as tf
from tensorflow import keras
#from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense#, Activation
#from tensorflow.keras.utils import get_custom_objects
#from tensorflow.keras import backend as K

def line(x): return 4*x + 3

def lineNet(left,right,sigma,Ntrain,Nvalid,batch,ep,seed, Npred=-1):
    if Npred<0: Npred=Ntrain
    np.random.seed(seed)
    x_train=np.random.uniform(left,right,Ntrain)
    x_valid=np.random.uniform(left,right,Nvalid)
    x_valid.sort()
#    y_target= line(x_valid)

    y_train=np.random.normal(line(x_train), sigma)
    y_valid=np.random.normal(line(x_valid), sigma)

    model= tf.keras.Sequential()
    model.add(Dense(1, input_shape=(1,)))
    model.compile(optimizer='sgd', loss='mse', metrics=['mse'])
    model.summary()

    x_pred=np.random.uniform(left,right,Npred)
    history=model.fit(x=x_train, y=y_train, batch_size=batch, epochs=ep, shuffle=True, validation_data=(x_valid, y_valid))
    #score=model.evaluate(x_valid, y_valid, batch_size=32, verbose=1)
    #score_target=model.evaluate(x_valid, y_target, batch_size=32, verbose=1)
    #print(np.shape(model.predict(x_pred)))

    return x_valid, y_valid, x_pred, model.predict(x_pred), model.get_weights(), np.array([history.history['loss'], history.history['val_loss']])
