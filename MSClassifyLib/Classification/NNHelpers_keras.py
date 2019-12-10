"""
This file contains all custom layers and losses of the keras models. To guarantee correct
saving and loading of your own models, add your custom layers and losses here and create an 
entry in the @get_custom_layer_dict() and @get_custom_losses_dict() functions below.
"""

from tensorflow.keras.models import Model
from tensorflow.keras import backend as K
from tensorflow.keras.layers import Activation, Add, Conv2D
from tensorflow.keras.layers import Dropout, Dense, Flatten, LocallyConnected2D, Input

from tensorflow.keras.layers import BatchNormalization
import numpy as np

#####################################################################
############### Section for Iteration tools #########################
#####################################################################

# Creates a simple iterable, can be used as a batch iterator
def batchGenerator(data, batch_size=64):
    for i in range(
        np.ceil(data.shape[0] / np.float32(batch_size)
               ).astype('int32')):
        start = batch_size * i
        end = batch_size * (i+1)
        yield (data[start:end], start)

# Only for Keras functions with one single input tensor,
# where the output is as large as the input.
def applyFunctionBatchwise(K_function, data, batch_size=64):
    output = np.zeros_like(data, dtype=np.float32)
    for batch, start in batchGenerator(data, batch_size=batch_size):
        size = batch.shape[0]
        output[start:start+size] = K_function([batch])[0]
    return output

##################################################################
#############       Section for Custom Layers        #############
##################################################################
"""
Your custom layer can contain a number of standard keras layers
or inherit its attributes from keras.layers.
"""

def ResidualLayer(x, name, num_filters, filter_size=(5,1),
                  stride=(1,1), activation='relu',batch_norm=True):
    
    input_layer = x

    # Get the number of input feature maps.
    if K.image_data_format() == 'channels_first':
        num_input_filters = K.shape(input_layer)[1]
    else:
        num_input_filters = K.shape(input_layer)[-1]

    # The convolutional path
    layer = Activation(activation, name=name+'_act1')(input_layer)
    
    layer = Conv2D(filters=num_filters, kernel_size=filter_size,
        strides=(1,1), activation='linear', padding='same',
        use_bias=False, kernel_initializer='he_uniform',
        name=name+'_conv1')(layer)
    if batch_norm:
        layer = BatchNormalization(name=name+'_bn1', fused=False)(layer)

    layer = Activation(activation, name=name+'__act2')(layer)

    layer = Conv2D(filters=num_filters, kernel_size=filter_size,
        strides=stride, activation='linear', padding='same',
        use_bias=False, kernel_initializer='he_uniform',
        name=name+'_conv2')(layer)
    if batch_norm:
        layer = BatchNormalization(name=name+'_bn2', fused=False)(layer)


    # The skip path
    if (num_filters != num_input_filters) or (stride!=(1,1)):
        input_layer = Conv2D(filters=num_filters, kernel_size=(1,1),
        strides=stride, activation='linear', padding='same',
        use_bias=False, kernel_initializer='he_uniform',
        name=name+'_conv_skip')(input_layer)

    # Addition of the convolutional and the skip part
    layer = Add(name=name+'_add')([input_layer,layer])
    
    return layer


######### Add a reference to your custom layer in this dict #########
def get_custom_layer_dict():
    
    custom_layer_dict = {'ResidualLayer' : ResidualLayer}

    return custom_layer_dict

##################################################################
#############       Section for Custom Losses        #############
##################################################################
"""
Your custom loss must comply to the keras loss standard. The only
arguments the calculation function can have are the keras tensors
y_true, y_pred. See https://keras.io/losses/ for details.

To use more arguments, you can use a wrapper function.
"""

def margin_loss(y_true, y_pred):
    """
    Implements the margin loss used by Sabour et al. in 
    "Dynamic Routing Between Capsules" 
    (https://arxiv.org/abs/1710.09829)
    
    Arguments
    ----
        y_true:     Keras Tensor
            Contains the label informations    
            y_true.shape = [batchsize, n_classes]
        
        y_pred:     Keras Tensor
            Prediction of the model
            y_pred.shape = [batchsize, n_classes]
    
    Return
    ----
        A scalar loss value.
    """
    L = y_true * K.square(K.maximum(0., 0.9 - y_pred)) + \
        0.5 * (1 - y_true) * K.square(K.maximum(0., y_pred - 0.1))

    return K.mean(K.sum(L, 1))


######### Add a reference to your custom loss in this dict #########
def get_custom_loss_dict():
    
    custom_layer_dict = {'margin_loss' : margin_loss}

    return custom_layer_dict