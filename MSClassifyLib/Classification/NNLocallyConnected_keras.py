"""
Keras Implementation of a IsotopeNet for mass spectra.
"""

from tensorflow.keras.models import Model
from tensorflow.keras import backend as K
from tensorflow.keras.layers import Activation, Add, BatchNormalization, Conv2D
from tensorflow.keras.layers import Dropout, Dense, Flatten, LocallyConnected2D, Input
from tensorflow.keras.layers import ZeroPadding2D

from NNHelpers_keras import ResidualLayer


def buildNet(inputShape, numUniqueClasses):
    """
    :param inputShape: The shape of the network input. The first element is
        always set to None to allow for arbitrary batch sizes.
    :param numUniqueClasses: Number of different classes.
    """

    model_input = Input(shape=inputShape, name='model_input')

    layers = ResidualLayer(model_input, name='res1', num_filters=8,
        filter_size=(3,1))

    # First spatial reduction
    layers = ResidualLayer(layers, name='res2', num_filters=8,
        filter_size=(3,1), stride=(5,1))

    layers = ResidualLayer(layers, name='res3', num_filters=8,
        filter_size=(3,1))

    # Second spatial reduction
    layers = ResidualLayer(layers, name='res4', num_filters=1,
        filter_size=(3,1), stride=(3,1))

    layers = Activation(activation='relu', name='act1')(layers)

    layers = Dropout(rate=0.3, name='dropout1')(layers)

    layers = ZeroPadding2D((2,0))(layers)

    # Note: deleted the bias and moved the order of layer-BN-NL to match
    # the original lasagne version, which was done with the batch_norm method.
    layers = LocallyConnected2D(filters=1, kernel_size=(5,1),
        kernel_initializer='he_uniform', activation='linear', use_bias=False,
        name='loc_con1')(layers)

    layers = BatchNormalization(name='bn1')(layers) 

    layers = Activation(activation = 'relu', name='act_after_loc')(layers)

    layers = Flatten(name='flatten1')(layers)

    layers = Dense(units=numUniqueClasses, activation='linear', 
        name='model_logits')(layers)
    
    model_output = Activation(activation = 'softmax', name='model_output')(layers)

    model = Model(inputs=model_input, outputs=model_output)
                               
    return model
