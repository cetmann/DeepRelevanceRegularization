"""
Implementation of a IsotopeNet for mass spectra.
"""

import lasagne

InputLayer = lasagne.layers.InputLayer
Conv2DLayer = lasagne.layers.Conv2DLayer
Pool2DLayer = lasagne.layers.Pool2DLayer
DenseLayer = lasagne.layers.DenseLayer
DropoutLayer = lasagne.layers.DropoutLayer
BatchNormLayer = lasagne.layers.BatchNormLayer
NonlinearityLayer = lasagne.layers.NonlinearityLayer
linear = lasagne.nonlinearities.linear
softmax = lasagne.nonlinearities.softmax    
batch_norm = lasagne.layers.batch_norm
ElemwiseSumLayer = lasagne.layers.ElemwiseSumLayer
He = lasagne.init.HeNormal
GlobalPoolLayer = lasagne.layers.GlobalPoolLayer
           

import theano.tensor as T

from lasagne import init
from lasagne import nonlinearities
import NNHelpers

rectify = NNHelpers.relu;


__all__ = [
    "LocallyConnected2DLayer",
]


def ResidualLayer(inputLayer,num_filters,filter_size=(5,1),
                  stride=(1,1),nonlinearity=rectify):
    layer = inputLayer
    layer = NonlinearityLayer(layer,nonlinearity=nonlinearity)
    layer = Conv2DLayer(layer,num_filters,filter_size,stride=(1,1),
                        nonlinearity=linear,pad='same',W=He('relu'),b=None)
    layer = BatchNormLayer(layer)
    layer = NonlinearityLayer(layer,nonlinearity=nonlinearity)
    layer = Conv2DLayer(layer,num_filters,filter_size,stride=stride,
                        nonlinearity=linear,pad='same',W=He('relu'),b=None)
    layer = BatchNormLayer(layer)              
    if (num_filters != inputLayer.output_shape[1]) or (stride!=(1,1)):
        inputLayer = Conv2DLayer(inputLayer,num_filters,filter_size=(1,1),
                                 stride=stride,nonlinearity=linear,
                                 pad='same',W=He('relu'),b=None)
    layer = ElemwiseSumLayer([inputLayer,layer])
    
    return layer


nonlinearity = rectify;

def buildNet(inputShape, numUniqueClasses):
    """
    The buildNet function returns a lasagne.layers.Layer object, which
    contains the computational graph of the neural network.
    
    :param inputShape: The shape of the network input. The first element is
        always set to None to allow for arbitrary batch sizes.
    :param numUniqueClasses: Number of different classes.
    """
    layers = InputLayer((None,) + inputShape[1:4])
    layers = ResidualLayer(layers, 8, 
                               filter_size = (3,1))
    layers = ResidualLayer(layers, 8, 
                               filter_size = (3,1), stride= (5,1))
    layers = ResidualLayer(layers, 8, 
                               filter_size = (3,1))
    layers = ResidualLayer(layers, 1, 
                               filter_size = (3,1), stride= (3,1))
    layers = NonlinearityLayer(layers, nonlinearity = nonlinearity)
    layers = DropoutLayer(layers,p=.3)   
    layers = batch_norm(NNHelpers.LocallyConnected2DLayer(layers,1,(5,1),
                        W=He('relu'),
                        nonlinearity=nonlinearity))            
    layers = DenseLayer(layers,num_units=numUniqueClasses,
                                   nonlinearity=linear)    
    layers = NonlinearityLayer(layers, nonlinearity=softmax)                                   
    return layers