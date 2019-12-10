"""
Implementation of a residual network for mass spectra.
"""

import lasagne

InputLayer = lasagne.layers.InputLayer
Conv2DLayer = lasagne.layers.Conv2DLayer
Pool2DLayer = lasagne.layers.Pool2DLayer
DenseLayer = lasagne.layers.DenseLayer
DropoutLayer = lasagne.layers.DropoutLayer
BatchNormLayer = lasagne.layers.BatchNormLayer
NonlinearityLayer = lasagne.layers.NonlinearityLayer
rectify = lasagne.nonlinearities.rectify
linear = lasagne.nonlinearities.linear
softmax = lasagne.nonlinearities.softmax    
batch_norm = lasagne.layers.batch_norm
ElemwiseSumLayer = lasagne.layers.ElemwiseSumLayer
He = lasagne.init.HeNormal
GlobalPoolLayer = lasagne.layers.GlobalPoolLayer

def ResidualLayer(inputLayer,num_filters,filter_size=(5,1),stride=(1,1),nonlinearity=rectify):
    layer = inputLayer
    layer = NonlinearityLayer(layer,nonlinearity=nonlinearity)
    layer = Conv2DLayer(layer,num_filters,filter_size,stride=(1,1),
                        nonlinearity=linear,pad='same',W=He('relu'),b=None)
    layer = BatchNormLayer(layer)
    layer = NonlinearityLayer(layer,nonlinearity=nonlinearity)
    layer = BatchNormLayer(layer)
    layer = Conv2DLayer(layer,num_filters,filter_size,stride=stride,
                        nonlinearity=linear,pad='same',W=He('relu'),b=None)
    if (num_filters != inputLayer.output_shape[1]) or (stride!=(1,1)):
        inputLayer = Conv2DLayer(inputLayer,num_filters,filter_size=(1,1),
                                 stride=stride,nonlinearity=linear,pad='same',
                                 W=He('relu'),b=None)
    layer = ElemwiseSumLayer([inputLayer,layer])
    
    return layer
           

def buildNet(inputShape, numUniqueClasses):
    """
    The buildNet function returns a lasagne.layers.Layer object, which
    contains the computational graph of the neural network.
    
    :param inputShape: The shape of the network input. The first element is
        always set to None to allow for arbitrary batch sizes.
    :param numUniqueClasses: Number of different classes.
    """
    numFeatureMaps = 16
    maxNumFeatureMaps = 128
    filterWidth = 5
    filterSize = (filterWidth,1)
    poolWidth = 3
    poolSize = (poolWidth,1)
    dataWidth = inputShape[2]

    layers = InputLayer((None,) + inputShape[1:4])
    
    while dataWidth > 30:
        
        # 
        layers = ResidualLayer(layers, numFeatureMaps,
                               filter_size = filterSize)
                               
        if numFeatureMaps < maxNumFeatureMaps:
            numFeatureMaps *= 2
            
        layers = ResidualLayer(layers, numFeatureMaps,
                               filter_size = filterSize, stride = poolSize)
        dataWidth /= poolWidth
        

    
    
    layers = ResidualLayer(layers,256,
                               filter_size = filterSize, stride = poolSize)
    layers = GlobalPoolLayer(layers)
    layers = DenseLayer(layers,num_units=numUniqueClasses,
                                   nonlinearity=softmax,W=He())
    
    return layers