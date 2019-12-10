"""
Standard architecture for the NNClassifier_theano. Is is a 1-dimensional equivalent
of the 'LeNet' architecture by Yann LeCun. This code can serve as a starting
point for new architectures.
"""

import lasagne

InputLayer = lasagne.layers.InputLayer
batch_norm = lasagne.layers.batch_norm
DenseLayer = lasagne.layers.DenseLayer
DropoutLayer = lasagne.layers.DropoutLayer
Conv2DLayer = lasagne.layers.Conv2DLayer
Pool2DLayer = lasagne.layers.Pool2DLayer

rectify = lasagne.nonlinearities.rectify
softmax = lasagne.nonlinearities.softmax

He = lasagne.init.HeNormal

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
        
        # Perform a 'valid' convolution, which decreases the size of the
        # data by floor(filterWidth/2). Apply batch normalization.
        layers = batch_norm(Conv2DLayer(layers, numFeatureMaps, filterSize,
                                        pad = 'valid',W = He('relu'), b=None,
                                        nonlinearity = rectify))
        dataWidth -= filterWidth/2
        
        # Perform a max-pooling, which decreases the size of the data by a
        # factor of 3.
        layers = Pool2DLayer(layers,poolSize)
        dataWidth /= poolWidth
        
        if numFeatureMaps < maxNumFeatureMaps:
            numFeatureMaps *= 2
    
    
    # Two dense layers (with dropout) followed by a softmax layer.
    layers = batch_norm(DenseLayer(layers, 256, W = He('relu'), b=None,
                                   nonlinearity = rectify))
    layers = DropoutLayer(layers)
    
    layers = batch_norm(DenseLayer(layers, 256, W = He('relu'), b=None,
                                   nonlinearity = rectify))
    layers = DropoutLayer(layers)
    
    
    layers = DenseLayer(layers, numUniqueClasses, nonlinearity = softmax)
    
    return layers