"""
Standard architecture for the NNClassifier_keras. Is is a 1-dimensional equivalent
of the 'LeNet' architecture by Yann LeCun. This code can serve as a starting
point for new architectures.
"""


from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Conv2D, Dropout, BatchNormalization, MaxPooling2D, Flatten, Activation


def buildNet(inputShape, numUniqueClasses):
    """
    The buildNet function returns a keras model object, which
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
    dataWidth = inputShape[1]

    model_input = Input(shape=inputShape, name='model_input')

    counter = 0
    
    while dataWidth > 30:
        
        # Build the connection to the input layer in the first run 
        if counter == 0:
            layers = model_input

        # Perform a 'valid' convolution, which decreases the size of the
        # data by floor(filterWidth/2). Apply batch normalization.
        layers = Conv2D(filters=numFeatureMaps, kernel_size=filterSize,
            padding='valid', activation='relu', kernel_initializer='he_uniform', 
            name = 'conv'+str(counter+1))(layers)

        layers = BatchNormalization(name='batch_norm'+str(counter+1))(layers)
        
        dataWidth -= filterWidth/2
        
        # Perform a max-pooling, which decreases the size of the data by a
        # factor of 3.
        layers = MaxPooling2D(pool_size=poolSize, name='max_pool'+str(counter+1))(layers)

        dataWidth /= poolWidth
        
        if numFeatureMaps < maxNumFeatureMaps:
            numFeatureMaps *= 2

        counter += 1
    

    if counter == 0:
        layers = model_input

    # Flatten the output of the convolutional layers
    layers = Flatten(name ='flatten1')(layers)

    # Two dense layers (with dropout) followed by a softmax layer.
    layers = Dense(units=256, activation='relu', 
        kernel_initializer='he_uniform', name = 'dense1')(layers)

    layers = BatchNormalization(name='batch_norm'+str(counter+1))(layers)
    
    layers = Dropout(rate=0.5, name='dropout1')(layers)

    layers = Dense(units=256, activation='relu', 
        kernel_initializer='he_uniform',name = 'dense2')(layers)

    layers = BatchNormalization(name='batch_norm'+str(counter+2))(layers)
    
    layers = Dropout(rate=0.5, name='dropout2')(layers)
    
    layers = Dense(units=numUniqueClasses, activation='linear', 
        name='model_logits')(layers)

    model_output = Activation(activation = 'softmax', name='model_output')(layers)

    model = Model(inputs=model_input, outputs=model_output)
    
    return model