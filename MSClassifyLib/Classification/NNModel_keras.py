#!/usr/bin/env python
'''
This file implements a custom class for mass spectrometry classification
with Keras. This contains the learned model and methods for
prediction as well as sensitivity analysis.
'''

import numpy as np
from tensorflow.keras import backend as K
from NNHelpers_keras import applyFunctionBatchwise

class NNModel_keras:
    """
    This class contains the learned model and methods for prediction as well
    as sensitivity analysis.
    
    :param network: keras.Model instance.
    :func classify: Classify a given dataset.
    :func scores: Output the predicted class probabilities.
    :func logits: Returns a numpy array of class logits.
    """
    def __init__(self, network):
        self.network = network


    def classify(self, data):
        """
        Returns a numpy array of class labels.
        
        :param data: Data that is to be labeled.
        """
        return np.argmax(self.network.predict(data, batch_size=32, verbose=1), axis=-1)


    def scores(self, data):
        """
        Returns a numpy array of class probabilities.
        
        :param data: Data of which the class probabilities are to be inferred.
        """
        return self.network.predict(data, batch_size=32, verbose=1)
        

    def logits(self, data):
        """
        Returns a numpy array of class logits.
        
        :param data: Data of which the class logits are to be inferred.
        """       
        logitLayer = self.network.get_layer(index=-2)
        logit_output = K.function([self.network.input], [logitLayer.output[0]])

        num_data = data.shape[0]
        logit_out = list()

        for i in range(num_data):
            logit_out.append(logit_output([np.expand_dims(data[i], axis=0)]))

        return np.array(logit_out)


    def save(self, path):
        """
        Save the Keras model to a .h5 file.
        This includes parameters, configuration and architecture.

        :param path: Saving path WITHOUT the .h5 ending
        """
        self.network.save(path+'_model.h5')


    def sensitivtyNLL(self, data, labels, batchSize=128):
        """
        Returns the sensitivity of the categorical crossentropy with respect
        to the input data.
        
        :param data: Input data.
        :param labels: Respective labels.
        :param batchSize: The network iterates through the dataset with 
            batches, whose batch size is given by this parameter.
        """

        # Gradient of the categorical crossentropy of the model output regarding
        # the model input
        grads = K.gradients(K.categorical_crossentropy(target=labels,
                                output=self.network.output, from_logits=False), 
                            self.network.input)[0]

        # Define a Keras function to calculate the gradient for a given input
        grads_func = K.function([self.network.input],
                                [grads])

        # Run the calculation for the gradients
        sens = applyFunctionBatchwise(grads_func, data, batch_size=batchSize)
        #sens = grads_func([data])

        return sens


    def sensitivityBinaryCrossentropy(self, data, labels, batchSize=128):
        """
        Returns the sensitivity of the categorical crossentropy with respect
        to the input data.
        
        :param data: Input data.
        :param labels: Respective labels.
        :param batchSize: The network iterates through the dataset with 
            batches, whose batch size is given by this parameter.
        """

        # Gradient of the binary crossentropy of the model output regarding
        # the model input
        grads = K.gradients(K.binary_crossentropy(target=labels,
                                output=self.network.output, from_logits=False), 
                            self.network.input)[0]

        # Define a Keras function to calculate the gradient for a given input
        grads_func = K.function([self.network.input], [grads])

        # Run the calculation for the gradients
        
        sens = grads_func([data])

        return sens
        

    def sensitivityPerOutput(self, data, layerIndex=-1, batchSize=128):
        """
        Returns a list of sensitivities of the output layer's activations with
        respect to the input data. This layer respresents the class
        probabilities.
        
        :param data: Input data.
        :param batchSize: The network iterates through the dataset with 
            batches, whose batch size is given by this parameter.
        """
        # Get the desired layer and save its output shape
        output_layer = self.network.get_layer(index=layerIndex)
        output_shape = output_layer.output_shape
        
        # The gradient of each output will be stored in a list
        sens = list()

        # To calculate the gradients per output, the layer output shape
        # should be (batch, neurons). Otherwise we will flatten it
        if len(output_shape) > 2:
            output_layer = K.batch_flatten(output_layer)

        print(output_shape)
        # Calculate the gradient for every output
        for i in range(output_shape[1]):
            grads = K.gradients(output_layer.output[:,i], self.network.input)[0]

            # Define a Keras function to calculate the gradient for a given input
            grads_func = K.function([self.network.input], [grads])

            # Run the calculation for the gradients
            sens.append(applyFunctionBatchwise(grads_func, data, batch_size=batchSize))

        return sens


##########################################################################################
########################      Work in Progress        ####################################
##########################################################################################
        
    def sensitivityWrtFeaturesBinary(self, data, wrtLayerNum, logitLayer=True, classes=None, argMax=1, batchSize = 128):
        """
        Returns the sensitivity of the categorical crossentropy with respect
        to the features of the specified layer. The correct label is estimated 
        by the neural network, so that this method also works for unlabeled 
        data.
        
        :param data: Input data.
        :param wrtLayer: specifies the target layer of the sensitivity
        :param batchSize: The network iterates through the dataset with 
            batches, whose batch size is given by this parameter.
        """
        Warning('The function @sensitivityWrtFeaturesBinary() is not yet implemented for the Keras backend!')
        sens = 0
        return sens

    def reconstructInput(self, samples, 
                         desiredClass, 
                         basis=None,
                         weighting=0.5,
                         argMax=1,
                         complement=False,
                         stepSize=0.1, 
                         fromLayer=1, 
                         numIterations=2,
                         coeffFeatSelection=1,
                         gradCoeff=0.01,
                         predCoeff=0.01, 
                         regCoeffL2=0.01, 
                         regCoeffL1=0.01, 
                         positiveProj=True, 
                         verbose=1000, 
                         orig=True, 
                         maldi=True,
                         fista=True,
                         smoothFista=False,
                         smooth=False,
                         logging=1000,
                         stepSizeControl=0,
                         bfgs=False):
        """
        Reconstructs input x from features of the network with FISTA

        :param samples: input samples to reconstruct the features
        :param desiredClass: specifies the desired prediction of the constructed
                             input
        :param weighting: decides whether weighting by sensitivity wrt to 
                          features is done
        :param stepSize: step size of gradient descent method
        :param fromLayer: specifies the index of the layer from which 
                          reconstruction is done
        :param numIterations: number of iterations of gradient descent
        :param predCoeff: coefficient for prediction weighting
        :param regCoeffL2: coefficient for L2-penalty (not used for reconstrion in basis)
        :param regCoeffL1: coefficient for L1-penalty
        :param positiveProj: projection onto positive values if true
        :param verbose: logging each verbose epoch 
        """
        Warning('The function @reconstructInput() is not yet implemented for the Keras backend!')
        xRecon = 0
        return xRecon   

    def reconstructInputinBasis(self, samples, basis, weighting=True,
                                stepSize=0.1, fromLayer=1, 
                                numIterations=2, regCoeffL1=0.01, 
                                positiveProj=True, verbose=True):
        """
        Reconstructs input x from features of the network with FISTA

        :param samples: input samples to reconstruct the features
        :param basis: spanning set for reconstruction
        :param stepSize: step size of gradient descent method
        :param fromLayer: specifies the index of the layer from which 
                          reconstruction is done
        :param numIterations: number of iterations of gradient descent
        :param regCoeffL1: coefficient for L1-penalty
        :param positiveProj: projection onto positive values if true
        """
        Warning('The function @reconstructInputinBasis() is not yet implemented for the Keras backend!')
        xRecon = 0
        return xRecon                    

    def computeJacobian(self, sample, layerNum, divideFeat=1):
        """
        Computation of jacobian of output wrt to input. 
        This method is computes the jacobian separately wrt to each
        feature map to reduce memory consumption. 
        Additionally, the computation for each feature map can be divided.   
        Note: The jacobian varies for due to nonlinearity of the network.
        
        param: sample jacobian for given datapoint
        param: layerNum index of output layer
        param: divideFeat parts used for dividing each output feature map
        """
        Warning('The function @computeJacobian() is not yet implemented for the Keras backend!')
        jacAll = 0
        return jacAll
