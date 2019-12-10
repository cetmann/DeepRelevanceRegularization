#!/usr/bin/env python
'''
This file implements a custom class for mass spectrometry classification
with Theano/Lasagne/Nolearn. This contains the learned model and methods for
prediction as well as sensitivity analysis.
'''

from lasagne.layers import get_output, get_output_shape
from lasagne.objectives import categorical_crossentropy
import numpy as np
import theano 
import theano.tensor as T
from nolearn.lasagne import BatchIterator
from scipy import sparse

class NNModel_theano:
    """
    This class contains the learned model and methods for prediction as well
    as sensitivity analysis.
    
    :param network: nolearn.NeuralNet instance.
    :func classify: Classify a given dataset.
    :func scores: Output the predicted class probabilities.
    :func sensitivityBinaryCrossentropy: Output the estimated NLL sensitivity
        (unsupervised).
    :func sensitivityNLL: Output the NLL sensitivity (supervised).
    :func sensitivityPerSoftmax: Output the sensitivities of the last layer.
    """
    def __init__(self, network):
        self.network = network

    def classify(self, data):
        """
        Returns a numpy array of class labels.
        
        :param data: Data that is to be labeled.
        """
        return self.network.predict(data)

    def scores(self, data):
        """
        Returns a numpy array of class probabilities.
        
        :param data: Data of which the class probabilities are to be inferred.
        """
        return self.network.predict_proba(data)
        
    def logits(self, data):
        """
        Returns a numpy array of class logits.
        
        :param data: Data of which the class probabilities are to be inferred.
        """       
        logitLayer = self.network.layers_[-2]
        return self.network.get_output(logitLayer,data)

    def sensitivtyNLL(self, data, labels,  batchSize = 128):
        """
        Returns the sensitivity of the categorical crossentropy with respect
        to the input data.
        
        :param data: Input data.
        :param labels: Respective labels.
        :param batchSize: The network iterates through the dataset with 
            batches, whose batch size is given by this parameter.
        """
        sens = np.zeros(data.shape)
        labelMatrix = T.ivector('labelVector')
        # Compute number of batches
        numBatches = int(np.ceil(float(sens.shape[0]) / float(batchSize)))
        startBatch = 0
        inputLayer = self.network.layers_[0].input_var
        output = get_output(self.network.layers_[-1], deterministic=True)
        score = categorical_crossentropy(output, labelMatrix).sum()
        calculatedGradients = theano.grad(score,inputLayer)
        for i in range(numBatches):
            endBatch = startBatch + batchSize
            if endBatch >= sens.shape[0]:
                endBatch = sens.shape[0]
                batchSize = endBatch - startBatch
            inputData = data[startBatch:endBatch].reshape(batchSize, 
                                data.shape[1], data.shape[2], data.shape[3])
            inputLabels = labels[startBatch : endBatch]
            sens[startBatch:endBatch] = \
                calculatedGradients.eval({inputLayer: inputData, 
                                          labelMatrix: inputLabels.astype('int32')})
            startBatch = endBatch        
        return sens



    def sensitivityBinaryCrossentropy(self, data,  batchSize = 128):
        """
        Returns the sensitivity of the categorical crossentropy with respect
        to the input data.
        
        :param data: Input data.
        :param labels: Respective labels.
        :param batchSize: The network iterates through the dataset with 
            batches, whose batch size is given by this parameter.
        """
        sens = np.zeros(data.shape)
        labelMatrix = T.ivector('labelVector')
        # Compute number of batches
        numBatches = int(np.ceil(float(sens.shape[0]) / float(batchSize)))
        startBatch = 0
        inputLayer = self.network.layers_[0].input_var
        output = get_output(self.network.layers_[-1], deterministic=True)
        score = categorical_crossentropy(output, labelMatrix).sum()
        calculatedGradients = theano.grad(score,inputLayer)
        for i in range(numBatches):
            endBatch = startBatch + batchSize
            if endBatch >= sens.shape[0]:
                endBatch = sens.shape[0]
                batchSize = endBatch - startBatch
            inputData = data[startBatch:endBatch].reshape(batchSize, 
                                data.shape[1], data.shape[2], data.shape[3])
            pred = output.eval({inputLayer: inputData}).argmax(axis=1)
            sens[startBatch:endBatch] = \
                calculatedGradients.eval({inputLayer: inputData, 
                                          labelMatrix: pred.astype('int32')})
            startBatch = endBatch        
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
        inputLayer = self.network.layers_[0].input_var
        output = get_output(self.network.layers_[layerIndex], deterministic=True)
        

        # Compute number of classes (using a batch to circumvent memory issues)
        if batchSize >= data.shape[0]:
            endBatch = data.shape[0]
        else:
            endBatch = batchSize
            
        inputData = data[0:endBatch].reshape(endBatch,
                                data.shape[1], data.shape[2], data.shape[3])
        numClasses = output.eval({inputLayer: inputData}).shape[1]
        
        # Compute number of batches
        numBatches = int(np.ceil(float(data.shape[0]) / float(batchSize)))
        
        
        classIndex = T.iscalar('classIndex')
        calculatedGradients =  theano.grad(output[:,classIndex].sum(),
                                           inputLayer)  
        
        # List of sensitivity per class probability.
        sens = [];
        # save original batch size
        origBatchSize = batchSize
        for j in range(numClasses):
            currentGrad = np.zeros(data.shape)
            startBatch = 0
            for i in range(numBatches):
                endBatch = startBatch + batchSize
                if endBatch >= data.shape[0]:
                    endBatch = data.shape[0]
                    batchSize = endBatch - startBatch
                inputData = data[startBatch:endBatch].reshape(batchSize,
                                data.shape[1], data.shape[2], data.shape[3])
                currentGrad[startBatch:endBatch] =\
                    calculatedGradients.eval({inputLayer: inputData, 
                                              classIndex: j})
                startBatch = endBatch 
                if np.sum(currentGrad[0]) == 0:
                    print(i)
            sens.append(currentGrad)
            # reset batch size to original batch size
            batchSize = origBatchSize
        return sens

        
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
        layerShape = get_output_shape(self.network.layers_[wrtLayerNum])
        sensShape = (data.shape[0], layerShape[1], layerShape[2], layerShape[3])
        sens = np.zeros(sensShape)
        
        # Compute number of batches
        numBatches = int(np.ceil(float(sens.shape[0]) / float(batchSize)))
        
        startBatch = 0
        for i in range(numBatches):
            endBatch = startBatch + batchSize
            if endBatch >= sens.shape[0]:
                endBatch = sens.shape[0]
                batchSize = endBatch - startBatch
   
            inputData = data[startBatch:endBatch].reshape(batchSize, 
                                data.shape[1], data.shape[2], data.shape[3])
            inputLayer = self.network.layers_[0].input_var

            if logitLayer: # gradient from logits
                logits, activation = get_output([self.network.layers_[-2], 
                                                 self.network.layers_[wrtLayerNum]],
                                                deterministic=True)
            else: # gradient from softmax-activations
                logits, activation = get_output([self.network.layers_[-1], 
                                                 self.network.layers_[wrtLayerNum]],
                                                deterministic=True)
                                                
            if classes is None:
                # take derivative from larrgest logit
                if -argMax+1 == 0:
                    logitMax = logits.eval({inputLayer: inputData}).argsort(axis=1)[:,-argMax:]
                else:
                    logitMax = logits.eval({inputLayer: inputData}).argsort(axis=1)[:,-argMax:-argMax+1]        
                sumLogitsMax = logits[np.arange(data.shape[0]), logitMax].sum()
            else:
                # take derivative from pre-defined logit
                sumLogitsMax = logits[np.arange(data.shape[0]), classes].sum()                                 
            sens[startBatch:endBatch]=  theano.grad(sumLogitsMax,
                activation).eval({inputLayer: inputData})
            startBatch = endBatch
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
        # select input layer to get features  
        inputLayer = self.network.layers_[0].input_var
        # get output features from selected layer
        output = get_output(self.network.layers_[fromLayer], deterministic=True)
        features = output.eval({inputLayer: samples})
        # selection of feature maps
        if coeffFeatSelection < 1:
            mask = np.ones((features.shape[0], features.shape[1]))
            sens = self.sensitivityWrtFeaturesBinary(samples, wrtLayerNum=fromLayer, logitLayer=False)
            product = sens * features
            for i in range(samples.shape[0]):    
                #normFeat = np.linalg.norm(features[i], ord=2, axis=(1,2))
                normFeat = np.mean(np.abs(product[i]), axis=(1,2))
                argInd = np.argsort(normFeat)[::-1]
                numFeat = np.floor(coeffFeatSelection * features.shape[1]).astype('int32')
                zeroIndices = argInd[numFeat:-1]
                mask[i, zeroIndices] = 0
            #mask = mask.astype('bool')
            output = output[mask.nonzero()]
            features = features[mask.nonzero()]
        
        # create initial input (uniform noise in [0,1])
        shapeInput = samples.shape
        # use original input
        if orig and not type(basis)==np.ndarray:
            xRecon = samples  +  0.000000001 * np.random.randn(shapeInput[0], 
                shapeInput[1], shapeInput[2], shapeInput[3]).astype('float32') 
        elif not orig and not type(basis)==np.ndarray:
            # different starting values for vision and maldi data
            if maldi:
                xRecon = 1e-4 *np.random.uniform(size=(shapeInput[0], shapeInput[1], shapeInput[2], 
                                    shapeInput[3])).astype('float32')
            else:
                xRecon = 300*np.random.uniform(size=(shapeInput[0], shapeInput[1], shapeInput[2], shapeInput[3])).astype('float32') - 150
        else: # reconstruction in basis
            # random init for basis coefficients
            numBasis = basis.shape[3]
            xRecon = np.random.uniform(size=(samples.shape[0], numBasis)).astype('float32')
            xRecon = xRecon / np.sum(xRecon, 1, keepdims=True)
            # define product of basis coefficients and basis
            z = theano.tensor.matrix()
            self.network.layers_[0].input_var = theano.tensor.tensordot(z, basis, [[1], [3]])
            output = get_output(self.network.layers_[fromLayer], deterministic=True)
            inputLayer = self.network.layers_[0].input_var
            # set l2-penalty to zero
            regCoeffL2 = 0
#        xRecon = np.ones((shapeInput[0], shapeInput[1], shapeInput[2], 
#                                shapeInput[3])).astype('float32')
        # normalize init with l1-norm per first dimension (TIC normalization)
        #xRecon = xRecon #/ np.sum(xRecon, (1,2,3), keepdims=True)  
              
        ### definition of terms for cost function ###
        difference = output - features
        # incorporating gradient (wrt input) difference in feature layer
        logits, activation = get_output([self.network.layers_[-2], 
                                         self.network.layers_[fromLayer]],
                                         deterministic=True)
        # take logits from argmax of orig sample
        logitMax = logits.eval({inputLayer: samples}).argmax(axis=1)          
        sumLogitsMax = logits[np.arange(samples.shape[0]), logitMax].sum()                                         
        sensRecon = theano.grad(sumLogitsMax, activation)
        sensSamples = theano.grad(sumLogitsMax, activation).eval({inputLayer: samples})    
        gradDifference = sensRecon - sensSamples
        scalingFactorGrad = np.linalg.norm(sensSamples)
        # incorporating prediction as penalty term
        predPenalty = categorical_crossentropy(
            get_output(self.network.layers_[-1], deterministic=True), desiredClass).sum()
        gradTerm = theano.tensor.sqr(theano.tensor._tensor_py_operators.norm(gradDifference,2)) / scalingFactorGrad
        ##### definition of cost function   
        # weighting by scores
        if weighting > 0 and weighting < 1:
            # weighting of discrepancy
            if smooth:
                # smoothgrad
                variance = 0.1 * (np.max(samples) - np.min(samples))  
                m = 30
                smoothGradient = np.zeros(features.shape)
                for i in range(1, m+1):
                    print(i)
                    currSamples = samples + variance * np.random.randn(samples.shape[0], 
                                                                       samples.shape[1],samples.shape[2], samples.shape[3]).astype('float32')
                    sensCurr = self.sensitivityWrtFeaturesBinary(currSamples, wrtLayerNum=fromLayer, argMax=argMax, logitLayer=False)
                    smoothGradient += 1. / m  * sensCurr
                sensFeat = np.abs(smoothGradient)
            else:
                # (use sensitivity wrt to current layer as weighting factor)
                sensFeat = np.abs(self.sensitivityWrtFeaturesBinary(samples, fromLayer, argMax=argMax, logitLayer=False))
            ########## TODO: change?
#            zerosMask = weight <= weighting * np.max(np.abs(weight))
#            weight[zerosMask] = 0
            ##########
            weight = np.ones(features.shape)
            for i in range(samples.shape[0]):
                argInd = np.argsort(sensFeat[i], axis=None)[::-1]
                numSel = np.floor(weighting * features.shape[1] * features.shape[2] * features.shape[3]).astype('int32')
                if complement: # reconstruction from "unimportant" features
                    zeroIndices = argInd[:numSel]
                else: # reconstruction from "important" features
                    zeroIndices = argInd[numSel:-1]
                currMask = np.ones(argInd.shape)
                currMask[zeroIndices] = 0
                weight[i] = currMask.reshape((features.shape[1], features.shape[2], features.shape[3]))           
            wDifference = weight * difference
            # scale discrepancy term by weighted features
            scalingFactorDisc = np.linalg.norm(weight * features)
            cost = (theano.tensor.sqr(
                        theano.tensor._tensor_py_operators.norm(wDifference, 2)) / scalingFactorDisc)
        else:      
            # scale discrepancy term by features
            scalingFactorDisc = np.linalg.norm(features)
            cost = (theano.tensor.sqr(
                        theano.tensor._tensor_py_operators.norm(difference, 2)) / scalingFactorDisc)
        # l2 penalty
        cost += regCoeffL2 * theano.tensor.sqr(theano.tensor._tensor_py_operators.norm(inputLayer, 2))
        # prediction term
        if predCoeff > 0:
            cost += predCoeff * predPenalty
        # gradient term
        if gradCoeff > 0:
            # weigth gradient term with above selection
            if weighting > 0 and weighting < 1:
                gradDifference = weight * gradDifference
            gradTerm = theano.tensor.sqr(theano.tensor._tensor_py_operators.norm(gradDifference,2)) / scalingFactorGrad
            cost += gradCoeff * gradTerm  
        # compute scaling term in soft-threshold
        softThres = regCoeffL1 * stepSize
        
        # init auxilary point y (for FISTA)
        yAux = xRecon
        # init step parameter for FISTA
        tFista = 1

        # compilation
        if not type(basis)==np.ndarray:
            # compile gradient function
            gradientX = theano.grad(cost, inputLayer)
            #### TODO: remove
            gradientGradTerm = theano.tensor.sqr(theano.tensor._tensor_py_operators.norm(theano.grad(gradTerm, inputLayer), 2))
            gradientGradTermFunction = theano.function([inputLayer], gradientGradTerm)
            #### 
            gradientFunction = theano.function([inputLayer], gradientX)
            # compile logging functions
            third =  theano.tensor.sqr(theano.tensor._tensor_py_operators.norm(inputLayer, 2))
            forth = theano.tensor._tensor_py_operators.norm(inputLayer, 1)
            fifth = categorical_crossentropy(
                            get_output(self.network.layers_[-1], deterministic=True), desiredClass).sum()
            loggingFunction = theano.function([inputLayer], [cost, third, forth, fifth, gradTerm])
        else:
            gradientX = theano.grad(cost, z)
            gradientFunction = theano.function([z], gradientX)
            # compile logging functions
            third =  theano.tensor.sqr(theano.tensor._tensor_py_operators.norm(z, 2))
            forth = theano.tensor._tensor_py_operators.norm(z, 1)
            fifth = categorical_crossentropy(
                            get_output(self.network.layers_[-1], deterministic=True), desiredClass).sum()
            loggingFunction = theano.function([z], [cost, third, forth, fifth, gradTerm])
        
        if bfgs:
            from scipy.optimize import minimize
            costTheanoFun = theano.function([inputLayer], cost)
            def costFun(inputValue):
                inputValue = inputValue.reshape(samples.shape).astype('float32')
                return costTheanoFun(inputValue).astype('float64')
            def gradPyFun(inputValue):
                inputValue = inputValue.reshape(samples.shape).astype('float32')
                gradX = gradientFunction(inputValue).astype('float64')
                return np.ndarray.flatten(gradX)

#            costFun = lambda x, costTheanoFun, gradientFunction: costTheanoFun(x)
#            gradPyFun = lambda x, costTheanoFun, gradientFunction: gradientFunction(x)
            print('bfgs')
            xRecon = minimize(costFun, xRecon, method='L-BFGS-B', jac=gradPyFun, options={'maxiter':numIterations , 'disp':True})
            return xRecon
        
        # prepare loggin
        if logging:
            logDictionary = dict()
            logDictionary['objective'] = list()
            logDictionary['discrepancy'] = list()
            logDictionary['l2-term'] = list()
            logDictionary['l1-term'] = list()
            logDictionary['pred-term'] = list()
            logDictionary['grad term'] = list()
            logDictionary['gradient norm'] = list()
        # gradient descent on initial image combined with Soft-Thresholding
        # and extension by FISTA
        for i in range(1, numIterations+1):
            if stepSizeControl and i % stepSizeControl == 0:
                stepSize *= 0.5
                softThres = regCoeffL1 * stepSize
                            
            # save old iterate for update of auxilary variable
            xReconOld = xRecon
            # compute gradient (without l1-term)
            if fista:
                # gradient step using auxilary variable
                if smoothFista and i > 700:
                    # DUMB SMOOTHING
                    noiseLevel = 0.1
                    variance = noiseLevel * (np.max(samples) - np.min(samples))  
                    m = 20
                    gradientX = np.zeros(yAux.shape).astype('float32')
                    for smoothIt in range(1, m+1):
                        currYAux = yAux + variance * np.random.randn(yAux.shape[0], 
                                                                           yAux.shape[1],yAux.shape[2], yAux.shape[3]).astype('float32')
                        gradCurr = gradientFunction(currYAux)
                        gradientX += 1. / m  * gradCurr
                else: 
                    gradientX = gradientFunction(yAux)
                xRecon = yAux - 2*stepSize*gradientX
                # soft thresholding
                xRecon = np.sign(xRecon) * np.maximum(
                    np.abs(xRecon) - softThres, 0.)
                # projection onto positive values (if positiveProj is true)
                if positiveProj:
                    xRecon[xRecon <0] = 0
                # update step parameter t for FISTA
                tFistaNew = (1 + np.sqrt(1 + 4* (tFista**2))) / 2
                # update auxilary variable
                yAux = xRecon + ((tFista - 1) / tFistaNew) * (xRecon - xReconOld)
                # set step parameter t
                tFista = tFistaNew
            else:
                gradientX = gradientFunction(xRecon)
                xRecon = xRecon - 2*stepSize*gradientX
                # soft thresholding
                xRecon = np.sign(xRecon) * np.maximum(
                    np.abs(xRecon) - softThres, 0.)
                # projection onto positive values (if positiveProj is true)
                if positiveProj:
                    xRecon[xRecon <0] = 0                
            # print objective function output and norm of gradient
            # (every verbose iterations)
            if logging > 0 and i % logging == 0 or i == 1:                                   
                cost, third, forth, fifth, sixth = loggingFunction(xRecon)
                eight = gradientGradTermFunction(xRecon)                             
                first = cost + regCoeffL1 * forth
                second = first - regCoeffL2 * third - regCoeffL1 * forth - predCoeff * fifth - gradCoeff * sixth              
                seventh = np.linalg.norm(gradientX)
                logDictionary['objective'].append(first)
                logDictionary['discrepancy'].append(second)
                logDictionary['l2-term'].append(third)
                logDictionary['l1-term'].append(forth)
                logDictionary['pred-term'].append(fifth)
                logDictionary['grad term'].append(sixth)
                logDictionary['gradient norm'].append(seventh)                
                if verbose > 0 and i % verbose == 0 or i == 1:
                    text = ' objective: {} , discrepancy: {}, l2-term: {}, l1-term: {}, pred-term: {}, grad term: {}, gradient norm: {}, norm grad term: {}'
                    print(i)
                    print(text.format(first, second, third, forth, fifth, sixth, seventh, eight))
        if logging:
            return xRecon, logDictionary
        else:
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

        numBasis = basis.shape[3]
        xRecon = np.random.uniform(size=(samples.shape[0], numBasis)).astype('float32')
        #xRecon = np.ones((samples.shape[0], numBasis)).astype('float32')
        xRecon = 1*(xRecon / np.sum(xRecon, 1, keepdims=True) )
        
        # select input layer to get features  
        inputLayer = self.network.layers_[0].input_var
        
        # get output features from selected layer
        output = get_output(self.network.layers_[fromLayer], deterministic=True)
        features = output.eval({inputLayer: samples})
        
        z = theano.tensor.matrix()
        self.network.layers_[0].input_var = theano.tensor.tensordot(z, basis, [[1], [3]])
        output = get_output(self.network.layers_[fromLayer], deterministic=True)
        # definition of cost function
        difference = output - features
        # check if weighting is used
        if weighting:
            # use sensitivity wrt to current layer as weighting factor
            weight = np.abs(self.sensitivityWrtFeaturesBinary(samples, fromLayer))
            wDifference = weight * difference
            # scale discrepancy term by weighted features
            scalingFactorDisc = np.linalg.norm(weight * features)
            # discrepancy term + l2-penalty
            cost = (theano.tensor.sqr(
                theano.tensor._tensor_py_operators.norm(wDifference, 2)) / scalingFactorDisc)
        else:
            # scale discrepancy term by features
            scalingFactorDisc = np.linalg.norm(features)
            # discrepancy term + l2-penalty
            cost = (theano.tensor.sqr(
                theano.tensor._tensor_py_operators.norm(difference, 2)) / scalingFactorDisc)
            
        # compute scaling term in soft-threshold
        softThres = regCoeffL1 * stepSize
        
        # init auxilary point y (for FISTA)
        yAux = xRecon
        # init step parameter for FISTA
        tFista = 1
        
        # auxilary for verbose
        if verbose:
            forth = 0        
        
        # gradient descent on initial image combined with Soft-Thresholding
        # and extension by FISTA
        for i in range(1, numIterations+1):
            print(i)
            # save old iterate for update of auxilary variable
            xReconOld = xRecon
            
            # compute gradient (without l1-term)
            gradientX = theano.grad(cost, z).eval({z: yAux})
            # gradient step using auxilary variable
            xRecon = yAux - 2*stepSize*gradientX
            # soft thresholding
            xRecon = np.sign(xRecon) * np.maximum(
                np.abs(xRecon) - softThres, 0.)
            # projection onto positive values (if positiveProj is true)
            if positiveProj:
                xRecon[xRecon <0] = 0
            # update step parameter t for FISTA
            tFistaNew = (1 + np.sqrt(1 + 4* (tFista**2))) / 2
            # update auxilary variable
            yAux = xRecon + ((tFista - 1) / tFistaNew) * (xRecon - xReconOld)
            # set step parameter t
            tFista = tFistaNew
                
            # print objective function output and norm of gradient
            if verbose:
                if regCoeffL1 != 0:
                    forth = np.linalg.norm(xRecon, ord=1)
                first = cost.eval({z: xRecon}) + regCoeffL1 * forth
                second = first -  regCoeffL1 * forth
                fifth = np.linalg.norm(gradientX) 
                text = 'objective: {} , discrepancy: {}, l1-term: {}, gradient norm: {}'
                print(text.format(first, second, forth, fifth))
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
        layers= self.network.layers_ 
        outputLayer = layers[layerNum]
        inputLayer = layers[0]
        networkInput = inputLayer.input_var
        networkOutput = get_output(outputLayer, deterministic=True)
        # get shape of output
        outShape = networkOutput.eval({networkInput:sample}).shape
        # get shape of input
        inputShape = sample.shape
        # intialize list of jacobians
        jacobianList = list()
        # iteration over channels
        for j in range(0, outShape[1]):
            # iteration over channel parts (only along first dimension)
            countPart = 0
            sizePart = outShape[2] / divideFeat
            # initialize jacobian for current channel
            outDim = outShape[2] * outShape[3]
            jacobianOutput = np.zeros((outDim, inputShape[0], 
                                   inputShape[1], inputShape[2], inputShape[3]))  
            for k in range(0, divideFeat):
                # consider rest
                if k == divideFeat -1:
                    rest = outShape[2] % divideFeat
                else:
                    rest = 0           
                # divide output channel
                currOutput = networkOutput[0, j, countPart * sizePart : 
                                           (countPart + 1) * sizePart + rest, :]
                # compute jacobian for current output part
                jacobianCurrOutput = theano.gradient.jacobian(currOutput.flatten(), networkInput)
                jacFunction = theano.function([networkInput], jacobianCurrOutput)
                # evaluate function
                jacExample = jacFunction(sample)
                # stack jacobian 
                jacobianOutput[countPart * sizePart :
                              (countPart + 1) * sizePart + rest] = jacExample
                countPart += 1
            # reshape jacobian
            imageSize = sample.shape[1] * sample.shape[2] * sample.shape[3]
            jacobianOutput = jacobianOutput.reshape((jacobianOutput.shape[0], imageSize))
            # convert jacobian to sparse matrix
            jacobianOutput = sparse.csc_matrix(jacobianOutput)
            # append jacobian of channel to list of jacobians
            jacobianList.append(jacobianOutput)
        jacAll = sparse.vstack(jacobianList)
        return jacAll
