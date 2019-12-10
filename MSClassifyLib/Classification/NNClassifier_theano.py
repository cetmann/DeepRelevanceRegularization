#!/usr/bin/env python
'''
This file implements a custom class for mass spectrometry classification
with Theano/Lasagne/nolearn. After training, a NNModel_theano instance is created.
'''

import numpy as np

# NN-specific libraries
import theano
import theano.tensor as T
from theano.sandbox.rng_mrg import MRG_RandomStreams as RandomStreams
#from theano.sandbox.cuda.dnn import dnn_conv
from theano.tensor.shared_randomstreams import RandomStreams as rStream
#from theano.compile.nanguardmode import NanGuardMode

import lasagne
from lasagne import regularization
from lasagne.objectives import aggregate, categorical_crossentropy
from lasagne.updates import nesterov_momentum
from lasagne.random import get_rng
from lasagne.layers import get_output, InputLayer
from lasagne.utils import one_hot

import nolearn
from nolearn.lasagne import NeuralNet, BatchIterator, TrainSplit, objective

# Custom output class
from NNModel_theano import NNModel_theano

# Libraries for data processing
from sys import stdout
from collections import OrderedDict
import time

# Functions for the modifiedObjective method
oneHot = lambda pred,label : one_hot(label,np.array(pred).shape[-1])
lossFn = lambda pred,label : categorical_crossentropy(pred,oneHot(pred,label))


def categorical_crossentropy_logdomain(log_predictions, targets):
    return -T.sum(targets * log_predictions, axis=1)

# Custom objective function capable of adding total variation regularization.
def modifiedObjective(layers,
              loss_function,
              target,
              aggregate=aggregate,
              deterministic=False,
              l1=0,
              l2=0,
              logitSens=0,
              probSens=0,
              lossSens=0,
              std=None,
              get_output_kw=None):
    """
    Modified implementation of the NeuralNet objective.

    :param layers: The underlying layers of the NeuralNetwork
    :param loss_function: The callable loss function to use
    :param target: the expected output
    :param aggregate: the aggregation function to use
    :param deterministic: Whether or not to get a deterministic output
    :param l1: Optional l1 regularization parameter
    :param l2: Optional l2 regularization parameter
    :param lossSens: Optional loss sensitivity regularization parameter
    :param lossSens: Optional loss sensitivity regularization parameter
    :param lossSens: Optional loss sensitivity regularization parameter
    :param get_output_kw: optional kwargs to pass to
                          :meth:`NeuralNetwork.get_output`
    :return: The total calculated loss
    """
    if get_output_kw is None:
        get_output_kw = {}
    output_layer = layers[-1]
    logit_layer = layers[-2]
    input_layer = layers[0]
    network_input = input_layer.input_var
    network_output = get_output(output_layer, deterministic=deterministic,
                                **get_output_kw)
    logit_output = get_output(logit_layer, deterministic=deterministic,
                                **get_output_kw)

    L = loss_function(network_output,
                      lasagne.utils.one_hot(target,
                                            output_layer.output_shape[1]))
    loss = aggregate(L)

    if l1:
        loss += regularization.regularize_layer_params(
            layers.values(), regularization.l1) * l1
            
    if l2:
        loss += regularization.regularize_layer_params(
            layers.values(), regularization.l2) * l2
    
    # logit sensitivity        
    if logitSens:
        logit = T.sum(logit_output * lasagne.utils.one_hot(target,
                                            output_layer.output_shape[1]),
                                            axis=1)
        G_logit = T.grad(T.sum(logit), network_input)

        if std is not None:
            G_logit = std * G_logit
        
        # Sparse saliency regularization
        absG_logit = T.abs_(G_logit)
        sumAbsG_logit = T.sum(absG_logit,axis = (1,2,3))
        loss += aggregate(sumAbsG_logit) * logitSens
    
    # probability sensitivity    
    if probSens:
        prob = T.sum(network_output * lasagne.utils.one_hot(target,
                                            output_layer.output_shape[1]),
                                            axis=1)
        G_prob = T.grad(T.sum(prob), network_input)

        if std is not None:
            G_prob = std * G_prob
        
        # Sparse saliency regularization
        absG_prob = T.abs_(G_prob)
        sumAbsG_prob = T.sum(absG_prob,axis = (1,2,3))
        loss += aggregate(sumAbsG_prob) * probSens        
     
    # Loss sensitivity
    if lossSens:        
        G_loss = theano.grad(T.sum(L),network_input)
        if std is not None:
            G_loss = std * G_loss
        absG_loss = T.abs_(G_loss)
        loss += aggregate(T.sum(absG_loss,axis = (1,2,3))) * lossSens


        # Double Backpropagation, uncomment if desired
        #sqG = G**2
        #sumSqG = T.sum(sqG,axis = (1,2,3))
        #loss += aggregate(sumSqG) * tv
    return loss

class NNClassifier_theano:
    """
    NNClassifier_theano is a custom class for mass spectra classification. It is
    essentially a wrapper for nolearn and processes the hyperparameters
    given by the neuralNetInterface.

    :param architecture: lasagne.layers.Layer object
    :param hyperparameter: dictionary of hyperparameter-value-pairs
    """
    def __init__(self, architecture, hyperparameter = {}):
        self.archi = architecture
        self.hyperp = hyperparameter
        self._srng = RandomStreams(get_rng().randint(1, 2147462579)) # for adaptive noise
        self._srng2 = rStream(2147462579)

        # Create nolearn ModifiedNeuralNet object
        self.classifier  = ModifiedNeuralNet(
            layers=self.archi,
            max_epochs=self.hyperp.setdefault('epochs',100),
            update=self.hyperp.setdefault('optimizer',lasagne.updates.adam),
            update_learning_rate=self.hyperp.setdefault('learningRate',0.001),
            objective = modifiedObjective,
            objective_logitSens = self.hyperp.setdefault('logitSens',0.),
            objective_probSens = self.hyperp.setdefault('probSens',0.),
            objective_lossSens = self.hyperp.setdefault('lossSens',0.),
            objective_std = self.hyperp.setdefault('trainingDataStd',None),
            objective_loss_function=categorical_crossentropy,
            verbose=0,
            batch_iterator_train = DataAugmentationBatchIterator(
                self.hyperp.setdefault('batchSize',64),
                disturbLabelRate=self.hyperp.setdefault('disturbLabelRate',0),
                sdWidth=self.hyperp.setdefault('sdWidth',0),
                sdNumber=self.hyperp.setdefault('sdNumber',0),
                shuffle=True),
            batch_iterator_test = nolearn.lasagne.BatchIterator(
                self.hyperp.setdefault('batchSize',64),shuffle=False),\
            train_split = TrainSplit(eval_size=self.hyperp.setdefault(
                'validationSetRatio',.1)),
            objective_l1 = self.hyperp.setdefault('l1',0.),
            objective_l2 = self.hyperp.setdefault('l2',0.01),
            on_training_started=[nolearn.lasagne.PrintLayerInfo()],
            on_epoch_finished=[getIndividualLosses,
                               printError,
                               addEndTimeToHistory,
                               printAdaptiveNoise,
                               saveBestValidNet])
        self.classifier.initialize()

    def trainModel(self, data, classes):
        validationSetRatio = self.hyperp.setdefault('validationSetRatio',.1)
        if validationSetRatio != 0:
            initialShuffleIndices = np.arange(len(data),dtype='int32')
            np.random.shuffle(initialShuffleIndices)
            data = data[initialShuffleIndices]
            classes = classes[initialShuffleIndices]
            
            numOfValPoints = np.floor(validationSetRatio * len(data))
            validationData = data[:np.int32(numOfValPoints)]
            validationClasses = classes[:np.int32(numOfValPoints)]
            
            data = data[np.int32(numOfValPoints):]
            classes = classes[np.int32(numOfValPoints):]
        
        y = T.ivector()
        
        verbose = self.hyperp.setdefault('verbose',1.)        
        
        maxEpochs = self.hyperp.setdefault('epochs',100)
        optimizer = self.hyperp.setdefault('optimizer',lasagne.updates.adam)
        learningRate=self.hyperp.setdefault('learningRate',0.001)
        
        if self.hyperp['useSensRegControl']:
            regMultiplier = theano.shared(0.)
        else:
            regMultiplier = theano.shared(1.)    
        
        # RENAME THIS GLOBALLY
        std = self.hyperp.setdefault('trainingDataStd',None)
        
        biTrain = nolearn.lasagne.BatchIterator(
                self.hyperp.setdefault('batchSize',64),shuffle=True)
        biVal = nolearn.lasagne.BatchIterator(
                self.hyperp.setdefault('batchSize',64),shuffle=False)
        #biTest = nolearn.lasagne.BatchIterator(self.hyperp.setdefault('batchSize',64),shuffle=False)
       
        layers = self.classifier.layers_
        outputLayer = layers[-1]
        inputLayer = layers[0]
        networkInput = inputLayer.input_var
      
        ###### FOR DEBUGGING, DELETE LATER ########
        print(self.hyperp)       
       
        # computation of losses for training (not deterministic)
        classificationLoss, totalLoss, l1Loss, l2Loss, logitSensLoss, logitDiffSensLoss, logitSqSensLoss, probSensLoss, lossSensLoss = self.computeLosses(
            y, std, regMultiplier, deterministic=False)
            
        # computation of losses for output and cleaning (deterministic)
        classificationLossVal, totalLossVal, l1LossVal, l2LossVal, logitSensLossVal, logitDiffSensLossVal,logitSqSensLossVal, probSensLossVal, lossSensLossVal = self.computeLosses(
            y, std, regMultiplier, deterministic=True)

        params = lasagne.layers.get_all_params(outputLayer,
                                               trainable=True)          
        updates = optimizer(totalLoss, params, learning_rate = learningRate)        
                          
#        testPrediction  = lasagne.layers.get_output(outputLayer,
#                                                    deterministic = True)
#        testAccuracy = T.mean(T.eq(T.argmax(testPrediction, axis=1), y), 
#                          dtype=theano.config.floatX)  
        
        # adaptive noise regularization #
        advNoise = self.hyperp.setdefault('adversarialNoise',0.)
        if advNoise > 0:            
            bXAdap = self.__adversarialNoise(std, advNoise)
            print('Compiling DeepFool graph...')
            deepFoolFunction = theano.function([networkInput], 
                                                bXAdap)
                       
        print('Compiling neural network graph...')       
        trainLoopFunction = theano.function([networkInput, y], 
                                                 totalLoss, 
                                                 updates = updates)  
                                                 
        print('Compiling output graph...')                                                 
        lossOutputs = theano.function([networkInput, y], 
                                       [classificationLossVal,l1LossVal, l2LossVal, 
                                        logitSensLossVal,logitSqSensLossVal, probSensLossVal,
                                        lossSensLossVal])  

        if verbose == 0:
            simpleOutput = theano.function([networkInput, y], 
                                         totalLoss)   
                                
        print('Beginning training loop...')        
        self.classifier.train_history_ = list()        
        for epoch in range(maxEpochs):
            ########################
            #    TRAINING DATA     #
            ########################
            trainBatches = biTrain(data, classes)
            numOfTrainingPoints = len(data)
            out = ""
            for bX,by in trainBatches:
                if advNoise > 0:  # adaptive noise on the network inputs
                    sizeBatch = bX.shape[0]
                    # DeepFool
                    bXFool = np.copy(bX)
                    # create dataset with normal and adversarial examples
                    # Parameter: p
                    p = 1
                    randFool = np.random.binomial(1, p, sizeBatch).astype('bool')
                    bXFool = bXFool[randFool,:,:,:]
                    # call deepFool on subset of samples
                    bXadap = deepFoolFunction(bXFool)
                    # stack noise samples and unchanged samples
                    bXNormal = bX[~randFool,:,:,:]
                    bXCombined = np.vstack((bXadap, bXNormal))
                    # stack labels
                    byFool = by[randFool]
                    byNormal = by[~randFool]
                    byCombined = np.hstack((byFool, byNormal)) 
                    # call training function
                    totalError = trainLoopFunction(bXCombined,byCombined)
                else:
                    totalError = trainLoopFunction(bX,by)
                
            cL=l1L=l2L=loL=loSqL=prL=lsL=0 
            
            if verbose == 1:
                for bX,by in trainBatches:
                    stdout.flush()
                    batchSize = bX.shape[0]
                    cLt,l1Lt,l2Lt,loLt,loSqLt,prLt,lsLt = lossOutputs(bX,by)
                    cL += cLt * batchSize
                    l1L += l1Lt * batchSize
                    l2L += l2Lt * batchSize
                    loL += loLt * batchSize
                    loSqL += loSqLt * batchSize
                    prL += prLt * batchSize
                    lsL += lsLt * batchSize
                
                cL/=numOfTrainingPoints
                l1L/=numOfTrainingPoints
                l2L/=numOfTrainingPoints
                loL/=numOfTrainingPoints
                loSqL/=numOfTrainingPoints
                prL/=numOfTrainingPoints
                lsL/=numOfTrainingPoints
                
                out = "Training set  : cL: %5.5f, l1: %5.5f, l2: %5.5f, logitSens: %5.5f, logitSqSens: %5.5f, probSens: %5.5f, lossSens: %5.5f" % (cL,l1L,l2L,loL,loSqL,prL,lsL)         
                
            if verbose == 0:
                out = "Total loss (last batch): %5.5f" % (totalError)
                stdout.flush()
                
            logDict = dict()
            logDict['cL'] = cL
            logDict['l1'] = l1L
            logDict['l2'] = l2L
            logDict['loL'] = loL
            logDict['loSqL'] = loSqL
            logDict['prL'] = prL
            logDict['lsL'] = lsL
            
            
            print("Epoch ", epoch+1, ":")
            print(out)
            
            if self.hyperp['useSensRegControl']:
                if epoch > .25*maxEpochs and epoch < .75*maxEpochs:
                    t = float(epoch - maxEpochs/4.) / (maxEpochs/2.)
                    regMultiplier.set_value(t)
                if epoch >= .75*maxEpochs:
                    regMultiplier.set_value(1)

            ########################
            #    VALIDATION DATA   #
            ########################
            if validationSetRatio != 0:
                valBatches = biVal(validationData, validationClasses)
                numOfTrainingPoints = len(validationData)
                out = ""
                    
                cL=l1L=l2L=loL=loSqL=prL=lsL=0 
                
                for bX,by in valBatches:
                    stdout.flush()
                    batchSize = bX.shape[0]
                    cLt,l1Lt,l2Lt,loLt,loSqLt,prLt,lsLt = lossOutputs(bX,by)
                    cL += cLt * batchSize
                    l1L += l1Lt * batchSize
                    l2L += l2Lt * batchSize
                    loL += loLt * batchSize
                    loSqL += loSqLt * batchSize
                    prL += prLt * batchSize
                    lsL += lsLt * batchSize
                
                cL/=numOfTrainingPoints
                l1L/=numOfTrainingPoints
                l2L/=numOfTrainingPoints
                loL/=numOfTrainingPoints
                loSqL/=numOfTrainingPoints
                prL/=numOfTrainingPoints
                lsL/=numOfTrainingPoints

                logDict['val_cL'] = cL
                logDict['val_l1'] = l1L
                logDict['val_l2'] = l2L
                logDict['val_loL'] = loL
                logDict['val_loSqL'] = loSqL
                logDict['val_prL'] = prL
                logDict['val_lsL'] = lsL
                
                out = "Validation set: cL: %5.5f, l1: %5.5f, l2: %5.5f, logitSens: %5.5f, logitSqSens: %5.5f, probSens: %5.5f, lossSens: %5.5f" % (cL,l1L,l2L,loL,loSqL,prL,lsL)         
                print(out)
                
                
            self.classifier.train_history_.append(logDict)
            
            # Check if the current training epoch is the best model
            #getBestValidNet(self.classifier)
            
            
        return NNModel_theano(self.classifier)
               
    def __adversarialNoise(self,  std, advNoise):
        """
        TODO: currently only works for 2 classes !!!
        """
        layers = self.classifier.layers_    
        outputLayer = layers[-2]
        inputLayer = layers[0]
        networkInput = inputLayer.input_var 
        networkOutput = get_output(outputLayer, deterministic=True)
        # difference of logits (only works for 2-classes !!!!)
        fX = networkOutput[:,0] -networkOutput[:,1]    
        g_logit = T.grad(T.sum(fX), networkInput)
        # scale gradient by std
        if std is not None:
            g_logit = std * g_logit
        # l2-norm squared    
        grad_l2 = T.sum(T.sqr(g_logit), axis = (1,2,3)) + 1e-12
        # scale again by gradient due to divison by squred l2-norm of gradient
        if std is not None:
            sqrG_logit = std * g_logit  
        r = - (sqrG_logit.dimshuffle(1,2,3,0) * (fX / grad_l2)).dimshuffle(3,0,1,2)
        # random scaling of noise
        randTensor = self._srng.uniform(size=(networkInput.shape[0],), low=-advNoise/2. ,high=advNoise)       
        addedTerm = (r.dimshuffle(2,3,1,0) * randTensor).dimshuffle(3,2,0,1)
        # projection onto non-negative values
        bXFool = T.maximum(networkInput + addedTerm, 0)
        return bXFool                       
                                     
        
    def computeLosses(self, y, std, regMultiplier, deterministic):

        logitSens = self.hyperp.setdefault('logitSens',0.)
        logitDiffSens = self.hyperp.setdefault('logitDiffSens',0.)
        logitSqSens = self.hyperp.setdefault('logitSqSens',0.)
        probSens = self.hyperp.setdefault('probSens',0.)
        lossSens = self.hyperp.setdefault('lossSens',0.)
        l1 = self.hyperp.setdefault('l1',0.)
        l2 = self.hyperp.setdefault('l2',0.)
         
        layers = self.classifier.layers_
        
        lossFunction = lasagne.objectives.categorical_crossentropy
        
        aggregate = T.mean # otherwise lasagne.objectives.aggregate
        
        outputLayer = layers[-1]
        logitLayer = layers[-2]
        inputLayer = layers[0]
        networkInput = inputLayer.input_var
        networkOutput = get_output(outputLayer, deterministic=deterministic)
        logitOutput = get_output(logitLayer, deterministic=deterministic)

        ######################################################################
        # Very weird thing:
        # lossSensitivity gradients can only be computed if the one-hot encoded
        # version of the loss function is used. BUT that version lacks a
        # stability optimization in Theano that leads to NaNs during training.
        # This is why both versions need to be employed here. 
        
        L = lossFunction(networkOutput,y)
    
        y_oneHot = lasagne.utils.one_hot(y,outputLayer.output_shape[1])
        L_oneHot = lossFunction(networkOutput, y_oneHot)

        #######################################################################
        
        classificationLoss = aggregate(L)        
        

        l1Loss = regularization.regularize_layer_params(
            layers.values(), regularization.l1)
         

        l2Loss = regularization.regularize_layer_params(
            layers.values(), regularization.l2)
    
        # logit sensitivity      
        logit = T.sum(logitOutput * y_oneHot, axis=1)
        G_logit = T.grad(T.sum(logit), networkInput)

        if std is not None:
            G_logit = std * G_logit
        
        # Sparse logit saliency regularization
        absG_logit = T.abs_(G_logit)
        sumAbsG_logit = T.sum(absG_logit,axis = (1,2,3))
        logitSensLoss = aggregate(sumAbsG_logit) 

        # Squared logit saliency regularization
        sqG_logit = G_logit**2
        sumSqG_logit = T.sum(sqG_logit,axis = (1,2,3))
        logitSqSensLoss = aggregate(sumSqG_logit) 
    
        # probability sensitivity  
        prob = T.sum(networkOutput * y_oneHot, axis=1)
        G_prob = T.grad(T.sum(prob), networkInput)

        if std is not None:
            G_prob = std * G_prob
        
        # Sparse probability saliency regularization
        absG_prob = T.abs_(G_prob)
        sumAbsG_prob = T.sum(absG_prob,axis = (1,2,3))
        probSensLoss = aggregate(sumAbsG_prob) 
        
        # Loss sensitivity
        G_loss = theano.grad(T.sum(L_oneHot),networkInput)
        if std is not None:
            G_loss = std * G_loss
        absG_loss = T.abs_(G_loss)
        sumAbsG_loss = T.sum(absG_loss,axis = (1,2,3))
        lossSensLoss = aggregate(sumAbsG_loss)


        ####### !!!!!!!!!!!!!!!!!!! EXPERIMENTAL !!!!!!!!!!!!!!!!!! ##########
        #### !!!! only makes sense for 2-class problems in this case !!!! ####

        # Clumsy way to regularize logit differences
        # It works by replacing the matrix of one-hot encoded labels by one
        # whose first column is ones and the rest is minus ones. After summing 
        # over each row, we are left with the difference of the logit of the
        # first class and the (sum of the) other class(es).

        plusMinusOneMatrix = 2 * lasagne.utils.one_hot(1,outputLayer.output_shape[1]) - T.ones_like(y_oneHot)
        logitDiff = T.sum(logitOutput * plusMinusOneMatrix, axis=1)
        G_logitDiff = T.grad(T.sum(logitDiff), networkInput)

        if std is not None:
            G_logitDiff = std * G_logitDiff

        absG_logitDiff = T.abs_(G_logitDiff)
        sumAbsG_logitDiff = T.sum(absG_logitDiff,axis = (1,2,3))
        logitDiffSensLoss = aggregate(sumAbsG_logitDiff) 

        
        # Sum up
        totalLoss = classificationLoss
        if l1 : totalLoss += regMultiplier * l1 * l1Loss
        if l2 : totalLoss += regMultiplier * l2 * l2Loss
        if logitSens : totalLoss += regMultiplier * logitSens * logitSensLoss
        if logitDiffSens : totalLoss += regMultiplier * logitDiffSens * logitDiffSensLoss  
        if logitSqSens : totalLoss += regMultiplier * logitSqSens * logitSqSensLoss
        if probSens : totalLoss += regMultiplier * probSens * probSensLoss
        if lossSens : totalLoss += regMultiplier * lossSens * lossSensLoss     
        
        return classificationLoss, totalLoss, l1Loss, l2Loss, logitSensLoss, logitDiffSensLoss, logitSqSensLoss, probSensLoss, lossSensLoss

# Normally, nolearn only saves the total loss (avg NLL + l1/l2-Regularization).
# This method splits this loss into its components.
def getIndividualLosses(net,history):
    reg = lasagne.regularization.regularize_layer_params
    layers = net.get_all_layers()
    if net.objective_l1:
        net.train_history_[-1]['l1Loss'] = reg(layers,
            lasagne.regularization.l1).eval() * net.objective_l1
    else:
        net.train_history_[-1]['l1Loss'] = 0

    if net.objective_l2:
        net.train_history_[-1]['l2Loss'] = reg(layers,
            lasagne.regularization.l2).eval() * net.objective_l2
    else:
        net.train_history_[-1]['l2Loss'] = 0

    if net.objective_tv:
        net.train_history_[-1]['tvLoss'] = 0
    else:
        net.train_history_[-1]['tvLoss'] = 0

# This method adds the system time of the end of each epoch to the history log.
def addEndTimeToHistory(net, history):
    net.train_history_[-1]['endTime'] = time.time()

# Add to on_epoch_finished in order to let the learning rate decay with each
# iteration.
def decay(net, history):
    decayRate = .01
    net.update_learning_rate = net.update_learning_rate * (1-decayRate)


# Prints the error in each epoch.
def printError(net,history):
    it = len(history)
    tL = history[-1]['train_loss']
    vL = history[-1]['valid_loss']
    vA = history[-1]['valid_accuracy']

    l1 = net.train_history_[-1]['l1Loss']
    l2 = net.train_history_[-1]['l2Loss']
    tv = net.train_history_[-1]['tvLoss']

    # The easiest way to obtain the average categorical crossentropy error per
    # epoch is to simply substract the l1 and l2 penalty terms from the total
    # loss.
    # nolearn calculates the training loss as the running average over all
    # batches per epoch. This leads to the validation loss often appearing to
    # be lower than the training loss.
    true_tL = tL - l1 - l2 - tv
    true_vL = vL - l1 - l2 - tv

    net.train_history_[-1]['trainNLL'] = true_tL
    net.train_history_[-1]['validNLL'] = true_vL

    out = "Epoch %4d , NLL(train): %5.5f" % (it,true_tL)

    if true_vL and not np.isnan(true_vL):
        out = out + ", NLL(valid): %5.5f" % (true_vL)
        out = out + ", Accuracy(valid): %5.5f" % (vA)
    if l1 and not np.isnan(l1):
        out = out + ", l1: %5.5f" % (l1)
    if l2 and not np.isnan(l2):
        out = out + ", l2: %5.5f" % (l2)
    print('\r' + out)

    # Without this, the output isn't show immediately
    stdout.flush()

# Print the adaptive noise in each epoch
def printAdaptiveNoise(net, history):
    try:
        noiseLevel = net.get_all_layers()[1].stochasticNoiselevel.eval()
        net.train_history_[-1]['adaptiveNoiseLevel'] = noiseLevel
        adaptiveNoiseNorm = np.linalg.norm(noiseLevel)
        net.train_history_[-1]['noise level norm'] = adaptiveNoiseNorm
        print('\r' + 'Noise level norm ' + str(adaptiveNoiseNorm))
    except:
        pass

# Save the best net (lowest NLL on validation set) during training
def saveBestValidNet(net, history):
    if net.train_split.eval_size > 0:
        # get NLL loss and penalty term losses
        vL = history[-1]['validNLL']
        l1 = net.train_history_[-1]['l1Loss']
        l2 = net.train_history_[-1]['l2Loss']
        tv = net.train_history_[-1]['tvLoss']
        # compute NLL loss
        true_vL = vL - l1 - l2 - tv
        try:
            # current iteration is best model
            if net.train_history_[0]['bestValidationNLL'] > true_vL:
                net.train_history_[0]['bestModel'] = net.\
                    get_all_params_values()
                net.train_history_[0]['bestValidationNLL'] = true_vL
                print('Best model saved!')
        except:  # bestValidationNLL variable not instantiated
            net.train_history_[0]['bestValidationNLL'] = np.inf
            
            
# Get the best parameter configuration of the current training process
def getBestValidNet(net):
    # get NLL loss and penalty term losses
    print(net.train_history_.keys())
    cL    = net.train_history_[-1]['validNLL']
#    l1L   = net.train_history_[-1]['val_l1']
#    l2L   = net.train_history_[-1]['val_l2']
#    loL   = net.train_history_[-1]['val_loL']
#    loSqL = net.train_history_[-1]['val_loSqL']
#    prL   = net.train_history_[-1]['val_prL']
#    lsL   = net.train_history_[-1]['val_lsL']  
    # compute NLL loss  
    # TODO: Check if all losses needed
    #NLL_val = cL - l1L - l2L - loL - loSqL - prL - lsL
    NLL_val = cL
    
    try:
        # current iteration is best model
        if net.train_history_[0]['bestValidationNLL'] > NLL_val:
            net.bestModel = net.get_all_params_values()
            net.train_history_[0]['bestValidationNLL'] = NLL_val
            cur_epoch = len(net.train_history_)
            print('Found new best model! Epoch: ' + str(cur_epoch))
    except:  # bestValidationNLL variable not instantiated
        net.train_history_[0]['bestValidationNLL'] = np.inf
              

class ModifiedNeuralNet(NeuralNet):

    def __init__(
        self,
        layers,
        update=nesterov_momentum,
        loss=None,  # BBB
        objective=objective,
        objective_loss_function=None,
        batch_iterator_train=BatchIterator(batch_size=128),
        batch_iterator_test=BatchIterator(batch_size=128),
        regression=False,
        max_epochs=100,
        train_split=TrainSplit(eval_size=0.2),
        custom_scores=None,
        scores_train=None,
        scores_valid=None,
        X_tensor_type=None,
        y_tensor_type=None,
        use_label_encoder=False,
        on_batch_finished=None,
        on_epoch_finished=None,
        on_training_started=None,
        on_training_finished=None,
        more_params=None,
        check_input=True,
        verbose=0,
        **kwargs
        ):

        NeuralNet.__init__(self,
            layers,
            update,
            loss,
            objective,
            objective_loss_function,
            batch_iterator_train,
            batch_iterator_test,
            regression,
            max_epochs,
            train_split,
            custom_scores,
            scores_train,
            scores_valid,
            X_tensor_type,
            y_tensor_type,
            use_label_encoder,
            on_batch_finished,
            on_epoch_finished,
            on_training_started,
            on_training_finished,
            more_params,
            check_input,
            verbose,
            **kwargs)


    def _create_iter_funcs(self, layers, objective, update, output_type):
        y_batch = output_type('y_batch')

        objective_kw = self._get_params_for('objective')

        loss_train = objective(
            layers, target=y_batch, **objective_kw)
        loss_eval = objective(
            layers, target=y_batch, deterministic=True, **objective_kw)

        output_layer = self._output_layers
        predict_proba = get_output(output_layer, None, deterministic=True)
        if not self.regression:
            predict = predict_proba[0].argmax(axis=1)
            accuracy = T.mean(T.eq(predict, y_batch))
        else:
            accuracy = loss_eval

        scores_train = [
            s[1](predict_proba, y_batch) for s in self.scores_train]
        scores_valid = [
            s[1](predict_proba, y_batch) for s in self.scores_valid]

        all_params = self.get_all_params(trainable=True)
        grads = T.grad(loss_train, all_params)
        for idx, param in enumerate(all_params):
            grad_scale = getattr(param.tag, 'grad_scale', 1)
            if grad_scale != 1:
                grads[idx] *= grad_scale
        update_params = self._get_params_for('update')
        updates = update(grads, all_params, **update_params)

        input_layers = [layer for layer in layers.values()
                        if isinstance(layer, InputLayer)]

        X_inputs = [theano.In(input_layer.input_var, name=input_layer.name)
                    for input_layer in input_layers]
        inputs = X_inputs + [theano.In(y_batch, name="y")]

        train_iter = theano.function(
            inputs=inputs,
            outputs=[loss_train] + scores_train,
            updates=updates,
            allow_input_downcast=True,
            )
        eval_iter = theano.function(
            inputs=inputs,
            outputs=[loss_eval, accuracy] + scores_valid,
            allow_input_downcast=True,
            )
        predict_iter = theano.function(
            inputs=X_inputs,
            outputs=predict_proba,
            allow_input_downcast=True,
            )

        return train_iter, eval_iter, predict_iter


# This layer realizes a pre-defined filtering with several filter:
# - gaussian smoothing with different variances
# - constant smoothing (5x1 and 3x1)
# - identity filtering
# (Note: constant filter 2 and identity only work for size 5x1)
class SmoothingLayer(lasagne.layers.Layer):
    """
    Convolutional layer with pre-defined smoothing filters as well as an
    identity mapping.
    :param incoming: incoming layer
    :param shape: TODO
    :param trainingSize: TODO
    """
    def __init__(self, incoming, border_mode='valid', **kwargs):
        self.border_mode = border_mode
        self.num_filter = 6
        super(SmoothingLayer, self).__init__(incoming, **kwargs)
        self.ds = (5, 1)
        # container for all kernels (6 different kinds here)
        cker = np.asanyarray(np.zeros((self.num_filter, 1, self.ds[0], self.ds[1])),
                             dtype='float32')
        # gaussian smoothing (small variance 0.5)
        from scipy.signal import get_window
        gauss1 = get_window(('gaussian', 0.5), self.ds[0])
        cker[0, 0] = np.reshape(np.asanyarray(gauss1 / np.sum(gauss1),
                                dtype='float32'), (self.ds[0],1))
        # gaussian smoothing (medium variance 1)
        gauss2 = get_window(('gaussian', 1), self.ds[0])
        cker[1, 0] = np.reshape(np.asanyarray(gauss2 / np.sum(gauss2),
                                 dtype='float32'), (self.ds[0],1))
        # gaussian smoothing (large variance 2)
        gauss3 = get_window(('gaussian', 2), self.ds[0])
        cker[2, 0] = np.reshape(np.asanyarray(gauss3 / np.sum(gauss3),
                                dtype='float32'), (self.ds[0],1))
        # constant smoothing (size as pre-defined)
        constant1 = np.ones(self.ds)
        cker[3, 0] = np.asanyarray(constant1 / self.ds[0], dtype='float32')
        # constant smoothing smaller
        constant2 = np.reshape(np.array([0.,  1.,  1.,  1.,  0.]),
                              (self.ds[0],1))
        cker[4, 0] = np.asanyarray(constant2 / 3., dtype='float32')
        # identity filter
        identity = np.reshape(np.array([0.,  0.,  1.,  0.,  0.]),
                              (self.ds[0],1))
        cker[5, 0] = np.asanyarray(identity, dtype='float32')

        self.sym_k = theano.shared(cker, 'kernels')

    def get_output_shape_for(self, input_shape):
        output_shape = list(input_shape)
        output_shape[1] = self.num_filter
        output_shape[2] = input_shape[2] - self.ds[0] + 1
        output_shape[3] = input_shape[3] - self.ds[1] + 1
        return tuple(output_shape)


    def get_output_for(self, input, **kwargs):
        return dnn_conv(input, self.sym_k, border_mode=self.border_mode)
        
# This layer realizes a pre-defined filtering with one filter:
# - constant smoothing (3x1)
class SmoothingLayerOnes(lasagne.layers.Layer):
    """
    Convolutional layer with pre-defined smoothing filters as well as an
    identity mapping.
    :param incoming: incoming layer
    :param shape: TODO
    :param trainingSize: TODO
    """
    def __init__(self, incoming, border_mode='valid', **kwargs):
        self.border_mode = border_mode
        self.num_filter = 1
        super(SmoothingLayerOnes, self).__init__(incoming, **kwargs)
        self.ds = (3, 1)
        # container for all kernels (6 different kinds here)
        cker = np.asanyarray(np.zeros((self.num_filter, 1, self.ds[0], self.ds[1])),
                             dtype='float32')
        # constant smoothing (size as pre-defined)
        constant1 = np.ones(self.ds)
        cker[0, 0] = np.asanyarray(constant1 / self.ds[0], dtype='float32')


        self.sym_k = theano.shared(cker, 'kernels')

    def get_output_shape_for(self, input_shape):
        output_shape = list(input_shape)
        output_shape[1] = self.num_filter
        output_shape[2] = input_shape[2] - self.ds[0] + 1
        output_shape[3] = input_shape[3] - self.ds[1] + 1
        return tuple(output_shape)


    def get_output_for(self, input, **kwargs):
        return dnn_conv(input, self.sym_k, border_mode=self.border_mode)

# Custom batch iterator for data augmentation
class DataAugmentationBatchIterator(BatchIterator):
    """
    Custom batch iterator for data augmentation.

    :param batch_size: batch size
    :param shuffle: whether to shuffle the data before every epoch
    :param disturbLabelRate: 0..1, rate of randomly assigned labels
    :param sdWidth: width of the structured dropout patches
    :param sdNumber: number of structured dropout patches
    :param seed: random seed
    """
    def __init__(self, batch_size, shuffle=False, disturbLabelRate=0,
                 sdWidth=30,sdNumber=0,seed=42):
        self.disturbLabelRate = disturbLabelRate
        self.sdWidth = sdWidth
        self.sdNumber = sdNumber
        BatchIterator.__init__(self, batch_size,shuffle, seed)

    def transform(self, X, y):
        yb = y.copy()
        Xb = X.copy()

        # Implements DisturbLabel.
        # This randomly assigns wrong labels to points in the data set
        # as a regularization method.
        # arXiv:1605.00055
        if y is not None:
            maxLabel = max(y)
            m = np.random.randint(y.shape[0],size=np.floor(y.shape[0]*\
            self.disturbLabelRate).astype(int))
            m = np.unique(m)
            yb[m] = np.random.randint(maxLabel+1,size=m.shape[0])

        # Implements 'structured dropout'.
        # This sets sdNumber random portions of width sdWith to 0.
        # This is a form of data augmentation.
        for k in range(self.sdNumber):
            if self.sdWidth is not 0:
                for i in range(Xb.shape[0]):
                    startPoint = np.random.randint(0,Xb.shape[2]-self.sdWidth)
                    Xb[i,0,startPoint:startPoint+self.sdWidth,0] = 0
        return Xb, yb


# This method applies a random elastic transformation to the data.
from scipy.ndimage.interpolation import map_coordinates
from scipy.ndimage.filters import gaussian_filter1d
def elasticTransform(spectrum, alpha, sigma, random_state=None):
    """
    This method applies a random elastic transformation to the spectrum.
    This is achieved by creating a random vector field of shifts, which is
    smoothed ahead of the transformation.

    :param spectrum: input spectrum
    :param alpha: strength of the smoothing operation
    :param sigma: variance of the randomized initial vector field
    """
    if random_state is None:
        random_state = np.random.RandomState(None)

    shape = (spectrum.shape[0],spectrum.shape[2])

    dy = gaussian_filter1d((random_state.rand(*shape) * 2 - 1),
                           sigma, mode="constant", cval=0) * alpha
    x, y = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), indexing='ij')

    indices = np.reshape(x, (-1, 1)), np.reshape(y+dy, (-1, 1))
    coordinates = map_coordinates(spectrum.reshape(shape), indices,
                                  order=1).reshape((shape[0],1,shape[1],1))
    return coordinates


class ElasticBatchIterator(BatchIterator):
    """
    Batch iterator which implements the elastic transformation.

    :param batch_size: batch size
    :param shuffle: whether to shuffle the data before every epoch
    """
    def __init__(self, batch_size, shuffle=False, seed=42):
        BatchIterator.__init__(self, batch_size,shuffle, seed)

    def transform(self, X, y):
        Xb = elasticTransform(X,.9,.5)
        return Xb, y
