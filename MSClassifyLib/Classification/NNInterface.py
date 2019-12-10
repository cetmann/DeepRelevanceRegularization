#!/usr/bin/env python
"""
This file implements the interface between the MSClassifyLib (MATLAB) and
the Python neural network libraries (Theano, tensorflow.keras, Lasagne, Nolearn etc.).
Furthermore, a log-file of the training procedure is created.
"""

import numpy as np

# The following libraries and methods are used for the file transfer between
# MATLAB and Python.
import sys
import os
#import pickle, cPickle
import h5py
try:
    import cPickle
except ImportError:
    import pickle as cPickle

import argparse


if __name__ == '__main__':

    # The standard recursion limit may be too low, so set this higher.
    sys.setrecursionlimit(50000)

    # Avoid bytecode clutter
    sys.dont_write_bytecode = True

    # Default values
    createClassDictionary = False
    createNewNetwork = False
    detailedPrediction = False

    # Parse arguments from the MATLAB classifier
    parser = argparse.ArgumentParser(
        description=('Interface between MATLAB and Python for neural network '
                     'classificiation of MALDI mass spectra.'))
    parser.add_argument('-neuralNetFile', nargs=1,
                        default='neuralNetFile.nnet',
                        help=('Name of the file the serialized neural '
                              'network is later stored in.'))
    parser.add_argument('-dataFile', nargs=1,
                        default='dataFile.mat',
                        help=('Name of the temporary file the data is '
                              'stored in.'))
    parser.add_argument('-resultFile', nargs=1,
                        default='resultFile.mat',
                        help=('Name of the temporary file the results '
                              'are stored in.'))
    parser.add_argument('-backend', nargs=1,
                        default='theano',
                        help=('Name of the backend. Default: "theano"'))
    parser.add_argument('-neuralNetArchitecture', nargs=1,
                        default='auto',
                        help='Which architecture to use. Default: "auto"')
    parser.add_argument('-mode', nargs=1, default='train',
                        help='"train", "predict", "sensitivity", "reconstruction" or "reconstructionBasis"')
                        
    ### Begin: reconstruction arguments ###
    parser.add_argument('-weighting', nargs=1, default=['1'],
                        help='Decides whether weighting by the classifier is performed')
    parser.add_argument('-reconLayer', nargs=1, default=['1'],
                        help='The layer from which reconstruction is performed')
    parser.add_argument('-reconStepSize', nargs=1, default=['0.1'],
                        help='The step size of gradient descent algorithm.') 
    parser.add_argument('-reconNumIterations', nargs=1, default=['10'],
                        help='The number of iterations of gradient descent algorithm.')
    parser.add_argument('-predCoeff', nargs=1, default=['0.1'],
                        help='The coefficient for NLL-penalty.')                    
    parser.add_argument('-reconRegCoeffL1', nargs=1, default=['0.1'],
                        help='The coefficient for L1-penalty.')
    parser.add_argument('-reconRegCoeffL2', nargs=1, default=['0.1'],
                        help='The coefficient for L2-penalty.') 
    ### End: reconstruction arguments ###

    parser.add_argument('-weightClasses', action='store_true',
                        help=('Whether to weight the classes based on prevalence'))
    parser.add_argument('-sensType', nargs=1, default='binary',
                        help='"NLL", "perClass" or "binary"')
    parser.add_argument('-overwrite', action='store_true',
                        help=('Whether to overwrite a network trained under '
                        'the same neuralNetFile file name. If not, the '
                        'network is continued to be trained.'))
    parser.add_argument('-useValidationSet', action='store_true',
                        help=('Whether to use 10 per cent of the training '
                        'data for validation.'))
    parser.add_argument('-useSensRegControl', action='store_true',
                        help=('Whether to turn on the saliency regularization '
                        'slowly after 25 per cent of the number of epochs.'))
    parser.add_argument('-saveBestValidModel', action='store_true',
                        help=('Whether to save and restore best net on '
                        'validation set.'))
    parser.add_argument('-learningRate', nargs=1, default=['0.001'],
                        help='The learning rate during training.')
    parser.add_argument('-lrSchedule', nargs=1, default=['0'],
                        help='The period when the learning rate is divided by 10.')
    parser.add_argument('-validationSetRatio', nargs=1, default=['0.0'],
                        help=('The relative amount of training samples which '
                        'are used for validation.'))
    parser.add_argument('-disturbLabelRate', nargs=1, default=['0.0'],
                        help=('The relative amount of labels that is changed '
                        'randomly during training.'))
    parser.add_argument('-adversarialNoise', nargs=1, default=['0.0'],
                        help=('The scaling factor of the adversarial noise.')) 
    parser.add_argument('-l1', nargs=1, default=['0.0'],
                        help=('The value of the l1-regularization parameter.'))
    parser.add_argument('-l2', nargs=1, default=['0.0'],
                        help='The value of the l2-regularization parameter.')
    parser.add_argument('-logitSens', nargs=1, default=['0.0'],
                        help=('The value of the logit sensitivity '
                        'regularization parameter.'))
    parser.add_argument('-logitDiffSens', nargs=1, default=['0.0'],
                        help=('The value of the logit difference sensitivity '
                        'regularization parameter.'))
    parser.add_argument('-logitSqSens', nargs=1, default=['0.0'],
                        help=('The value of the squared logit sensitivity '
                        'regularization parameter.'))
    parser.add_argument('-probSens', nargs=1, default=['0.0'],
                        help=('The value of the probability sensitivity '
                        'regularization parameter.'))
    parser.add_argument('-lossSens', nargs=1, default=['0.0'],
                        help=('The value of the loss sensitivity '
                        'regularization parameter.'))
    parser.add_argument('-trainingDataStd', nargs=1, default=['data_std'],
                        help=('Standardization factor for some regularization '
                        'methods.'))
    parser.add_argument('-lossFunc', nargs=1, default=['categorical_crossentropy'],
                        help=('The main loss for the model.'))
    parser.add_argument('-epochs', nargs=1, default=['150'],
                        help='The number of epochs.')
    parser.add_argument('-batchSize', nargs=1, default=['128'],
                        help='The size of each mini-batch.')
    parser.add_argument('-sdWidth', nargs=1, default=['0'],
                        help=('The number of connected m/z-bins which are set '
                        'to 0 with structured dropout.'))
    parser.add_argument('-sdNumber', nargs=1, default=['0'],
                        help='This parameter determines how many times per '
                        'spectrum structured dropout is used.')
    parser.add_argument('-verbose', nargs=1, default=['1.0'],
                        help='This parameter determines the amount of '
                        'information displayed.')

    args = parser.parse_args()
    neuralNetFile = args.neuralNetFile[0]
    dataFile = args.dataFile[0]
    resultFile = args.resultFile[0]
    neuralNetArchitecture = args.neuralNetArchitecture[0]
    backend = args.backend[0]
    mode = args.mode[0]
    typeSens = args.sensType[0]
    overwrite = args.overwrite

    if args.useValidationSet:
        validationSetRatio = np.float(args.validationSetRatio[0])
    else:
        validationSetRatio = .0

    # Set hyperparameters from the input
    hyperparameters = dict()
    hyperparameters['verbose'] = np.float(args.verbose[0])
    hyperparameters['lrSchedule'] = np.float(args.lrSchedule[0])
    hyperparameters['learningRate'] = np.float(args.learningRate[0])
    hyperparameters['disturbLabelRate'] = np.float(args.disturbLabelRate[0])
    hyperparameters['adversarialNoise'] = np.float(args.adversarialNoise[0])
    hyperparameters['l1'] = np.float(args.l1[0])
    hyperparameters['l2'] = np.float(args.l2[0])
    hyperparameters['logitSens'] = np.float(args.logitSens[0])
    hyperparameters['logitDiffSens'] = np.float(args.logitDiffSens[0])
    hyperparameters['logitSqSens'] = np.float(args.logitSqSens[0])
    hyperparameters['probSens'] = np.float(args.probSens[0])
    hyperparameters['lossSens'] = np.float(args.lossSens[0])
    hyperparameters['batchSize'] = np.int(np.float(args.batchSize[0]))
    hyperparameters['epochs'] = np.int(np.float(args.epochs[0]))
    hyperparameters['sdWidth'] = np.int(np.float(args.sdWidth[0]))
    hyperparameters['sdNumber'] = np.int(np.float(args.sdNumber[0]))
    hyperparameters['validationSetRatio'] = validationSetRatio
    hyperparameters['useSensRegControl'] = args.useSensRegControl
    hyperparameters['trainingDataStd'] = args.trainingDataStd[0]
    hyperparameters['lossFunc'] = args.lossFunc[0]


    # set parameters for reconstruction
    paramRecon = dict()
    paramRecon['weighting'] = bool(np.int(args.weighting[0]))
    paramRecon['reconLayer'] = np.int(np.float(args.reconLayer[0]))
    paramRecon['reconStepSize'] = np.float(args.reconStepSize[0])
    paramRecon['reconNumIterations'] = np.int(np.float(args.reconNumIterations[0]))
    paramRecon['predCoeff'] = np.float(args.predCoeff[0])
    paramRecon['reconRegCoeffL1'] = np.float(args.reconRegCoeffL1[0])
    paramRecon['reconRegCoeffL2'] = np.float(args.reconRegCoeffL2[0])


    # Make backend specific imports 
    if backend == 'keras':
        import tensorflow as tf
        from tensorflow.keras import backend as K
        K.set_image_data_format('channels_first')
        from tensorflow.keras.models import load_model


    # Load data from the dataFile
    if 'dataFile' in globals():
        dataSet = h5py.File(dataFile)
    else:
        raise ValueError('No dataFile was provided.')

    # Load the class labels and convert them to int32 (which is required for
    # GPU parallelization with CUDA).
    if mode == 'train' or mode == 'reconstruction' or typeSens == 'NLL':
        classes = dataSet['classes'].value.astype(np.int32).reshape(-1)
    else:
        classes = None
    # load basis for reconstruction option
    basisRecon = False
    if mode == 'reconstructionBasis': # EDIT THIS OUT TODO       
        basisRecon = dataSet['basisRecon'].value.astype(np.float32)

    # The standard implementation of the rectified linear unit
    # (theano.tensor.nnet.rectify) may lead to NaN gradients. This is because
    # with cuDNN, the computational graph for the subderivative relu'(x)=sgn(x)
    # is transformed into x/abs(x), which leads to NaN for x=0. The following
    # implementation circumvents this. It needs to be loaded before unpacking
    # the data contained in the neuralNetFile, otherwise cPickle throws 
    # an error.
    if backend == 'theano':
        from theano import tensor
        def relu(x):
            return tensor.maximum(x, 0.)



    ######################## Data loading/processing #########################

    def loading_Model(backend, neuralNetFileContents, neuralNetFile):
        if backend == 'theano':
            neuralNetworkModel = neuralNetFileContents['neuralNetworkModel']
        
        elif backend == 'keras':
            from NNHelpers_keras import get_custom_layer_dict
            from NNClassifier_keras import MaldiRegularizer, get_main_loss

            custom_elements = get_custom_layer_dict()
            custom_regularizer = {'MaldiRegularizer' : MaldiRegularizer}

            custom_loss = {'calc_loss': get_main_loss(hyperparameters.get('lossFunc'))}

            custom_metrics = {str(custom_loss['calc_loss']):custom_loss['calc_loss'],
                                'logitSens':custom_loss['calc_loss'],
                                'logitSqSens':custom_loss['calc_loss'],
                                'probSens':custom_loss['calc_loss'],
                                'lossSens':custom_loss['calc_loss'],
                                'l1Loss':custom_loss['calc_loss'],
                                'l2Loss':custom_loss['calc_loss']}
            
            neuralNetworkModel = load_model(os.path.basename(neuralNetFile) + '_model.h5',
                custom_objects={**custom_elements, **custom_regularizer, **custom_loss, **custom_metrics})

        return neuralNetworkModel


    if mode == 'train':
        neuralNetworkModel = None
        if os.path.exists(neuralNetFile):
            if overwrite:
                os.remove(neuralNetFile)
                os.mknod(neuralNetFile)
                createClassDictionary = True
                createNewNetwork = True
            else:
                f = open(neuralNetFile, 'rb')
                neuralNetFileContents = cPickle.load(f)
                f.close()
                newClasses = neuralNetFileContents['newClasses']
                oldClasses = neuralNetFileContents['oldClasses']
                numUniqueClasses = \
                    neuralNetFileContents['numUniqueClasses']
                neuralNetworkModel = loading_Model(backend, 
                    neuralNetFileContents, neuralNetFile)
                
                createClassDictionary = False
                createNewNetwork = False

        else:
            os.mknod(neuralNetFile)
            createClassDictionary = True
            createNewNetwork = True

    else:
        if not os.path.exists(neuralNetFile):
            raise ValueError('neuralNetFile does not exist.')

        f = open(neuralNetFile, 'rb')
        neuralNetFileContents = cPickle.load(f)
        f.close()
        newClasses = neuralNetFileContents['newClasses']
        oldClasses = neuralNetFileContents['oldClasses']
        numUniqueClasses = neuralNetFileContents['numUniqueClasses']
        neuralNetworkModel = loading_Model(backend, 
                    neuralNetFileContents, neuralNetFile)
        createClassDictionary = False
        createNewNetwork = False

    # Classes are required to be in 0..N-1 classes. The next lines create 
    # dictionary mappings that can convert the old class names to new ones and 
    # vice versa.
    if createClassDictionary:
        uniqueClasses = np.unique(classes)
        numUniqueClasses = len(uniqueClasses)
        newClasses = dict()
        oldClasses = dict()
        for i in range(len(uniqueClasses)):
            newClasses.update({uniqueClasses[i] : i})
            oldClasses.update({i : uniqueClasses[i]})
        classes = np.array([newClasses[classItem] for classItem in classes],
                           dtype=np.int32)


    # Load the data and convert them to float32 (which is required for
    # GPU parallelization with CUDA).
    data = dataSet['data'].value.astype(np.float32)
    
    # Transposing is necessary when using HDF5 between MATLAB and Python.
    data = data.T
    
    # Reshape the data to the following format, which is standard for 
    # convolutional neural networks for image processing:
    # (number of Samples, number of feature maps, c1, c2)
    data = data.reshape((data.shape[0], 1, -1, 1))

    dataShape = data.shape
    dataSet.close()


    # Set the standardization factor
    # It is used in some regularization methods
    if mode == 'train':   
        if backend == 'theano' and not  hyperparameters['trainingDataStd'] == 'data_std':
            print('Changed your regularization standardization factor to: ' + 'data_std'
                + ' since you are using the theano backend.')
            hyperparameters['trainingDataStd'] = 'data_std'

        # Use the standard deviation of the data as a factor.
        # This is the only supported method for the theano backend!
        if hyperparameters['trainingDataStd'] == 'data_std':
            trainingDataStd = data.std(0)
            hyperparameters['trainingDataStd'] = trainingDataStd

        trainingDataStd = hyperparameters['trainingDataStd']
    
    elif mode == 'predict':
        # Load the standardization factor/ methods from the save file
        trainingDataStd = neuralNetFileContents['trainingDataStd']

    elif mode == 'sensitivity':
        trainingDataStd = np.ones_like(data.std(0))

    # TODO Structure this in train, predict, etc.
    # Start the work with the model. This is the first time we have to
    # use specific code regarding the chosen backend
    if backend == 'theano':
        import NNWorker_theano
        #raise DeprecationWarning('Most of the packages for the theano backend \
        #are no longer maintained! Consider switching to the keras backend in the future.')
        NNWorker_theano.startWork(createNewNetwork, neuralNetArchitecture, dataShape, numUniqueClasses,
            hyperparameters, mode, data, classes, args, newClasses, oldClasses, relu,
            trainingDataStd, neuralNetFile, resultFile, typeSens, paramRecon, basisRecon, 
            neuralNetworkModel)
        
    elif backend == 'keras':
        import NNWorker_keras
        NNWorker_keras.startWork(createNewNetwork, neuralNetArchitecture, dataShape, numUniqueClasses,
            hyperparameters, mode, data, classes, args, newClasses, oldClasses, trainingDataStd,
            neuralNetFile, resultFile, typeSens, paramRecon,
            neuralNetworkModel)
        K.clear_session()
        
    else:
        raise ValueError('Your chosen backend is not supported! Unknown: ' + str(backend))



    
