"""
This file contains all routines for training, predicting etc. with
the classic 'theano' backend.
"""

# Import python packages
import theano
import numpy as np
import os, sys
import h5py
try:
    import cPickle
except ImportError:
    import pickle as cPickle

# Import self written python file(s)
from NNClassifier_theano import NNClassifier_theano

def startWork(createNewNetwork, neuralNetArchitecture, dataShape, numUniqueClasses,
    hyperparameters, mode, data, classes, args, newClasses, oldClasses, relu,
    trainingDataStd, neuralNetFile, resultFile, typeSens, paramRecon, basisRecon,
    neuralNetworkModel):

    print('You are using the theano backend for your experiments.')


    ####################### Network Import ###########################

    if createNewNetwork:
        if neuralNetArchitecture == 'auto':
            neuralNetArchitecture = 'NNStandardArchitecture_theano.py'

        # In order to be able to import the neuralNetArchitecture from a *.py
        # file, the directory the file is in needs to be added to the system
        # path.
        architectureDirectoryName = os.path.dirname(neuralNetArchitecture)
        sys.path.append(architectureDirectoryName)
        architectureFileName = os.path.basename(neuralNetArchitecture)
        architectureLibName = os.path.splitext(architectureFileName)[0]

        # Check whether there is an '__init__.py' file in
        # architectureDirectoryName. If not, create it. This is needed for
        # the import.
        if not os.path.exists((architectureDirectoryName + '/__init__.py')):
            open('__init__.py', 'a').close()

        # Import the defined network.
        import importlib
        importedArchitecture = importlib.import_module(architectureLibName)
        # First coordinate = None makes arbitrary batch sizes possible.
        netInputShape = (None,) + dataShape[1:4]
        layers = importedArchitecture.buildNet(netInputShape, numUniqueClasses)
        neuralNetworkClassifier = NNClassifier_theano(layers, hyperparameters)



    ######################### Training #########################

    if mode == 'train':
        # Sometimes, theano's compile lock is set.
        theano.gof.compilelock.set_lock_status(False)

        print('Beginning neural network training...')
        neuralNetworkModel = neuralNetworkClassifier.trainModel(data, classes)

        # load best model from training (if possible)
        if args.saveBestValidModel:
            try:
                neuralNetworkModel.network.load_params_from(\
                    neuralNetworkModel.network.bestModel)
                print('Best model is loaded')
            except Exception:
                print('Saving and restoring best model unsuccessful!')


        # Save Neural Network and class dictionaries
        neuralNetFileDump = {'neuralNetworkModel' : neuralNetworkModel,
                             'newClasses' : newClasses,
                             'oldClasses' : oldClasses,
                             'numUniqueClasses' : numUniqueClasses,
                             'relu' : relu,
                             'trainingDataStd' : trainingDataStd}

        fileObject = open(neuralNetFile, 'wb')
        cPickle.dump(neuralNetFileDump, fileObject)
        fileObject.close()


        # The following lines create a log file for the training process.
        NLL = list()
        l1 = list()
        l2 = list()
        logitSensitivity = list()
        logitSqSensitivity = list()
        probabilitySensitivity = list()
        lossSensitivity = list()

        for listitem in neuralNetworkClassifier.classifier.train_history_:
            NLL.append(listitem['cL'])
            l1.append(listitem['l1'])
            l2.append(listitem['l2'])
            logitSensitivity.append(listitem['loL'])
            logitSqSensitivity.append(listitem['loSqL'])
            probabilitySensitivity.append(listitem['prL'])
            lossSensitivity.append(listitem['lsL'])

        logDictionary = {'NLL' : NLL,
                         'l1' : l1,
                         'l2' : l2,
                         'logitSensitivity' : logitSensitivity,
                         'logitSqSensitivity' : logitSqSensitivity,
                         'probabilitySensitivity' : probabilitySensitivity,
                         'NLLSensitivity' : lossSensitivity}

        if args.useValidationSet:
            val_NLL = list()
            val_l1 = list()
            val_l2 = list()
            val_logitSensitivity = list()
            val_logitSqSensitivity = list()
            val_probabilitySensitivity = list()
            val_lossSensitivity = list()
            
            for listitem in neuralNetworkClassifier.classifier.train_history_:
                val_NLL.append(listitem['val_cL'])
                val_l1.append(listitem['val_l1'])
                val_l2.append(listitem['val_l2'])
                val_logitSensitivity.append(listitem['val_loL'])
                val_logitSqSensitivity.append(listitem['val_loSqL'])
                val_probabilitySensitivity.append(listitem['val_prL'])
                val_lossSensitivity.append(listitem['val_lsL'])
                
            logDictionary['val_NLL'] = val_NLL
            logDictionary['val_l1'] = val_l1
            logDictionary['val_l2'] = val_l2
            logDictionary['val_logitSensitivity'] = val_logitSensitivity
            logDictionary['val_logitSqSensitivity'] = val_logitSqSensitivity
            logDictionary['val_probabilitySensitivity'] = val_probabilitySensitivity
            logDictionary['val_lossSensitivity'] = val_lossSensitivity

#        if noiseLevellogging:
#            logDictionary['adaptiveNoiseLevel'] = adaptiveNoiseLevel

        print('Saving training log file...')
        
        # Use the filename of the neuralNetFile (minus possible file extension)
        # and add '_log.mat'
        logFileName = os.path.splitext(neuralNetFile)[0] + '_log.mat'
        trainLogFileObject = h5py.File(logFileName,'w')
        for key in logDictionary.keys():
            trainLogFileObject[key] = logDictionary[key]
        trainLogFileObject.close()



    ######################### Prediction #########################

    if mode == 'predict':
        prediction = neuralNetworkModel.classify(data)
        prediction = np.array([oldClasses[classItem] for classItem
                               in prediction],
                           dtype=np.int)
        logits = neuralNetworkModel.logits(data)
        output = {'prediction' : prediction.T,\
                    'scores' : neuralNetworkModel.scores(data).T,\
                    'logits' : logits.T}

        if not os.path.exists(resultFile):
            os.mknod(resultFile)

        resultFileObject = h5py.File(resultFile,'w')
        print('Saving results...')
        for key in output.keys():
            resultFileObject[key] = output[key]
        resultFileObject.close()


    
    ######################### Sensitivity Analysis #########################
				
    if mode == 'sensitivity':
        # select type of sensitivity analysis
        if typeSens == 'NLL':
            # conversion of classes
            classesNew = np.array([newClasses[classItem] for classItem
                               in classes],
                           dtype=np.int)
            # sensitivity analysis
            sensitivityBC = neuralNetworkModel.sensitivtyNLL(data, classesNew)
            print('NLL used')
        elif typeSens == 'perClass':
            sensitivityBC = neuralNetworkModel.sensitivityPerOutput(data,-1)
            print('perClass used')
        elif typeSens == 'perClassLogit':
            sensitivityBC = neuralNetworkModel.sensitivityPerOutput(data,-2)
            print('perClassLogit used')
        elif typeSens == 'binary': # default
            sensitivityBC = neuralNetworkModel.\
                sensitivityBinaryCrossentropy(data)
            print('binary used')

        if sensitivityBC.__class__ == list:
            for listItem in sensitivityBC:
                listItem = listItem.T
        else:
            sensitivityBC = sensitivityBC.T
        output = {'sensitivity': sensitivityBC}

        if not os.path.exists(resultFile):
            os.mknod(resultFile)

        resultFileObject = h5py.File(resultFile,'w')
        print('Saving results...')
        for key in output.keys():
            resultFileObject[key] = output[key]
        resultFileObject.close()
        


    ####################### Input Reconstruction ###########################

    if mode == 'reconstruction':
        # conversion of classes
        classesNew = np.array([newClasses[classItem] for classItem
                               in classes],
                           dtype=np.int)
        # start reconstruction
        reconInput = neuralNetworkModel.reconstructInput(
            data,
            classesNew,
            weighting=paramRecon['weighting'],
            fromLayer=paramRecon['reconLayer'],
            stepSize= paramRecon['reconStepSize'],
            numIterations=paramRecon['reconNumIterations'],
            predCoeff=paramRecon['predCoeff'],
            regCoeffL1=paramRecon['reconRegCoeffL1'],
            regCoeffL2=paramRecon['reconRegCoeffL2'])
        output = {'reconstruction': reconInput}

        if not os.path.exists(resultFile):
            os.mknod(resultFile)

        resultFileObject = h5py.File(resultFile,'w')
        print('Saving results...')
        for key in output.keys():
            resultFileObject[key] = output[key]
        resultFileObject.close()



    ################## Input Reconstruction in Basis ####################

    if mode == 'reconstructionBasis':
        # transposed due to hdf5
        basisRecon = basisRecon.T
        reconInput = neuralNetworkModel.reconstructInputinBasis(
            data,
            basisRecon,
            fromLayer=paramRecon['reconLayer'],
            stepSize= paramRecon['reconStepSize'],
            numIterations=paramRecon['reconNumIterations'], 
            regCoeffL1=paramRecon['reconRegCoeffL1'])
        output = {'reconstructionBasis': reconInput}

        if not os.path.exists(resultFile):
            os.mknod(resultFile)

        resultFileObject = h5py.File(resultFile,'w')
        print('Saving results...')
        for key in output.keys():
            resultFileObject[key] = output[key]
        resultFileObject.close()