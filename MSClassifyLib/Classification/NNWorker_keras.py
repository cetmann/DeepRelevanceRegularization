"""
This file contains all routines for training, predicting etc. with
the more modern 'keras' backend.
"""

import os, sys
import numpy as np
from tensorflow.keras import backend as K
import h5py

try:
    import cPickle
except ImportError:
    import pickle as cPickle

from NNClassifier_keras import NNClassifier_keras
from NNModel_keras import NNModel_keras

def startWork(createNewNetwork, neuralNetArchitecture, dataShape, numUniqueClasses,
    hyperparameters, mode, data, classes, args, newClasses, oldClasses, trainingDataStd,
    neuralNetFile, resultFile, typeSens, paramRecon,
    neuralNetworkModel):

    print('You are using the Keras backend for your experiments.')


    # Set the dimension ordering to 'channels first'
    if K.image_data_format() == 'channels_last':
        K.set_image_data_format('channels_first')
        print('Your Keras configuration was changed to use dim ordering: channels first ')



    ####################### Network Import ###########################

    if createNewNetwork:
        if neuralNetArchitecture == 'auto':
            # Load the standard NN architecture
            neuralNetArchitecture = 'NNStandardArchitecture_keras.py'

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
       
        # Build the network and print a summary 
        print('Data Shape:' + str(dataShape))
        netInputShape = dataShape[1:4]
        model = importedArchitecture.buildNet(netInputShape, numUniqueClasses)
        model.summary()

        # Create the classifier object
        neuralNetworkClassifier = NNClassifier_keras(architecture=model, 
            hyperparameter=hyperparameters, numUniqueClasses=numUniqueClasses,
            arguments=args)



    ######################### Training #########################
    
    if mode == 'train':
        print('Beginning neural network training...')
        neuralNetworkModel, train_history = neuralNetworkClassifier.trainModel(data, classes, neuralNetFile)

        # Save Neural Network and class dictionaries
        # TODO Replace pickle with .hdf5
        neuralNetworkModel.save(path=os.path.basename(neuralNetFile))
        neuralNetFileDump = {'newClasses' : newClasses,
                             'oldClasses' : oldClasses,
                             'numUniqueClasses' : numUniqueClasses,
                             'trainingDataStd' : trainingDataStd}

        fileObject = open(neuralNetFile, 'wb')
        cPickle.dump(neuralNetFileDump, fileObject)
        fileObject.close()


    ######################### Prediction #########################

    if mode == 'predict':
        # Define model
        neuralNetworkModel = NNModel_keras(network=neuralNetworkModel)

        # Perform predictions
        print('Create Prediction...')
        prediction = neuralNetworkModel.classify(data)
        prediction = np.array([oldClasses[classItem] for classItem
                               in prediction],
                           dtype=np.int)
        print('Prediction completed!')
        
        print('Get Logits...')
        logits = neuralNetworkModel.logits(data)
        print('Logits collected!')

        print('Get Scores...')
        scores = neuralNetworkModel.scores(data)
        print('Scores collected!')
        output = {'prediction' : prediction.T,\
                    'scores' : scores.T,\
                    'logits' : logits.T}
        #output = {'prediction' : prediction.T}

        # Save the prediction results
        if not os.path.exists(resultFile):
            os.mknod(resultFile)

        resultFileObject = h5py.File(resultFile,'w')
        print('Saving results...')
        for key in output.keys():
            resultFileObject[key] = output[key]
        resultFileObject.close()


    
    ######################### Sensitivity Analysis #########################
			
    if mode == 'sensitivity':
        neuralNetworkModel = NNModel_keras(neuralNetworkModel)
        # Select type of sensitivity analysis
        if typeSens == 'NLL':
            # Conversion of classes
            classesNew = np.array([newClasses[classItem] for classItem
                               in classes],
                           dtype=np.int)
            # Sensitivity analysis
            sensitivityBC = neuralNetworkModel.sensitivtyNLL(data, classesNew)
            print('NLL used')

        elif typeSens == 'perClass':
            print(data)
            sensitivityBC = neuralNetworkModel.sensitivityPerOutput(data,-1)
            print('perClass used')

        elif typeSens == 'perClassLogit':
            sensitivityBC = neuralNetworkModel.sensitivityPerOutput(data,-2)
            print('perClassLogit used')

        elif typeSens == 'binary': # default
            sensitivityBC = neuralNetworkModel.sensitivityBinaryCrossentropy(data)
            print('binary used')

        if sensitivityBC.__class__ == list:
            for listItem in sensitivityBC:
                # listItem is a one-element list
                # Transpose for Matlab compatibility
                listItem = listItem[0].T
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
        print('@reconstruction is not yet implemented in the Keras backend!')



    ################## Input Reconstruction in Basis ####################

    if mode == 'reconstructionBasis':
        print('@reconstructionBasis is not yet implemented in the Keras backend!')
        