#!/usr/bin/env python
'''
This file implements a custom class for mass spectrometry classification
with Keras. After training, a NNModel_keras instance is created.
'''

import numpy as np
import warnings


# NN packages
import tensorflow
from tensorflow.keras import losses
from tensorflow.keras import optimizers
from tensorflow.keras.optimizers import Adam as Adam
from tensorflow.keras.optimizers import SGD as SGD
from tensorflow.keras.regularizers import Regularizer
from tensorflow.keras import backend as K
from tensorflow.keras.utils import Sequence
from tensorflow.keras.callbacks import Callback, TensorBoard, CSVLogger
from tensorflow.keras.utils import to_categorical

# Custom output class
from NNModel_keras import NNModel_keras
from NNHelpers_keras import get_custom_loss_dict

# Libraries for data processing
import os
import sys 
from contextlib import redirect_stdout, redirect_stderr
from collections import OrderedDict
import time
import threading


K.set_image_data_format('channels_first')


class NoBatchTensorBoard(TensorBoard):
    def on_batch_end(self, batch, logs=None):
        pass

class NNClassifier_keras:
    """
    NNClassifier is a custom class for mass spectra classification. It is
    essentially a wrapper for Keras and processes the hyperparameters
    given by the neuralNetInterface.

    :param architecture: keras model object
    :param hyperparameter: dictionary of hyperparameter-value-pairs
    """

    def __init__(self, architecture, numUniqueClasses, hyperparameter = {}, arguments={}):
        self.model  = architecture
        self.args = arguments
        self.hyperp = hyperparameter
        self.numUniqueClasses = numUniqueClasses

        print('Number of unique classes: ' + str(self.numUniqueClasses))

        self.best_weights = []

        optimizer   = self.hyperp.setdefault('optimizer',Adam)
        lr          = self.hyperp.setdefault('learningRate',0.001)
        std_val     = self.hyperp.setdefault('trainingDataStd',None)

        ## Apply weight regularization to the model
        #self.weight_regularizer_dict = dict()
        #print('Apply regularizers...')
        #config = self.model.get_config()
        #layers = config.get('layers')
        #for layer in layers:
        #    # Get the type of the current layer
        #    class_name = layer.get('class_name')
        #    if 'Conv' in class_name or 'Dense' in class_name or 'Locally' in class_name or 'orm' in class_name:
        #        layer_name = layer.get('name')
        #        print('Added regularisation to: ' + layer_name)
        #        # Create a regularizer object and save it in a dict, so
        #        # we can access it later in the callback to get the
        #        # seperate loss values
        #        keras_l2_regularizer = tensorflow.keras.regularizers.l2
        #        regularizer = keras_l2_regularizer(self.hyperp.get('l2'))
        #        self.model.get_layer(layer_name).kernel_regularizer = regularizer
        #        self.model.get_layer(layer_name).bias_regularizer = regularizer
        #        if hasattr(self.model.get_layer(layer_name),'beta_regularizer'):
        #            self.model.get_layer(layer_name).beta_regularizer = regularizer
        #        self.weight_regularizer_dict[layer_name] = MaldiRegularizer(
        # #           l1=self.hyperp.get('l1'),
        #            l2=self.hyperp.get('l2'))

                #if hasattr(self.model.get_layer(layer_name),'gamma_regularizer'):
                #    self.model.get_layer(layer_name).gamma_regularizer = regularizer 


        # std_val may be a list with one entry
        if type(std_val) == list:
            std_val = std_val[0]

        # Create a custom loss objective
        self.custom_loss = MaldiLoss(model=self.model, 
            main_loss=self.hyperp.get('lossFunc'),
            l1=self.hyperp.get('l1'),
            l2=self.hyperp.get('l2'),
            logitSens=self.hyperp.get('logitSens'), 
            logitDiffSens=self.hyperp.get('logitDiffSens'),
            logitSqSens=self.hyperp.get('logitSqSens'), 
            probSens=self.hyperp.get('probSens'), 
            lossSens=self.hyperp.get('lossSens'), 
            std_val=std_val,
            verbose=self.hyperp.setdefault('verbose',1))

        # Create metric list from the custom loss object
        print('Added metrics')
        model_metrics = self.custom_loss.metric_list

        # Compile the model
        print('Compile model...')
        self.model.compile(loss=self.custom_loss.calc_loss, 
            optimizer=optimizer(lr=lr,epsilon=1e-8,amsgrad=True),
            metrics=model_metrics)


    def trainModel(self, data, classes, neuralNetFile):

        verbose = self.hyperp.setdefault('verbose',1.)        
        maxEpochs = self.hyperp.setdefault('epochs',100)
        batch_size = self.hyperp.setdefault('batchSize',128)
        n_classes = self.numUniqueClasses
        n_samples = data.shape[0]
        dim = data.shape[1:]

        # Define the validation generator
        validation_generator = None
        
        # Divide the data in train and validation
        validationSetRatio = self.hyperp.setdefault('validationSetRatio',.1)
        if validationSetRatio != 0:
            print(str(validationSetRatio) + ' of your data is for validation!')
            initialShuffleIndices = np.arange(len(data),dtype='int32')
            np.random.shuffle(initialShuffleIndices)
            data = data[initialShuffleIndices]
            classes = classes[initialShuffleIndices]
            
            numOfValPoints = np.floor(validationSetRatio * len(data))
            validationData = data[:np.int32(numOfValPoints)]
            validationClasses = classes[:np.int32(numOfValPoints)]

            # Define a generator for the validation data
            validation_generator = DataGenerator(data=validationData, 
                labels=validationClasses, dim=dim, batch_size=batch_size, 
                n_classes=n_classes, shuffle=True)
            
            data = data[np.int32(numOfValPoints):]
            classes = classes[np.int32(numOfValPoints):]

        classPrevalence = np.bincount(classes)
        self.classWeights = dict()
        for i in range(self.numUniqueClasses):
            if self.args.weightClasses:
                self.classWeights[i] = np.float(n_samples) / np.float(n_classes * classPrevalence[i])
            else:
                self.classWeights[i] = 1.

        print(self.classWeights)
        # Define the data generator here to include changes to the
        # data due to the possible train / val split
        training_generator = DataGenerator(data=data, labels=classes, 
            dim=dim, batch_size=batch_size, n_classes=n_classes, shuffle=True)
        
        
        # Define the callbacks
        # Our special Maldi Callback to monitor the loss and save the
        # best weights
        maldi_cb = MaldiCallback(Classifier_obj=self)

        # Tensorboard Callback
        nn_basename = os.path.split(neuralNetFile)[1]
        nn_basename = os.path.splitext(nn_basename)[0]
        logDir = os.path.dirname(neuralNetFile) + '/' + nn_basename + '_log'
        tboard_cb = NoBatchTensorBoard(log_dir=logDir, 
            histogram_freq=0, batch_size=batch_size, write_graph=False, 
            write_grads=False, write_images=False)


        # CSV-Logger Callback
        csv_name = logDir + '/' + nn_basename + '_csv_log'

        csv_logger = CSVLogger(filename=csv_name, separator=',', append=False)

        callbacks = [maldi_cb, tboard_cb, csv_logger]

        lr_schedule_epochs = self.hyperp.setdefault('lrSchedule', 0)
        if lr_schedule_epochs > 0:
            lr = self.hyperp.setdefault('learningRate', 0.001)
            def schedule(epoch, learning_rate):
                return 0.1 ** np.floor(epoch / lr_schedule_epochs) * lr

            lr_scheduler = tensorflow.keras.callbacks.LearningRateScheduler(
                schedule,
                verbose=np.int(verbose>0))
            callbacks.append(lr_scheduler)

        print('Beginning training loop...')    
        trainHistory = self.model.fit_generator(epochs=maxEpochs,
                                generator=training_generator,
                                validation_data=validation_generator,
                                use_multiprocessing=False,
                                workers=1,
                                max_queue_size=10,
                                verbose=np.int(verbose>0),
                                callbacks = callbacks,
                                class_weight = self.classWeights)


        # Load best weights into the model
        if self.args.saveBestValidModel:
            print("Applying early stopping...")
            self.model.set_weights(self.best_weights)

        return NNModel_keras(self.model), trainHistory



########################################################
############### The Keras Maldi Callback ###############
########################################################
class MaldiCallback(Callback):
    """
    Special callback to monitor and manipulate different 
    properties during the training of the model.
    You have access to the model and its parameters, as
    well as to the log of the loss

    Currently implemented:

    -   Storing the best weights according to the (val) loss
        value
    """

    def __init__(self, Classifier_obj, verbose=1):
        self.Classifier_obj = Classifier_obj
        self.verbose = verbose

        # Save the best weights according to the
        # smallest loss value
        if self.Classifier_obj.hyperp['validationSetRatio'] > 0.0 :
            self.monitor = 'val_loss'
        else :
            self.monitor = 'loss'

        self.monitor_op = np.less
        self.best = np.Inf


    def on_epoch_begin(self, epoch, logs=None):
        pass

    def on_epoch_end(self, epoch, logs=None):
        logs = logs or {}

        ############### Check for best weights ###############
        # Get the current loss-value
        current = logs.get(self.monitor)
        if current is None:
            warnings.warn('Can save best model only with %s available, '
                            'skipping.' % (self.monitor), RuntimeWarning)
        else:
            if self.Classifier_obj.args.saveBestValidModel:
                # Compare with the best result so far
                if self.monitor_op(current, self.best):
                    if self.verbose > 0:
                        print('\nEpoch %03d: %s improved from %0.5f to %0.5f,'
                                ' saving weights' %(epoch+1, self.monitor,
                                                    self.best, current))
                    self.best = current
                    self.Classifier_obj.best_weights = self.model.get_weights()
                else:
                    if self.verbose > 0:
                        print('\nEpoch %03d: %s did not improve from %0.5f' %
                                (epoch + 1, self.monitor, self.best))


    def on_batch_begin(self, batch, logs=None):
        pass


    def on_batch_end(self, batch, logs=None):
        pass


    def on_train_begin(self, logs=None):
        pass


    def on_train_end(self, logs=None):
        pass    



########################################################
############ The Keras Maldi Data Generator ############
########################################################
class DataGenerator(Sequence):
    """
    This class is a keras data generator that provides batches of
    data-label pairs during the execution.
    See below for possible data augmentation
    """

    def __init__(self, data, labels, dim, batch_size=128,
                 n_classes=10, shuffle=True):
        """
        Initialize the data generator.

        Parameter
        ----------
        data : numpy array
            Numpy array containing the data

        labels : numpy array
            Numpy array containing the labels in one hot format

        dim : int (tuple)
            Tuple of the data you provide to the model excluding 
            the batch dimension

        batch_size : int [std = 128]
            The batch size

        n_classes : int [std = 10]
            Number of different classes

        shuffle : bool [std = True]
            Shuffle the data after each epoch 
        """
        self.dim = dim
        self.batch_size = batch_size
        self.labels = labels
        self.data = data
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.on_epoch_end()

    def __len__(self):
        """
        Calculate the number of steps per epoch
        """
        len = int(np.floor(self.data.shape[0] / self.batch_size))
        return len

    def __getitem__(self, index):
        """
        Generate one batch of data
        """
        # Generate indexes of the batch
        indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

        # Generate data
        X, y = self.__data_generation(indexes)

        return X, y

    def on_epoch_end(self):
        """
        Shuffles indexes after each epoch
        """
        self.indexes = np.arange(self.data.shape[0])
        if self.shuffle == True:
            np.random.shuffle(self.indexes)

    def __data_generation(self, indexes):
        """
        Generates data containing batch_size samples

        Here you can implement data augmentation methods, that are
        applied during run time. BE CAREFUL: atm the data and labels
        are NOT copied via np.copy() and therefore you will directly 
        manipulate the original data. 
        """
        # Initialization
        X = np.empty((self.batch_size, *self.dim))
        y = np.empty((self.batch_size, self.n_classes))

        # Generate data
        for i, index in enumerate(indexes):
            # Store sample
            X[i] = self.data[index] 

            # Store class
            y[i] = to_categorical(self.labels[index], self.n_classes)

        return X, y



########################################################
########## The Keras Maldi Weight Regularizer ##########
########################################################
class MaldiRegularizer(Regularizer):
    """
    Implements the different WEIGHT regularizers for the network
    used in the Maldi context. The Parameter x in @__calc__
    can either be the weight, bias or activation of a layer.
    
    If you want to implement a regularisation method, that uses 
    label informations, do this in @MaldiLoss below.

    The following weight regularizers are currently fully implemented:
    L1, L2
    
    Arguments
    ----
        l1: Float 
            L1 regularization factor.

        l2: Float 
            L2 regularization factor.

    """

    def __init__(self, l1=0., l2=0.):
        self.l1 = K.cast_to_floatx(l1)
        self.l2 = K.cast_to_floatx(l2)

        # These are variable values, that can change during the
        # training process and act as a factor for the loss terms.
        # If you save a reference to the MaldiRegularizer-Object,
        # you can change them in a callback with K.set_value()
        self.l1_var = K.variable(l1)
        self.l2_var = K.variable(l2)

    def __call__(self, x):
        regularization = 0.0
        
        if self.l1:
            l1_loss = K.sum(self.l1_var * K.abs(x))
            regularization += l1_loss

        if self.l2:
            l2_loss = K.sum(self.l2_var * K.square(x))
            regularization += l2_loss

        return regularization

    def get_config(self):
        return {'l1': float(self.l1),
                'l2': float(self.l2)}



########################################################
############# The Keras Maldi Custom Loss ##############
########################################################
class MaldiLoss():
    """
    Implements all loss terms that rely on informations
    about the model output y_pred and the labels y_true
    
    Arguments
    ----
        model:  Keras model
            Supply the model object to access its input and output

        main_loss:  String
            This will be your main loss, e.g. categorical crossentropy
            We will first try to load a custom loss, then try to load 
            a keras standard loss and otherwise fall back to 
            categorical crossentropy

        logitSens: Float
            Logit Sensitivity regularization factor

        logitDiffSens: Float
            NOT YET IMPLEMENTED

        logitSqSens: Float
            Logit Squared Sensitivity regularization factor

        probSens: Float
            Probability Sensitivity regularization factor

        lossSens: Float
            Loss Sensitivity regularization factor

        std:

    """

    def __init__(self,
                 model,
                 main_loss,
                 l1=0.,
                 l2=0.,
                 logitSens=0.,
                 logitDiffSens=0.,
                 logitSqSens=0.,
                 probSens=0.,
                 lossSens=0.,
                 std_val=None,
                 verbose=1):
        
        # The main loss function
        self.main_loss = get_main_loss(main_loss)
        
        # The initialization values for the weights of the different
        # loss terms
        # These will decide if a loss term is active in the first place
        self.logitSens_on       = K.cast_to_floatx(logitSens)
        self.logitDiffSens_on   = K.cast_to_floatx(logitDiffSens)
        self.logitSqSens_on     = K.cast_to_floatx(logitSqSens)
        self.probSens_on        = K.cast_to_floatx(probSens)
        self.lossSens_on        = K.cast_to_floatx(lossSens)
        self.l1_loss_on         = K.cast_to_floatx(l1)
        self.l2_loss_on         = K.cast_to_floatx(l2)

        # These are variable values, that can change during the
        # training process and act as a factor for the loss terms.
        # If you save a reference to the MaldiLoss-Object, you can
        # change them in a callback with K.set_value()
        self.logitSens_var      = K.variable(logitSens, name='logitSens_var')
        self.logitDiffSens_var  = K.variable(logitDiffSens, name='logitDiffSens_var')
        self.logitSqSens_var    = K.variable(logitSqSens, name='logitSqSens_var')
        self.probSens_var       = K.variable(probSens, name='probSens_var')
        self.lossSens_var       = K.variable(lossSens, name='lossSens_var')

        # Same for L1/L2 loss
        self.l1_var = K.variable(l1, name='l1_var')
        self.l2_var = K.variable(l2, name='l2_var')

        # Save a reference to the model to access different layers and their input
        # and output        
        self.model = model

        self.verbose = verbose

        # Different methods for regularization standardization factors
        if isinstance(std_val, np.ndarray):
            self.std_val = K.cast_to_floatx(std_val)
            print("Scaling sensitivitites by the provided array.")
        elif std_val == 'model_input':
            self.std_val = model.input
            print("Scaling sensitivitites by input.")
        else:
            self.std_val = K.cast_to_floatx(1.0)
            print("The sensitivities are not scaled.")

        # Initialize some loss values
        self.main_loss_val = 0.0
        if self.logitSens_on:
            self.logitSens_loss = 0.0

        #if self.l1_loss_on:
        #    self.l1_loss = 0.0

        #if self.l2_loss_on:
        #    self.l2_loss = 0.0

        # Create a list of all active loss terms for the calculation of metrics
        self.metric_list = list()
        self.create_metric_list()


    def create_metric_list(self):
        self.metric_list.append(self.main_loss)

        self.metric_list.append(tensorflow.keras.metrics.categorical_accuracy)

        if self.verbose > 1:
            # Logit Sensitivity Loss
            if self.logitSens_on:
                self.metric_list.append(self.logitSens)

            # Logit Diff Sensitivity Loss
            if self.logitDiffSens_on:
                pass

            # Logit Squared Sensitivity Loss
            if self.logitSqSens_on:
                self.metric_list.append(self.logitSqSens)

            # Probability Sensitivity Loss
            if self.probSens_on:
                self.metric_list.append(self.probSens)

            # Main Loss Sensitivity Loss
            if self.lossSens_on:
                self.metric_list.append(self.lossSens)

            # Main L1 Weight Loss
            if self.l1_loss_on:
                self.metric_list.append(self.l1Loss)

            # Main L2 Weight Loss
            if self.l2_loss_on:
                self.metric_list.append(self.l2Loss)


    def calc_loss(self, y_true, y_pred):
        # y_true.shape = (batchsize, num_classes)

        # Initialize the loss value
        loss_val = 0.0

        # Calculate the main loss (e.g. categorical crossentropy)
        # This one will allways be calculated!
        # main_loss_val.shape = (batchsize,)
        self.main_loss_val = self.main_loss(y_true, y_pred)
        loss_val += self.main_loss_val

        # Logit Sensitivity Loss
        # Only used during training
        if self.logitSens_on:
            self.logitSens_loss  = self.logitSens(y_true, y_pred) 
            loss_val +=  self.logitSens_loss

        # Logit Diff Sensitivity Loss
        if self.logitDiffSens_on:
            pass

        # Logit Squared Sensitivity Loss
        # Only used during training
        if self.logitSqSens_on:
            logitSqSens_loss = self.logitSqSens(y_true, y_pred)
            loss_val += logitSqSens_loss

        # Probability Sensitivity Loss
        if self.probSens_on:
            probSens_loss = self.probSens(y_true, y_pred)
            loss_val += probSens_loss

        # Main Loss Sensitivity Loss
        if self.lossSens_on:
            lossSens_loss = self.lossSens(y_true, y_pred)
            loss_val += lossSens_loss

        if self.l1_loss_on:
            l1_loss_val = self.l1Loss(y_true, y_pred)
            loss_val+= l1_loss_val

        if self.l2_loss_on:
            l2_loss_val = self.l2Loss(y_true, y_pred)
            loss_val += l2_loss_val


        return loss_val



    def logitSens_main(self, y_true, y_pred):
        """
        Calculation of the logitSens_loss without 
        application of the weight factor
        """
        # Make sure that your logit layer is the
        # second last in your model! 
        # The output shape of this layer should be 
        # logitLayer.output.shape = (batchsize, num_classes)
        logitLayer = self.model.get_layer(index=-2)

        # Multiply logit output with the one hot label
        # logit.shape = (batchsize,)
        logit = K.sum(logitLayer.output * y_true, axis=-1)

        # Calculate the gradient
        # G_logit.shape = (batchsize, H, W, C)
        G_logit = K.gradients(K.sum(logit), self.model.input)[0]

        G_logit = self.std_val * G_logit


        return G_logit


    def logitSens(self, y_true, y_pred):
        # Scale with weight factor
        # Sparse saliency regularization
        # absG_logit.shape = (batchsize, H, W, C)
        G_logit = self.logitSens_main(y_true, y_pred)
        absG_logit = K.abs(G_logit)

        # absG_logit_flat.shape = (batchsize, H*W*C)
        absG_logit_flat = K.flatten(absG_logit)

        # logitSens_loss.shape = (batchsize,)
        logitSens_loss = K.sum(absG_logit_flat, axis=-1, keepdims=False)


        logitSens_loss =  self.logitSens_var * logitSens_loss
        return logitSens_loss


    def logitSqSens(self, y_true, y_pred):
        # Use the logitSens_loss if it was already computed
        G_logit = self.logitSens_main(y_true, y_pred)
        sqG_logit = K.square(G_logit)

        # absG_logit_flat.shape = (batchsize, H*W*C)
        sqG_logit_flat = K.flatten(sqG_logit)

        # logitSens_loss.shape = (batchsize,)
        logitSqSens_loss = K.sum(sqG_logit_flat, axis=-1, keepdims=False)
        logitSqSens_loss = self.logitSqSens_var * logitSqSens_loss
        return logitSqSens_loss


    def probSens(self, y_true, y_pred):
        # prob.shape = (batchsize,)
        prob = K.sum(self.model.output * y_true, axis=-1, keepdims=False)

        # G_prob.shape = (batchsize, H, W, C)
        G_prob = K.gradients(K.sum(prob, axis=-1), self.model.input)[0]

        G_prob = self.std_val * G_prob

        # absG_prob.shape = (batchsize, H, W, C)
        absG_prob = K.abs(G_prob)

        # absG_prob_flat.shape = (batchsize, H*W*C)
        absG_prob_flat = K.flatten(absG_prob)

        # probSens_loss.shape = (batchsize,)
        probSens_loss = K.sum(absG_prob_flat, axis=-1, keepdims=False)

        # Scale with weight factor
        probSens_loss = self.probSens_var * probSens_loss
        return probSens_loss


    def lossSens(self, y_true, y_pred):
        # Calculate the gradient of the main loss with respect to
        # the input of the model
        # G_loss.shape = (Batchsize, H, W, C)
        G_loss = K.gradients(self.main_loss_val, self.model.input)[0]
        
        G_loss = self.std_val * G_loss
        
        # We only want absolute values
        # absG_loss.shape = (batchsize, H, W, C)
        absG_loss = K.abs(G_loss)

        # absG_loss_flat.shape = (batchsize, H*W*C)
        absG_loss_flat = K.flatten(absG_loss)

        # lossSens_loss.shape = (batchsize,)
        lossSens_loss = K.sum(absG_loss_flat, axis=-1, keepdims=False)

        # Scale with weight factor
        lossSens_loss = self.lossSens_var * lossSens_loss

        return lossSens_loss

    def l1Loss(self, y_true, y_pred):
        l1Loss_val = 0.0
        regularizable_weights = [W for W in self.model.trainable_weights if 'kernel' in W.name or 'gamma' in W.name]
        for W in regularizable_weights:
            l1Loss_val+= self.l1_var * K.sum(K.abs(W))
        return l1Loss_val

    def l2Loss(self, y_true, y_pred):
        l2Loss_val = 0.0
        regularizable_weights = [W for W in self.model.trainable_weights if 'kernel' in W.name or 'gamma' in W.name]
        for W in regularizable_weights:
            l2Loss_val+= self.l2_var * K.sum(W**2)
        return l2Loss_val

    def get_config(self):
        return {'logitSens_on': float(self.logitSens_on),
                'logitDiffSens_on': float(self.logitDiffSens_on),
                'logitSqSens_on': float(self.logitSqSens_on),
                'probSens_on': float(self.probSens_on),
                'lossSens_on': float(self.lossSens_on),
                'logitSens_var': float(self.logitSens_var),
                'logitDiffSens_var': float(self.logitDiffSens_var),
                'logitSqSens_var': float(self.logitSqSens_var),
                'probSens_var': float(self.probSens_var),
                'lossSens_var': float(self.lossSens_var),
                'l1_var': float(self.l1_var),
                'l2_var': float(self.l2_var)}



def get_main_loss(main_loss_str):
    """
    Get the main loss function using a string.
    The search is structured in 3 steps:
        1. Search for the loss in the custom losses
        2. Search for the loss in the standard keras losses
        3. Fall back to categorical crossentropy if the loss does not exist

    Arguments
    ----
        main_loss_str:  String
            Name of the loss function you are looking
            for.

    """
    if type(main_loss_str) is list:
        main_loss_str = main_loss_str[0]

    custom_loss_dict = get_custom_loss_dict()
    loss = custom_loss_dict.get(main_loss_str, 'not_custom')

    if loss == 'not_custom':
        try:
            loss = losses.get(main_loss_str)
        except ValueError:
            print("We could not find your custom loss: " + str(main_loss_str))
            print("Fall back to categorical crossentropy")
            loss = losses.get('categorical_crossentropy')

    return loss
