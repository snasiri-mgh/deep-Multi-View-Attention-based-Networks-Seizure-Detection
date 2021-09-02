
import scipy.io as sio
import numpy as np
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Dense, Activation, Permute, Dropout, Concatenate, Average, Reshape, Multiply
from tensorflow.keras.layers import Conv2D, MaxPooling2D, AveragePooling2D, AveragePooling1D, Conv1D, MaxPooling1D
from tensorflow.keras.layers import SeparableConv2D, DepthwiseConv2D
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import SpatialDropout2D
from tensorflow.keras.regularizers import l1_l2
from tensorflow.keras.layers import Input, Flatten
from tensorflow.keras.constraints import max_norm
from tensorflow.keras import backend as K
import os
import time
from sklearn.metrics import confusion_matrix

from keras.backend import expand_dims
from sklearn.utils.class_weight import compute_class_weight
from tensorflow.keras.constraints import NonNeg
from sklearn.metrics import classification_report
from sklearn import manifold
import matplotlib.pyplot as plt
from keras.models import model_from_json
from attention_ViewSelector import Attention_ViewSelector
from tensorflow.keras.utils import plot_model

#%%
path_folder = './'

path_folder_save = './'

#%%
for Out_Subject in np.arange(0, 1):
    
    print('\nSubject_out: ' + str(Out_Subject) + ' started...\n')
    
    
    View1_data = np.load(path_folder_save + 'View1_Subject'+str(Out_Subject)+'.npy')
    View2_data = np.load(path_folder_save + 'View2_Subject'+str(Out_Subject)+'.npy')
    View3_data = np.load(path_folder_save + 'View3_Subject'+str(Out_Subject)+'.npy')
    labels_test = np.load(path_folder_save + 'Subject'+str(Out_Subject)+'_epochLabels.npy')

    
#%% EEGNET model #########################################3###############################################################################################################################:
    nb_classes = 2
    Chans = 22
    Samples = 750
    dropoutRate = 0.5
    kernLength = 32
    F1 = 4
    D = 1
    F2 = F1*D 
    dropoutType = 'Dropout'
    norm_rate = 0.25
    Concat_Dense = 16
        
    input1   = Input(shape = (Chans, Samples, 1))

    block1       = Conv2D(F1, (1, kernLength), padding = 'same',
                                   input_shape = (Chans, Samples, 1),
                                   use_bias = False)(input1)
    block1       = BatchNormalization(axis = 1)(block1)
    block1       = DepthwiseConv2D((Chans, 1), use_bias = False, 
                                   depth_multiplier = D,
                                   depthwise_constraint = max_norm(1.))(block1)
    block1       = BatchNormalization(axis = 1)(block1)
    block1       = Activation('elu')(block1)
    block1       = AveragePooling2D((1, 4))(block1)
    block1       = Dropout(dropoutRate)(block1)
        
    block2       = SeparableConv2D(F2, (1, 16), use_bias = False, padding = 'same')(block1)
    block2       = BatchNormalization(axis = 1)(block2)
    block2       = Activation('elu')(block2)
    block2       = AveragePooling2D((1, 4))(block2)
    block2       = Dropout(dropoutRate)(block2)
        
    flatten1      = Flatten()(block2)
    flatten1        = Dense(Concat_Dense)(flatten1)
    flatten1_temp = Reshape((1, flatten1.shape[-1]))(flatten1)
    
    dense_model1        = Dense(2, kernel_constraint = max_norm(norm_rate))(flatten1)
    softmax1      = Activation('softmax')(dense_model1)

    model1 = Model(inputs=input1, outputs=softmax1)
    
    
# CNN-image model #######################################################################################################3333###############################################################3:
        
    input2   = Input(shape = (View2_data.shape[1], View2_data.shape[2], View2_data.shape[3]))

    block1       = BatchNormalization(axis = 1)(input2)
    block1       = Conv2D(32, (7, 7), padding = 'valid',
                                   input_shape = (View2_data.shape[1], View2_data.shape[2], View2_data.shape[3]))(block1)    
    block1       = BatchNormalization(axis = 1)(block1)
    block1       = Activation('elu')(block1)
    block1       = Conv2D(32, (7, 7), padding = 'valid')(block1)
    block1       = BatchNormalization(axis = 1)(block1)
    block1       = Activation('elu')(block1)
    block1       = MaxPooling2D((2, 2))(block1)
    block1       = Dropout(dropoutRate)(block1)
    
    block2       = Conv2D(32, (5, 5), padding = 'valid')(block1)
    block2       = BatchNormalization(axis = 1)(block2)
    block2       = Activation('elu')(block2)
    block2       = Conv2D(Concat_Dense, (5, 5), padding = 'valid')(block2)
    block2       = BatchNormalization(axis = 1)(block2)
    block2       = Activation('elu')(block2)
    block2       = MaxPooling2D((2, 2))(block2)
    block2       = Dropout(dropoutRate)(block2)
        

    flatten2      = Flatten()(block2)
    flatten2_temp = Reshape((1, flatten2.shape[-1]))(flatten2)
    
    dense_model2        = Dense(2, kernel_constraint = max_norm(norm_rate))(flatten2)
    softmax2      = Activation('softmax')(dense_model2)

    model2 = Model(inputs=input2, outputs=softmax2)
        
## MLP model #######################################3333##############################################3###################################################################################3:
    nb_classes = 2
    Chans = 22
    Samples = 946
    dropoutRate = 0.5
    F2 = 128
    F3 = 64
    F4 = 64
    AP = 8
    dropoutType = 'Dropout'
    norm_rate = 0.1
        
        
    input3   = Input(shape = (Samples,))

    block1       = BatchNormalization(axis = 1)(input3)
    block1       = Dense(Concat_Dense, input_shape = (Chans,))(block1)
    block1       = BatchNormalization(axis = 1)(block1)
    block1       = Activation('elu')(block1)
    block1       = Dropout(dropoutRate)(block1)
    

    flatten3      = Flatten()(block1)
    flatten3_temp = Reshape((1, flatten3.shape[-1]))(flatten3)
    
    dense_model3        = Dense(2, kernel_constraint = max_norm(norm_rate))(flatten3)
    softmax3      = Activation('softmax')(dense_model3)

    model3 = Model(inputs=input3, outputs=softmax3)
    
# Merged Models ##################################################################################################################  ######################################################3 

    dense = Concatenate(axis=1)([flatten1_temp, flatten2_temp, flatten3_temp])  
    flatten4_1 = Attention_ViewSelector()(dense)
    
    dense = Concatenate(axis=1)([flatten1_temp, flatten3_temp, flatten2_temp])  
    flatten4_2 = Attention_ViewSelector()(dense)

    dense = Concatenate(axis=1)([flatten2_temp, flatten1_temp, flatten3_temp])  
    flatten4_3 = Attention_ViewSelector()(dense)

    dense = Concatenate(axis=1)([flatten2_temp, flatten3_temp, flatten1_temp])  
    flatten4_4 = Attention_ViewSelector()(dense)

    dense = Concatenate(axis=1)([flatten3_temp, flatten1_temp, flatten2_temp])  
    flatten4_5 = Attention_ViewSelector()(dense)

    dense = Concatenate(axis=1)([flatten3_temp, flatten2_temp, flatten1_temp])  
    flatten4_6 = Attention_ViewSelector()(dense)
    
    dense = Concatenate()([flatten4_1, flatten4_2, flatten4_3, flatten4_4, flatten4_5, flatten4_6])
    dense       = Dropout(dropoutRate)(dense)
    flatten4      = Flatten()(dense)

    dense        = Dense(2, kernel_constraint = max_norm(norm_rate))(flatten4)
    softmax      = Activation('softmax')(dense)
    model = Model(inputs=[input1, input2, input3], outputs=[softmax])
    
#    model.summary()
#%%#####################################################################################################################################################################################33
    model.compile(loss=['categorical_crossentropy'], optimizer='adam', 
              metrics = ['accuracy'])    
    
    model.load_weights(path_folder+'model_weights_sub'+str(Out_Subject)+'.h5')
#%% Temporal and saptial weights:
    W = model.get_weights()
    Temporal_weights = W[10]
    Spatial_weights = W[14]
    
    np.save(path_folder_save+'Temporal_weights_sub'+str(Out_Subject)+'.npy', Temporal_weights)
    np.save(path_folder_save+'Spatial_weights_sub'+str(Out_Subject)+'.npy', Spatial_weights)
#%%
    View3_data = np.squeeze(View3_data)
    labels_test_predict = model.predict([View1_data, View2_data, View3_data], batch_size=16)            
    y_true = np.argmax(labels_test, 1)
    y_pred =  np.argmax(labels_test_predict, 1)    
    y_pred_scores = np.max(labels_test_predict, 1) 
    C = confusion_matrix(y_true, y_pred, labels=[0,1])    
    print('\nSubject_out: ' + str(Out_Subject) +' >>>>>>>>>>>> Model Org: \n'+ str(C))        
    print(classification_report(y_true, y_pred))
    
