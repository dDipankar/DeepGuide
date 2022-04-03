from collections import OrderedDict
import os
import sys
import warnings

import argparse
import logging
import h5py as h5
import numpy as np
import pandas as pd 
import scipy.io

import six
from six.moves import range

from sklearn.metrics import roc_auc_score, confusion_matrix
from keras.preprocessing import sequence
from keras.optimizers import RMSprop,Adam, SGD
from keras.models import Sequential, Model
from keras.layers.core import  Dropout, Activation, Flatten
from keras.regularizers import l1,l2,l1_l2
from keras.constraints import maxnorm
#from keras.layers.recurrent import LSTM, GRU
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.layers import Conv1D, MaxPooling1D, Dense, LSTM, Bidirectional, BatchNormalization, MaxPooling2D, AveragePooling1D, Input, Multiply, Add, UpSampling1D
from sklearn.metrics import mean_squared_error as mse
import scipy.stats as st

from Bio.Seq import Seq
import re
from random import shuffle
from Bio import SeqIO
np.random.seed(1369)

def detect_guides_cas9(fastaseq):

    ff = fastaseq
    seq = Seq(fastaseq)
    ff_reverse = seq.reverse_complement()
    ff_c = str(ff_reverse)

    substring = 'GG'

    guides_28 = []
    guides_20 = []
    guides_20_r = []

    # for forward strand
    matches = re.finditer(substring, ff)
    matches_positions = [match.start() for match in matches]

    for ii in matches_positions:
        if ii-23>=0:
            guide = ff[ii-23:ii+5]
            gg = ff[ii-21:ii-1]
            guides_28.append(guide)
            guides_20.append(gg)

    ind = []
    for g in guides_20:
        if ff.find(g)>=0:
            ind.append(ff.find(g))

    # for reverse complimentary strand
    matches = re.finditer(substring, ff_c)
    matches_positions = [match.start() for match in matches]
    for ii in matches_positions:
        if ii-23>=0:
            guide = ff_c[ii-23:ii+5]
            gg = ff_c[ii-21:ii-1]
            guides_28.append(guide)
            guides_20_r.append(gg)

    d = {'Guide_28': guides_28}
    df = pd.DataFrame(d)
    return df

def one_hot_encoding(lines):
    data_n = len(lines) 
    SEQ = np.zeros((data_n, 28, 4), dtype=int)
    
    for l in range(0, data_n):
        #data = lines[l].split(',')
        seq = lines[l]
        for i in range(28):
            if seq[i] in "Aa":
                SEQ[l, i, 0] = 1
            elif seq[i] in "Cc":
                SEQ[l, i, 1] = 1
            elif seq[i] in "Gg":
                SEQ[l, i, 2] = 1
            elif seq[i] in "Tt":
                SEQ[l, i, 3] = 1

    return SEQ

def scores_guides_cas9(guides):

    seq= one_hot_encoding(guides)

    # model architecture

    SEQ = Input(shape=(28,4))
    conv_1 = Conv1D(activation="relu", padding="same", strides=1, filters=20, kernel_size=5, kernel_regularizer = l2(0.0001))(SEQ)
    bat_norm1 = BatchNormalization()(conv_1)
    pool = MaxPooling1D(pool_size=(2))(bat_norm1)
    conv_2 = Conv1D(activation="relu", padding="same", strides=1, filters=40, kernel_size=8, kernel_regularizer = l2(0.0001))(pool)
    bat_norm2 = BatchNormalization()(conv_2)
    pool_1 = AveragePooling1D(pool_size=(2))(bat_norm2)
    enc = pool_1
    dec_pool_1 =  UpSampling1D(size=2)(enc)
    dec_bat_norm2 = BatchNormalization()(dec_pool_1)
    dec_conv_2  = Conv1D(activation="relu", padding="same", strides=1, filters=40, kernel_size=8, kernel_initializer='glorot_uniform',kernel_regularizer = l2(0.0001))(dec_bat_norm2)
    dec_pool = UpSampling1D(size=2)(dec_conv_2)
    dec_conv_1 = Conv1D(activation="relu", padding="same", strides=1, filters=20, kernel_size=5, kernel_initializer='glorot_uniform',kernel_regularizer = l2(0.0001))(dec_pool)
    dec = Conv1D(activation="relu", padding="same", strides=1, filters=4, kernel_size=5, kernel_initializer='glorot_uniform',kernel_regularizer = l2(0.0001))(dec_pool)
    model_seq = Model(inputs = SEQ, outputs= dec) 
    flatten = Flatten()(enc)
    dropout_1 = Dropout(0.5)(flatten)
    dense_1 = Dense(80, activation='relu', kernel_initializer='glorot_uniform')(dropout_1)
    dropout_2 = Dropout(0.5)(dense_1)
    out = Dense(units=1,  activation="linear")(dropout_2) 
    model = Model(inputs = SEQ, outputs= out)

    model.load_weights("../trained_deep_models/cas9_seq_wtt.h5")
    pred_y = model.predict(seq)

    score = [-1*ele for ele in pred_y.flatten()]
    gg = [x[2:22] for x in guides]
    d = {'Guide': gg, 'Score': score}
    df = pd.DataFrame(d)
    return df

fastafile = sys.argv[1]
fasta_sequences = SeqIO.parse(open(fastafile),'fasta')
if os.path.isfile('../data/activity_score_cas9_seq.csv'):
    os.remove('../data/activity_score_cas9_seq.csv')

for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        df = detect_guides_cas9(sequence)
        df = scores_guides_cas9(df['Guide_28'].tolist())
        df.to_csv('../data/activity_score_cas9_seq.csv', mode='a',index = None)
        print('Completed')
        print('see results in data/activity_score_cas9_seq.csv file')

