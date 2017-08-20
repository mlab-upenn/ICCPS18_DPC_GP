# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 19:26:59 2017

@author: Achin Jain (achinj@seas.upenn.edu)

"""
    
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

# negative orderAR does forecast
def autoregress(x, orderAR):
    colName = x.name
    xInit = x
    x = pd.DataFrame()
    for ido in orderAR:
        temp = xInit
        temp = temp.shift(ido)
        if ido>0: temp.name = colName+'-%s' %ido
        if ido<0: temp.name = colName+'+%s' %(-ido)
        x = pd.concat([x, temp], axis=1)
    return x
    
# concatenates autoregresive/forecast terms in the features
def construct(features, data):
    x = pd.DataFrame()
    for colIdx, colName in enumerate(features.keys()):
        currentx = autoregress(data[colName], orderAR=features[colName])
        x = pd.concat([x, currentx], axis=1)
    return x

# normalize data between -1 and 1
def normalize(x):
    normalizer = MinMaxScaler(feature_range=(-1,1))
    x_normalized = normalizer.fit_transform(x)
    return x_normalized, normalizer
    
def load_data(dataType, buildingType, features_d, features_c, output, get_dummy=False):
    
    # training data
    allData = pd.read_csv('../data/' + dataType + '-' + buildingType  + '.csv')
    
    Xd = construct(features_d, data=allData)
    Xc = construct(features_c, data=allData)
    y = construct(output, data=allData)
    
    # remove NaNs
    maxAR = 0
    minAR = 0
    for idx, val in enumerate(features_d.values()):
        maxAR = np.max([maxAR,np.max(val)])
        minAR = np.min([minAR,np.min(val)])
    if minAR==0: minAR=-1
    Xd = Xd.iloc[maxAR:minAR]
    Xc = Xc.iloc[maxAR:minAR]
    y = y.iloc[maxAR:minAR]
    
    if get_dummy:
        tod_dummy = pd.get_dummies(Xd['TOD'])
        dow_dummy = pd.get_dummies(Xd['DOW'])
        Xd.drop('TOD',1)
        Xd.drop('DOW',1)
        Xd = pd.concat([Xd, tod_dummy, dow_dummy], axis=1)

    return Xd, Xc, y
    
if __name__ == '__main__':
    
    # example
    dataType = 'Train'
    buildingType = 'LargeHotel'
    features_d = {'TotalLoad':[1,2,3], 'Ambient':[0,1], 'Humidity':[0,1], 
    'TOD':[0],'DOW':[0]}
    features_c = {'ClgSP':[0,1]}
    output = {'TotalLoad':[0]}
    XdTrain, XcTrain, yTrain = load_data(dataType, buildingType, features_d, features_c, output)