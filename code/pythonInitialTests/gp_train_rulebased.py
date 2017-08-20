# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 16:49:34 2017

@author: Achin Jain (achinj@seas.upenn.edu)

"""

from __future__ import division

print('======================================================================')

print('importing libraries...')
import numpy as np
from collections import OrderedDict
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.metrics import explained_variance_score
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (RBF, RationalQuadratic,
                                              ConstantKernel, Sum, Product)
from utils import load_data, normalize
from plot_data import plot_gp
import time
import cPickle as pickle

print('----------------------------------------------------------------------')
print('simulation params: \n')

buildingType = 'LargeHotel'
print('building type: %s' %buildingType)

orderAR = 3
ctrlHzn = 1
print('order of AR: %s' %orderAR)
print('control horizon: %s' %ctrlHzn)

n_samples = 2000
n_test_samples = 2000

plot_results = False
save_results = True

print('----------------------------------------------------------------------')
print('loading data...')

distRange = list(range(-ctrlHzn+1,orderAR))
controlRange = list(range(-ctrlHzn+1,1))

features_d = OrderedDict([('Ambient',distRange), ('Humidity',distRange), 
              ('TOD',[0]), ('DOW',[0]), ('HtgSP',[0]), ('KitchenHtgSP',[0]),
              ('TotalLoad',range(1,orderAR))])
features_c = OrderedDict([('ClgSP',controlRange), ('KitchenClgSP',controlRange), 
              ('GuestClgSP',controlRange), ('SupplyAirSP',controlRange),
              ('ChwSP',controlRange)])
output = OrderedDict([('TotalLoad',[0])])

XdTrain, XcTrain, yTrain = load_data('rulebased', buildingType, features_d, features_c, output)
XdTest, XcTest, yTest = load_data('unconstrained', buildingType, features_d, features_c, output)

# select samples
XdTrain = XdTrain.iloc[0:n_samples,:]
XcTrain = XcTrain.iloc[0:n_samples,:]
yTrain = yTrain.iloc[0:n_samples,:]
n_features = XdTrain.shape[1]

XdTest = XdTest.iloc[n_test_samples::,:]
XcTest = XcTest.iloc[n_test_samples::,:]
yTest = yTest.iloc[n_test_samples::,:]

# normalize data
XdTrain, normalizer_d = normalize(XdTrain)
XdTest = normalizer_d.transform(XdTest)
XcTrain, normalizer_c = normalize(XcTrain)
XcTest = normalizer_c.transform(XcTest)
yTrain, normalizer_y = normalize(yTrain)
                             
print('----------------------------------------------------------------------')

print('training a Gaussian Process on selected features... \n')
print('using %s samples for training' %n_samples)

XTrain = np.concatenate((XdTrain,XcTrain), axis=1)
k1 = 1.0**2*RBF(length_scale=np.ones(XTrain.shape[1]), 
                length_scale_bounds=(1e-7, 1e7))
k2 = ConstantKernel(0.1, (1e-7, 1e7))
k3 = 1.0**2*RationalQuadratic(length_scale=1.0, 
                              length_scale_bounds=(1e-7, 1e7), alpha=0.1)
kernel = Product(Sum(k1,k2),k3)

# generate data and fit GP
model = GaussianProcessRegressor(kernel=kernel, normalize_y=True)
start = time.time()
model.fit(XTrain, yTrain.reshape(-1, 1))
end = time.time()
print('elapsed time: %ss \n' %(end - start))

print('final kernel: %s' %(model.kernel_))

print('----------------------------------------------------------------------')

print('calculating posterior on test data... \n')
start = time.time()
XTest = np.concatenate((XdTest,XcTest), axis=1)
ymu, ystd = model.predict(XTest, return_std=True)

# rescaling mu and sigma
ymu = normalizer_y.inverse_transform(ymu.reshape(-1, 1))
ymu = ymu.flatten()
ystd =ystd*(normalizer_y.data_max_-normalizer_y.data_min_)/2

end = time.time()
print('elapsed time: %ss' %(end - start))

print('----------------------------------------------------------------------')

print('results: \n')
ytrue = np.array(yTest).flatten()
ypred = ymu
MSE = mean_squared_error(ytrue, ypred, multioutput='raw_values')
R2Score = r2_score(ytrue, ypred, multioutput='raw_values')
EV = explained_variance_score(ytrue, ypred, multioutput='raw_values')

print('root mean square error: %s' %(np.sqrt(MSE)))
print('normalized mean square error: %s' %(np.sqrt(MSE)/np.array(np.mean(ytrue))))
print('R2 score: %s' %(R2Score))
print('explained variance: %s' %(EV))

print('----------------------------------------------------------------------')
#plt.close('all')

if plot_results:
    print('plotting results...')
    plot_gp(ytrue, ymu, ystd)
    
else:
    print('results not plotted.')

print('----------------------------------------------------------------------')

if save_results:
    print('saving results...')    
    filename = ('../models/gp-rulebased-%s-numsamples-%s.pickle' %(buildingType, n_samples))
    
    with open(filename, 'wb') as f:
        pickle.dump([model, buildingType, features_c, features_d, output,
                     normalizer_c, normalizer_d, normalizer_y], f)

else:
    print('results not saved.')

print('======================================================================')