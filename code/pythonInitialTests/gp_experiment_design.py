# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 19:09:30 2017

@author: Achin Jain (achinj@seas.upenn.edu)

"""

from __future__ import division

import numpy as np
from utils import load_data
import matplotlib.pyplot as plt
from plot_data import plot_experiment_design
import cPickle as pickle

building_type = 'LargeHotel'
n_samples = 2000

data_type = 'unconstrained'

load_model = False
#model_type = 'constrained' 
model_type = 'unconstrained' 
#model_type = 'rulebased'

# load model
if load_model:
    filename = ('../models/gp-%s-%s-numsamples-%s.pickle' %(model_type, building_type, n_samples))
    with open(filename, 'rb') as f:
        (model, buildingType, features_c, features_d, output, 
         normalizer_c, normalizer_d, normalizer_y) = pickle.load(f)
     

Xd_test, Xc_test, y_test = load_data(data_type, building_type, features_d, features_c, output)

# select a sample
idx = 3005

Xd_test = Xd_test.iloc[idx,:]
Xc_test = Xc_test.iloc[idx,:].as_matrix()
print(Xc_test)
y_test = y_test.iloc[idx,:].as_matrix()
Xd_test = normalizer_d.transform(Xd_test.reshape(1,-1))
#Xc_test = normalizer_c.transform(Xc_test.reshape(1,-1))


# vary control variables from min to max
ClgSP = np.linspace(22,32,50).reshape(-1,1)
KitchenClgSP = np.linspace(24,32,50).reshape(-1,1)
GuestClgSP = np.linspace(22,26,50).reshape(-1,1)
SupplyAirSP = np.linspace(12,14,50).reshape(-1,1)
ChwSP = np.linspace(3.7,9.7,50).reshape(-1,1)

#Xc_design = np.concatenate((ClgSP,KitchenClgSP,GuestClgSP,SupplyAirSP,ChwSP), axis=1)
Xc_design = np.concatenate((ClgSP, np.tile(Xc_test[1:], (50,1))), axis=1)
Xc_design = normalizer_c.transform(Xc_design)

Xd_design = np.tile(Xd_test, (50,1))

# predict using GP
X_design = np.concatenate((Xd_design, Xc_design), axis=1)
y_mean, y_std = model.predict(X_design, return_std=True)
y_mean = normalizer_y.inverse_transform(y_mean.reshape(-1, 1))
y_mean = y_mean.flatten()
y_std =y_std*(normalizer_y.data_max_-normalizer_y.data_min_)/2


plt.close('all')
plot_experiment_design(ClgSP.flatten(), y_mean, y_std, Xc_test[0], y_test[0])
#plot_experiment_design(KitchenClgSP.flatten(), y_mean, y_std)
#plot_experiment_design(GuestClgSP.flatten(), y_mean, y_std)
#plot_experiment_design(SupplyAirSP.flatten(), y_mean, y_std)
#plot_experiment_design(ChwSP.flatten(), y_mean, y_std)