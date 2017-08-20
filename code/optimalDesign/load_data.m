function [X, y] = load_data(file, orderAR, ctrlHzn)

load(['../data/' file '.mat']);

disturbances = [Ambient, Humidity];
control = [ClgSP, KitchenClgSP, GuestClgSP, SupplyAirSP, ChwSP];
output = TotalLoad;
proxy = [TOD, DOW];

X_d = [lagmatrix(disturbances, 0:1:orderAR-1), lagmatrix(output, 1), proxy];
X_c = control;
y = output;

X_d(1:orderAR-1,:) = [];
X_c(1:orderAR-1,:) = [];
y(1:orderAR-1,:) = [];
% X_d(end-ctrlHzn:end,:) = [];
% X_c(end-ctrlHzn:end,:) = [];
% y(end-ctrlHzn:end,:) = [];

X = [X_d, X_c];