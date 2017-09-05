function [X, y] = load_data(file, orderAR, ctrl_variables)

load(['../data/' file '.mat']);

disturbances = [Ambient, Humidity];

control = zeros(size(disturbances,1),numel(ctrl_variables));
for idc = 1:numel(ctrl_variables)
    control(:,idc) = eval(ctrl_variables{idc});
end
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