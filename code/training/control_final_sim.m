function [ sim_naive, sim_mc, sim_taylor, sim_taylor_sym ] = control_final_sim( gpmodel, valData, lag_y, lag_Ta, lag_humid, lag_dr, remove_nonworkdays, norm_y_min, norm_y_max )
% Performs various simulations for a given model.

% len = [1 96] + 95 + 96*8; % last number: 7 -> 7/7; 8 -> 8/7
len = [1 8] + 95 + 96*8; % last number: 7 -> 7/7; 8 -> 8/7

%% Construct a GP object for various simulation methods

% Check lag_y
lag_y = unique(floor(lag_y));
lag_output = numel(lag_y);
assert(max(lag_y) == lag_output && min(lag_y) == 1,...
    'lag_y must be [1:lag].');

mygp = simplegp.GP(gpmodel.hyp, gpmodel.mean, gpmodel.cov, gpmodel.lik, gpmodel.inputs, gpmodel.target);


%% Naive Multi-step simulation for 1 day
% Note that first day is Sunday
sim_naive = control_simulate(mygp, 'naive', valData, lag_y, lag_Ta, lag_humid, lag_dr, remove_nonworkdays, norm_y_min, norm_y_max, len);
t = [0:length(sim_naive.y)-1]' / 4; %time

f1=figure('Name', sprintf('Naive Simulation len = [%d, %d]', len(1), len(2)));
plotgp(f1, t, sim_naive.target, sim_naive.y, sqrt(sim_naive.s2));
xlim([0 24]);


%% MC Multi-step simulation for 1 day
% Note that first day is Sunday
sim_mc = control_simulate(mygp, 'mc', valData, lag_y, lag_Ta, lag_humid, lag_dr, remove_nonworkdays, norm_y_min, norm_y_max, len, 400);
t = [0:length(sim_mc.y)-1]' / 4; %time

f1=figure('Name', sprintf('MC Simulation len = [%d, %d]', len(1), len(2)));
plotgp(f1, t, sim_mc.target, sim_mc.y, sqrt(sim_mc.s2));
xlim([0 24]);


%% Taylor Multi-step simulation for 1 day
% Note that first day is Sunday
get_full_covariance_matrix = false;
sim_taylor = control_simulate(mygp, 'taylor', valData, lag_y, lag_Ta, lag_humid, lag_dr, remove_nonworkdays, norm_y_min, norm_y_max, len, get_full_covariance_matrix);
t = [0:length(sim_taylor.y)-1]' / 4; %time

f1=figure('Name', sprintf('Taylor Simulation len = [%d, %d]', len(1), len(2)));
plotgp(f1, t, sim_taylor.target, sim_taylor.y, sqrt(sim_taylor.s2));
xlim([0 24]);

%% Symbolic Taylor Multi-step simulation for 1 day
% Note that first day is Sunday
get_full_covariance_matrix = false;
import casadi.*
simlength = len(2) - len(1) + 1;
sym_u = casadi.MX.sym('u', simlength);
valData.u_dr = casadi.MX(valData.u_dr);
valData.u_dr(len(1):len(2)) = sym_u;
tic;
sim_taylor_sym = control_simulate(mygp, 'taylor', valData, lag_y, lag_Ta, lag_humid, lag_dr, remove_nonworkdays, norm_y_min, norm_y_max, len, get_full_covariance_matrix);
toc
sim_taylor_sym.sym_u = sym_u;

end
