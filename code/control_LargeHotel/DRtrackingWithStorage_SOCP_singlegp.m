% Demand tracking controller using GP, with a coupled battery / storage
% (tracking) slack variables are used.
%
% This code uses naive simulation of a single GP model.
%
% The memory requirement may be large because multiple GP objects are used
% simultaneously.
%
% Requires:
%   CasADi v3.0 or above.
%   CasADi optistack.
%   IPOPT with MA27 or other MA* solver installed.

% This approximates the two-sided chance constraints with an SOCP derived
% from this paper (lemma 16): http://www.optimization-online.org/DB_FILE/2016/02/5345.pdf
%
% The optimization variables are:
%   S  : DR signal
%   Pb : controlled battery power when not tracking
%   Delta: tracking slack delta(i) when tracking
%   P  : expected power
%   Sigma: variance (sigma, not sigma^2) of power
%   SOC: expected battery's SoC at each time step
%   this is not needed -> tP  : The t variables used in the above lemma for Pb chance constraints
%   tS  : The t variables used in the above lemma for SOC chance constraints

classdef DRtrackingWithStorage_SOCP_singlegp < OBNNode
    properties (Access=public)
        % The optimization problem using casadi/optistack
        % optvar = structure of optimization variables
        % optpar = structure of parameters
        % optprob = the NLP object, used to solve
        %optvar = struct();
        optpar = struct();
        % optprob
        
        % Parameters for the cost function that takes into account the
        % confidence in demand <= reference - beta*variance
        alpha
        beta
        
        timestep = 0.25;    % time step in hours
        gpmodels     % The GP models
        gpcommon        % Common values for GPs
        ref
        reflen
        curRef      % vector of current reference values (unnormalized), excluding NaN values
        curOutputIdx   % indices of predicted output values used in computing the cost (c.t. curRef above), i.e. indices of non-NaN ref values
        
        active_ahead
        horizon
        
        % Weights in cost function
        wvar
        wdelta
        wPb
        
        ramplimit = inf;    % Ramp limit on the CHW setpoint (normalized)
        
        prevSol = [];  % previous solution, used to warm start
        curHorizon      % the current horizon of the optimization
        
        %{
        ar_ta = [];
        ar_humid = [];
        ar_dr = [];
        ar_power = [];
        %}
        
        curStep = 1;    % the current step, used to get the reference value
        unconstrained_solved = 0;   % number of times the unconstrained NLP gives feasible solutions
        constrained_solved = 0;     % number of times the constrained NLP need to be solved

        batt_params     % battery parameters
        batt_soc        % SoC of battery (kWh)
        batt_used       % Whether to use the battery

        cur_batt_power = 0;  % current battery power, NaN if it takes the residual between ref and building's power
        cur_delta = NaN;      % If cur_batt_power is NaN, this is the tracking slack

        batt_power_tracking     % indices within the horizon that we track the reference -> battery should only take the residual
        batt_power_nontracking  % indices within the horizon that we don't track the reference -> battery power can be optimized

        % SignalsValues object to store all time-series variables in the
        % system.
        signals;
        
        last_active;        % whether the last step was active (to update the GPs)
        
        % These vectors keep track of past values of various variables used
        % to update the GPs.
        last_hour;
        last_dayofweek;
        
        % Thresholds for deciding whether to include a new point
        % If a new point either has error in mean prediction >= .mean, or a
        % standard deviation (sqrt of variance) >= .stdvar, the point will
        % be included. All values are denormalized (i.e., in kW).
        egp_thresholds = struct('mean', 0, 'stdvar', 0);
        
        %{
        dr_signal        % to save the DR signal time series
        power_expected  % expected power
        power_var       % variance of the prediction
        batt_power_hist      % save the battery power history
        batt_soc_hist      % save the battery SoC history
        tracking_delta      % save the history of tracking delta
        forecasts   % Forecasts of Ta and humid
        %}
    end
    
    methods (Access=public)
        function obj = DRtrackingWithStorage_SOCP_singlegp(gpmodels, normparams, ref, forecasts, active_ahead, horizon, ...
                wdelta, wPb, wvar,...
                ramplimit, ...
                batt_params, use_battery, ...
                name, varargin)
            % gpmodel: structure array of a GP model of fields:
            % covariance_function, mean_function, likelihood, model_inputs,
            % model_target, model_excepts, hyp, inputs, target (last two
            % are for the normalized training set).
            %
            % normparams : normalization parameters.
            %
            % ref: vector of reference values. An entry can be NaN, which
            % means this entry should be ignored in the tracking.
            %
            % forecasts: a structure containing forecasts for Ambient and
            % Humidity, which are vectors. Note that time is shifted by 1,
            % i.e. at time 1, we need forecasts at 2, 3, 4, ... which are
            % obtained from forecasts(1,2,3,...); so forecasts start from
            % the same index as the current time step.
            %
            % active_ahead is the number of steps the controller will look
            % ahead to see whether it should be activated. At time step k,
            % the controller will look active_ahead steps into the future:
            % if it sees any non-NaN reference values, the controller is
            % activated; otherwise, it's inactive.
            %
            % horizon >= active_ahead is the horizon of MPC.
            %
            % wdelta, wPb, wvar are the weights in the objective function
            % (see the objective function method for details).
            %
            % batt_params is a structure of the battery's parameters, inc.
            %   .power_weight = specifies probability of satisfying
            %       power constraint.
            %   .soc_weight = specifies probability of satisfying
            %       SoC constraint.
            %   .power_max = maximum power kW (for charging, positive).
            %   .power_min = minimum power kW (for discharging, negative).
            %   .soc_max, soc_min = max and min of SoC in kWh
            %   .timestep = time step for integrating the power, in hours.
            %
            % The objective function is:
            %   sum_i { (r_k - mu_k)^2 + wvar*var(y_k) }
            % where mu_k = E[y_k]; and i is for all r_k that is not NaN.
            
            assert(isstruct(normparams));
            assert(isvector(ref), 'ref must be a vector of reference values.');
            assert(isvector(forecasts.Ambient) && length(forecasts.Ambient) >= length(ref) && ...
                isvector(forecasts.Humidity) && length(forecasts.Humidity) >= length(ref),...
                'Forecasts must contain vectors similar to and at least as long as ref.');
            assert(active_ahead > 0 && active_ahead <= horizon, 'Must have 0 < active_ahead <= horizon.');
            
            assert(isscalar(wdelta) && isnumeric(wdelta) && isscalar(wPb) && isnumeric(wPb));
            
            assert(~isfinite(ramplimit) || ramplimit > 0, 'Invalid RAMPLIMIT.');
            
            % check batt_params
            assert(isstruct(batt_params) && ...
                all(isfield(batt_params, {'power_weight', 'soc_weight', 'power_max', 'power_min', 'soc_max', 'soc_min', 'timestep'})),...
                'Invalid battery parameter structure.');
            assert(batt_params.power_max > batt_params.power_min && ...
                batt_params.soc_max > batt_params.soc_min && ...
                batt_params.timestep > 0);
            
            % Process GP model
            assert(isstruct(gpmodels));
            
            gpmodel_commons = struct();
            mygpmodels = struct();
            
            % create the signal sources and add some signals
            signals = SignalsValues();
            
            % Note that for forecasts (e.g., fTa), fTa(k) is the first
            % forecast value we will need at step k, i.e., fTa(k) is
            % Ta(k+1) in the signals below. That's why we prepend the
            % forecast vectors with 0.
            signals.addMultiSignals(struct(...
                'Ambient', [0; forecasts.Ambient], ...
                'Humidity', [0; forecasts.Humidity], ...
                'TotalLoad', [], ...
                'TOD', [], ...
                'DOW', [], ...
                'ChwSP', [], ...
                'ClgSP', [], ...
                'GuestClgSP', [], ...
                'KitchenClgSP', [], ...
                'SupplyAirSP', []), ...
                normparams);

            signals.addSignal('power_expected');    % expected power
            signals.addSignal('power_var');         % variance of the prediction
            signals.addSignal('batt_power');        % save the battery power history
            signals.addSignal('batt_soc');          % save the battery SoC history
            signals.addSignal('tracking_delta');
            
            fprintf('Loading and pre-computing symbolic functions with CasADi for GP...\n');

            % Create the signal model for this GP
            mygpmodels.signalmodel = SignalsModel(signals);
            mygpmodels.signalmodel.setInputs(gpmodels.model_inputs);
            mygpmodels.signalmodel.setOutput(gpmodels.model_target);
            
            mygpmodels.gp = nextgp.GP(nextgp.GPData(gpmodels.hyp, gpmodels.mean_function, ...
                gpmodels.covariance_function, gpmodels.likelihood, ...
                gpmodels.inputs, gpmodels.target));
            
            % Functions to denormalize the output mean and variance
            gpmodel_commons.norm_y_min = normparams.(gpmodels.model_target).min;
            gpmodel_commons.norm_y_max = normparams.(gpmodels.model_target).max;
            DeltaYnorm = gpmodel_commons.norm_y_max - gpmodel_commons.norm_y_min;
            gpmodel_commons.norm_y_delta = DeltaYnorm/2;

            % DONE
            disp('Done pre-computing symbolic functions with CasADi for optimization.');
            
            obj = obj@OBNNode(name, varargin{:});
            
            % Assign properties
            obj.gpmodels = mygpmodels;
            obj.gpcommon = gpmodel_commons;
            
            obj.signals = signals;
            
            obj.ref = ref(:);
            obj.reflen = numel(ref);
            obj.active_ahead = active_ahead;
            obj.horizon = horizon;
            
            obj.wdelta = wdelta;
            obj.wPb = wPb;
            
            % Calculate the normalized ramp limit
            % -r/Delta <= x' - x <= r/Delta where Delta = (max-min)/2
            if isfinite(ramplimit)
                obj.ramplimit = ramplimit * 2 / (normparams.ChwSP.max - normparams.ChwSP.min);
            else
                obj.ramplimit = inf;
            end
            
            obj.batt_used = use_battery;
            
            obj.batt_params = batt_params;
            obj.batt_soc = (batt_params.soc_min + batt_params.soc_max)/2; % initialize SoC
            obj.batt_params.norm_power_delta = (obj.batt_params.power_max - obj.batt_params.power_min)/2;
            % The normalized Pb value corresponding to actual Pb = 0
            obj.batt_params.Pb_zero = -obj.batt_params.power_min/obj.batt_params.norm_power_delta - 1;
            % The normalized SoC bounds
            obj.batt_params.soc_max_norm = obj.batt_params.soc_max/obj.batt_params.timestep/gpmodel_commons.norm_y_delta;
            obj.batt_params.soc_min_norm = obj.batt_params.soc_min/obj.batt_params.timestep/gpmodel_commons.norm_y_delta;
            
            
            % Pre-Calculate the Phi-inverse values for chance constraints
            eps_power = (1 - batt_params.power_weight)/1.25;
            eps_soc = (1 - batt_params.soc_weight)/1.25;
            obj.optpar = struct(...
                'PhiEps', norminv(eps_power), ...
                'PhiEpsCompl', norminv(1-eps_power), ...
                'PhiHalf', norminv(eps_power/2), ...
                'PhiEpsSOC', norminv(eps_soc), ...
                'PhiEpsComplSOC', norminv(1-eps_soc), ...
                'PhiHalfSOC', norminv(eps_soc/2));
            
            % save the initial SoC
            obj.signals.setSignal('batt_soc', 1, obj.batt_soc);
            
            obj.last_active = false;
            nGPs = 1;
            obj.last_hour = zeros(1, nGPs);
            obj.last_dayofweek = zeros(1, nGPs);
            
            % Create the input, output, callback...
            obj.createOutputPort('setpoints', 'v', 'double');
            
            % Create input ports
            obj.createInputPort('bldgin', 'v', 'double');  % Ambient temperature
            
            obj.addCallback(@obj.calculateDRsignal, 'Y', 0);
            %obj.addCallback(@obj.initializeSim, 'INIT');
            obj.addCallback(@obj.finalizeSim, 'TERM');
        end
        
        function calculateDRsignal(self)            
            if mod(self.curStep, 10) == 0
                fprintf('Step %d ...\n', self.curStep);
            end
            
            bldgInputs = self.input('bldgin');
            
            power = bldgInputs(9);
            
            % Update signals with measurements
            if self.curStep > 1
                self.signals.setSignal('TotalLoad', self.curStep-1, power);
            end
            self.signals.setSignal('Ambient', self.curStep, bldgInputs(6));
            self.signals.setSignal('Humidity', self.curStep, bldgInputs(7));

            % NOTE: except for the time of day (in hours, (0, 24]), all
            % date values are shifted by one step (delayed by 1 step), so
            % at 24:00 = 0:00 next day, the values of month, day-of-month,
            % day-of-week and holiday are of the previous day.
            % Typically this does not affect the algorithm, if there is no
            % change in the output of the algorithm at 0:00.
            tod = bldgInputs(3);   % time of day in hour
            dayofweek = bldgInputs(4);
            self.signals.setSignal('TOD', self.curStep, tod);
            self.signals.setSignal('DOW', self.curStep, dayofweek);
            
            % Update the battery status NOW, before calling isActive()
            % NOTE that the power value is of the previous step, so we
            % don't update power and SoC until step 2.
            if self.curStep > 1
                % curRef = last reference value
                % cur_batt_power = NaN if battery should take the residual,
                %       otherwise it's fixed power
                % batt_soc = last battery's SoC
                if isnan(self.cur_batt_power)
                    self.cur_batt_power = self.curRef(1) - self.cur_delta - power;  % residual
                end
                % Saturate battery power
                batt_power = max(min(self.cur_batt_power, self.batt_params.power_max), self.batt_params.power_min);
                
                % Compute the next SoC with saturation
                last_soc = self.batt_soc;

%                 if tod >= 0 && tod <= 0.49                  % At 0 hour we will reset the Battery's SoC to middle level
%                     self.batt_soc = (self.batt_params.soc_min + self.batt_params.soc_max)/2;
%                 else
                    self.batt_soc = max(min(last_soc + self.batt_params.timestep * batt_power, ...
                        self.batt_params.soc_max), self.batt_params.soc_min);
%                end
                % The actual battery power, with both power and SoC constraints
                self.signals.setSignal('batt_power', self.curStep-1, (self.batt_soc - last_soc) / self.batt_params.timestep);
                self.signals.setSignal('batt_soc', self.curStep, self.batt_soc);
            end

            % Set default setpoints
            [ClgSP, HtgSP, KitchenClgSP, KitchenHtgSP, GuestClgSP, GuestHtgSP, SupplyAirSP, ChwSP] = self.defaultSetpoints(tod);
            
            % Save uncontrolled setpoints
            self.signals.setSignal('ClgSP', self.curStep, ClgSP);
            self.signals.setSignal('KitchenClgSP', self.curStep, KitchenClgSP);
            
            self.last_active = self.isActive();
            if self.last_active
                % The controller is active, and the reference is obtained
                % The optimization has the following optimization variables
                %   s = the sequence of DR signal values (length = curHorizon)
                %   Pb = the sequence of battery power values when we don't
                %           track the reference (when we track, Pb =
                %           residual).

                nPbVars = length(self.batt_power_nontracking);
                
                fprintf('At step %d: ', self.curStep);
                
                startnlptime = tic;
                % Create the variables
                optvar = self.create_optstructs(self.curHorizon, nPbVars, self.prevSol);
                
                % Simulate the system (naive)
                GSim = self.simulate_naive(optvar, tod, dayofweek);
                
                % Construct the battery constraints
                if self.batt_used
                    GBatt = self.battery_constraints(optvar);
                else
                    GBatt = {};
                end
                
                % Ramp limit constraints
                GRamp = self.rampConstraints(optvar);
                
                % Construct the objective function
                fobj = self.objfunc_tracking(optvar.Delta, optvar.P, optvar.Sigma, optvar.SOC);
                
                % Construct and solve the constrained NLP
                
                mysolver = 'ipopt';
                options = struct('quiet', true, ...
                                 'ipopt', struct('linear_solver', 'ma57', 'tol', 1e-6, 'acceptable_tol', 1e-4, 'max_iter', 400, 'print_level', 0));
                
                nlp = optisolve(fobj, [GSim, GBatt, GRamp], options, mysolver);
                
                %{
                mysolver = 'sqpmethod';
                options = struct('quiet', true, ...
                    'qpsol', 'qpoases', ...
                    'hessian_approximation', 'limited-memory', ...
                    'print_header', false, 'print_time', false, 'qpsol_options', struct('printLevel', 'none'));

                optvar = self.initialize_optvar(optvar, self.extract_solution(optvar));
                nlp = optisolve(fobj, [GSim, GBatt], options, mysolver);
                %}
                
                % Check the result
                status = nlp.solver.stats();
                if isstruct(status) && isfield(status, 'return_status')
                    fprintf('Constrained NLP: %g (%s)\n', toc(startnlptime), status.return_status);
                    %{
                    if ~strcmpi(status.return_status, 'Solve_Succeeded')
                        disp('Trying to recompute solution with fmincon...');
                        
                        % Reinitiliaze the variables with last solutions
                        % returned by Ipopt
                        optvar = self.initialize_optvar(optvar, self.extract_solution(optvar));
                        
                        mysolver = 'fmincon';
                        options = struct('quiet', true, ...
                            'fmincon', optimoptions('fmincon', 'TolX', 1e-6, 'Algorithm', 'sqp', 'MaxIter', 20));
                        
                        nlp_fmincon = optisolve(fobj, [GSim, GBatt], options, mysolver);
                        
                        if nlp_fmincon.fmincon_flag <= 0
                            % Not successful, try to solve again with Ipopt
                            disp('Trying to recompute solution with Ipopt...');
                            
                            % Reinitiliaze the variables with last solutions
                            % returned by fmincon
                            optvar = self.initialize_optvar(optvar, self.extract_solution(optvar));
                            
                            nlp.resolve();
                            
                            status = nlp.solver.stats();
                            if isstruct(status) && isfield(status, 'return_status')
                                fprintf('Final try with Ipopt: %s\n', status.return_status);
                            else
                                fprintf('Final try with Ipopt: unknown status\n');
                            end
                        end
                    end
                    %}
                else
                    fprintf('Constrained NLP: %g (unknown status)\n', toc(startnlptime));
                end
                
                self.prevSol = self.extract_solution(optvar);
                
                % Compute battery power for this step
                if self.batt_used
                    if nPbVars > 0 && self.batt_power_nontracking(1) == 1
                        % This step's battery power is fixed
                        % Note that the Pb variable in optimization is
                        % normalized into [-1,1], so we must de-normalize it.
                        result_Pb = self.prevSol.Pb;
                        self.cur_batt_power = (result_Pb(1)+1)*self.batt_params.norm_power_delta + self.batt_params.power_min;
                        self.cur_delta = NaN;
                    else
                        self.cur_batt_power = NaN;  % Just track the residual
                        solDelta = self.prevSol.Delta;
                        self.cur_delta = solDelta(1)*self.gpcommon.norm_y_delta;
                    end
                else
                    self.cur_batt_power = 0;
                    self.cur_delta = NaN;
                end
                
                % Delete the optivar objects in optvar
                self.delete_objects_in_struct(optvar);
                
                % Obtain the controlled setpoints and save and calculate the
                % denormalized values
                self.signals.setSignalUnnormalized('GuestClgSP', self.curStep, self.prevSol.GuestClgSP(1));
                GuestClgSP = self.signals.getSignalDenorm('GuestClgSP', self.curStep);
                
                self.signals.setSignalUnnormalized('SupplyAirSP', self.curStep, self.prevSol.SupplyAirSP(1));
                SupplyAirSP = self.signals.getSignalDenorm('SupplyAirSP', self.curStep);
                
                self.signals.setSignalUnnormalized('ChwSP', self.curStep, self.prevSol.ChwSP(1));
                ChwSP = self.signals.getSignalDenorm('ChwSP', self.curStep);

            else
                
                % Save controlled setpoints
                self.signals.setSignal('GuestClgSP', self.curStep, GuestClgSP);
                self.signals.setSignal('SupplyAirSP', self.curStep, SupplyAirSP);
                self.signals.setSignal('ChwSP', self.curStep, ChwSP);
                
                self.cur_batt_power = 0;  % Do not change the battery's state
                self.cur_delta = NaN;
            end
            
            % Send to building
            self.output('setpoints', [ClgSP, HtgSP, KitchenClgSP, KitchenHtgSP, GuestClgSP, GuestHtgSP, SupplyAirSP, ChwSP]);
            
            % Calculate and save the predicted mean and variance
            xtest = self.gpmodels.signalmodel.getIOVectors(self.curStep);
            if ~any(isnan(xtest))
                [y, vary] = self.gpmodels.gp.predict(xtest);
                
                % Denormalize and save them
                self.signals.setSignal('power_expected', self.curStep, postNorm(y, self.gpcommon.norm_y_min, self.gpcommon.norm_y_max));
                self.signals.setSignal('power_var', self.curStep, postNormVar(vary, self.gpcommon.norm_y_min, self.gpcommon.norm_y_max));
            end
                        
            self.signals.setSignal('tracking_delta', self.curStep, self.cur_delta);
            
            self.curStep = self.curStep + 1;
        end
                
        function finalizeSim(self)
            % Truncate vectors that store values
            self.signals.truncate(self.curStep-1);
            
            fprintf('Unconstrained NLP succeeded %d times; constrained NLP needed %d times.\n',...
                self.unconstrained_solved, self.constrained_solved);
        end
        
    end

    
    methods (Static)
        function optvar = initialize_optvar(optvar, X0)
            % Initialize the optivar objects in optvar by corresponding
            % fields in X0
            initialize_var('GuestClgSP');
            initialize_var('SupplyAirSP');
            initialize_var('ChwSP');
            initialize_var('Pb');
            initialize_var('Delta');
            initialize_var('P');
            initialize_var('Sigma');
            initialize_var('SOC');
            initialize_var('tS');
            
            function initialize_var(varname)
                if isfield(X0, varname) && ~isempty(X0.(varname)) && isfield(optvar, varname)
                    X0l = length(X0.(varname));
                    varlen = length(optvar.(varname));
                    if X0l > 1
                        initval = X0.(varname)(2:end);
                        if X0l - 1 >= varlen
                            optvar.(varname).setInit(initval(1:varlen));
                        else
                            % Repeat the last value
                            optvar.(varname).setInit([initval; repmat(initval(end),varlen-X0l+1,1)]);
                        end
                    else
                        % The previous value is a scalar -> repeat that
                        optvar.(varname).setInit(repmat(X0.(varname)(1),varlen,1));
                    end
                end
            end
        end
        
        function delete_objects_in_struct(S)
            % Delete all objects in a structure
            flds = fieldnames(S);
            for k = 1:numel(flds)
                if isa(S.(flds{k}), 'optivar')
                    delete(S.(flds{k}));
                end
            end
        end
        
        function X0 = extract_solution(optvar)
            % Extract the values of the solutions in the variables in
            % optvar
            X0 = structfun(@(v) v.value, optvar, 'UniformOutput', false);
        end
        
        function [ClgSP, HtgSP, KitchenClgSP, KitchenHtgSP, GuestClgSP, GuestHtgSP, SupplyAirSP, ChwSP] = defaultSetpoints(tod)
            tod = rem(tod, 24);
            % need this because some inputs will follow rule-based schedules
            if tod <= 7
                ClgSP = 30;
                KitchenClgSP = 30;
                GuestClgSP = 24;
                SupplyAirSP = 13;
                ChwSP = 6.7;
                HtgSP = 16;
                KitchenHtgSP = 16;
                GuestHtgSP = 21;
            else
                ClgSP = 24;
                KitchenClgSP = 26;
                GuestClgSP = 24;
                SupplyAirSP = 13;
                ChwSP = 6.7;
                HtgSP = 21;
                KitchenHtgSP = 19;
                GuestHtgSP = 21;
            end
        end
    end
    
    methods (Access=private)    
        % Create the optimization structures, which will be used during the
        % optimization
        %
        % This function tries to reuse variables if possible.
        %
        % horizon = current horizon
        % nPbVars = number of Pb variables when not tracking the ref
        % X0 = a structure of similar fields as optvar, which contains
        %       values for the variables from the previous iteration, which
        %       can be used for warmstarting the optimization.
        % nTracking = number of steps tracking the ref,
        %           should be = horizon - nPbVars
        function optvar = create_optstructs(self, horizon, nPbVars, X0)
            %optvar = self.optvar;   % Load the current variables
            optvar = struct();
            
            %if isfield(optvar, 'S') && length(optvar.S) >= horizon
            optvar.GuestClgSP = optivar(horizon, 1, 'GuestClgSP');
            setLb(optvar.GuestClgSP, -1);
            setUb(optvar.GuestClgSP, 1);

            optvar.SupplyAirSP = optivar(horizon, 1, 'SupplyAirSP');
            setLb(optvar.SupplyAirSP, -1);
            setUb(optvar.SupplyAirSP, 1);

            optvar.ChwSP = optivar(horizon, 1, 'ChwSP');
            setLb(optvar.ChwSP, -1);
            setUb(optvar.ChwSP, 1);
            
            if nPbVars > 0
                optvar.Pb = optivar(nPbVars, 1, 'Pb');
                setLb(optvar.Pb, -1);
                setUb(optvar.Pb, 1);
            end
            
            % The tracking slack (delta(i))
            optvar.Delta = optivar(length(self.batt_power_tracking), 1, 'Delta');
            % setInit(optvar.Delta, 0);
            
            optvar.P = optivar(horizon, 1, 'P');
            
            optvar.Sigma = optivar(horizon, 1, 'Sigma');
            setLb(optvar.Sigma, 0);
            
            optvar.SOC = optivar(horizon, 1, 'SOC');
            setUb(optvar.SOC, self.batt_params.soc_max_norm);
            setLb(optvar.SOC, self.batt_params.soc_min_norm);
            
            % optvar.tP = optivar(nTracking, 1, 'tP');
            
            optvar.tS = optivar(horizon, 1, 'tS');
            setLb(optvar.tS, 0);
            
            % Set the initial values of the variables if possible
            if exist('X0', 'var') && isstruct(X0)
                optvar = self.initialize_optvar(optvar, X0);
            end
        end
        
        function inputmatrix = construct_data(self, todseries, dowseries, spseries, k, optvar)
            % Construct the input matrix for simulation at step k (starting
            % from 0 = current time step).
            
            % Create the input substitutions
            inputsubs = struct;
            
            inputsubs.TotalLoad.start = self.curStep;
            inputsubs.TotalLoad.values = optvar.P;
            
            inputsubs.GuestClgSP.start = self.curStep;
            inputsubs.GuestClgSP.values = optvar.GuestClgSP;
            
            inputsubs.SupplyAirSP.start = self.curStep;
            inputsubs.SupplyAirSP.values = optvar.SupplyAirSP;
            
            inputsubs.ChwSP.start = self.curStep;
            inputsubs.ChwSP.values = optvar.ChwSP;
            
            ts = self.curStep:(self.curStep+self.curHorizon-1);
            
            self.signals.setSignal('TOD', ts, todseries);
            self.signals.setSignal('DOW', ts, dowseries);
            self.signals.setSignal('ClgSP', ts, spseries(:,1));
            self.signals.setSignal('KitchenClgSP', ts, spseries(:,2));
            
            assert(~any(isnan(todseries)) && ~any(isnan(dowseries)) && ~any(isnan(spseries(:,1))) && ~any(isnan(spseries(:,2))));
            
            % Get the regressor vector from the SignalsModel object            
            inputmatrix = self.gpmodels.signalmodel.getIOVectors(self.curStep+k, 'inputsubs', inputsubs, 'type', 'casadi.MX');
            
            for kk = 1:numel(inputmatrix)
                if isa(inputmatrix(kk), 'casadi.MX')
                    if is_constant(inputmatrix(kk))
                        assert(~isnan(full(to_DM(inputmatrix(kk)))), '%d is NaN', kk);
                    end
                else
                    assert(~isnan(inputmatrix(kk)), '%d is NaN', kk);
                end
            end
        end
        
        function G = simulate_naive(self, optvar, curHour, curDay)
            % Simulate the GP model with naive algorithm, using symbolic
            % variables.  Returns the constraints between variables that
            % represent the dynamics of the GP model over the horizon.

            need_var = isfield(optvar, 'Sigma');
            if need_var
                step_size = 2;
            else
                step_size = 1;
            end
            
            % Construct the future values for several inputs
            todseries = (0:(self.curHorizon-1))'*self.timestep + curHour;
            dowseries = zeros(size(todseries)) + curDay;
            over24 = todseries > 24;
            dowseries(over24) = dowseries(over24) + fix(todseries(over24)/24);
            dowseries = rem(dowseries, 7);
            todseries(over24) = rem(todseries(over24), 24);
            
            spseries = zeros(self.curHorizon, 5);
            for kk = 1:self.curHorizon
                [spseries(kk,1), ~, spseries(kk,2), ~, spseries(kk,3), ~, spseries(kk,4), spseries(kk,5)] = ...
                    self.defaultSetpoints(todseries(kk));
            end
            
            G = cell(1, step_size * self.curHorizon);
            
            % Perform naive simulation iteratively
            for k = 1:self.curHorizon
                % Construct the input for the GP for this step
                inputmatrix = self.construct_data(todseries, dowseries, spseries, k-1, optvar);
                
                % Compute the mean and/or variance and assign the
                % constraints
                [ymu, ys2] = self.gpmodels.gp.predict(inputmatrix);
                
                if isa(ymu, 'casadi.MX')
                    if is_constant(ymu)
                        assert(~isnan(full(to_DM(ymu))));
                    end
                else
                    assert(~isnan(ymu));
                end
                
                G{(k-1)*step_size+1} = (optvar.P(k) == ymu);
                if need_var
                    G{k*2} = (optvar.Sigma(k)^2 == ys2);
                end
                
                %{
                G{(k-1)*step_size+1} = (self.gpmodels(k).fymu(inputmatrix) == optvar.P(k));
                if need_var
                    G{k*2} = (self.gpmodels(k).fys2(inputmatrix) == optvar.Sigma(k)^2);
                end
                %}
            end
        end
        
        function cost = objfunc_tracking(self, Delta, P, Sigma, SOC) %#ok<INUSD>
            % The objective function.
            % drsignal is the DR signal values up to the horizon.
            %       sum_t (wdelta*delta(t)^2 + ...
            %              wPb*(r(t)-delta(t)-E[P(t)])^2 ) ...
            %       - wvar*p
            % where delta(t) is the tracking slack at time step t.
            % w* are the weights.
            
            % Now compute the cost
            include_var = abs(self.wvar) > 1e-6;
            include_delta = abs(self.wdelta) > 1e-6;
            include_Pb = abs(self.wPb) > 1e-6;
            cost = {};
            
            if include_delta && self.batt_used
                if isvector(Delta) && numel(Delta) > 1
                    cost = [cost, {sqrt(2*self.wdelta)*Delta}];
                else
                    cost = [cost, {self.wdelta * sum_square(Delta)}];
                end
            end
            
            if include_Pb
                norm_ref = (self.curRef - self.gpcommon.norm_y_min)/self.gpcommon.norm_y_delta - 1;    % normalized reference
                if self.batt_used
                    v = norm_ref - Delta - P(self.curOutputIdx);
                else
                    v = norm_ref - P(self.curOutputIdx);
                end
                if isvector(v) && numel(v) > 1
                    cost = [cost, {sqrt(2*self.wPb)*v}];
                else
                    cost = [cost, {self.wPb * sum_square(v)}];
                end
            end
            
            if include_var
                v = Sigma(self.curOutputIdx);
                if isvector(v) && numel(v) > 1
                    cost = [cost, {sqrt(2*self.wvar)*v}];
                else
                    cost = [cost, {self.wvar * sum_square(v)}];
                end
            end
        end
        
        %{
        function cost = objfunc_reduction(self, ymus, yvars)
            % The objective function, by simulating the system then
            % calculate the cost.
            % It uses naive multi-step simulation of given GP model.
            % drsignal is the DR signal values up to the horizon.
            %       sum_t alpha*(r(t)-E[P(t)])^2 + ...
            %           (1-alpha)*(E[P(t)]-r(t)+beta*sqrt(var[P(t)]))
            error('Currently not implemented.');
            
            % Now compute the cost
            cost = sum(self.alpha*(self.curRef - self.cur_opt_data.py).^2 + ...
                (1-self.alpha)*(self.cur_opt_data.py - self.curRef + self.beta*sqrt(self.cur_opt_data.pvar)));
        end
        %}
        
        function G = battery_constraints(self, optvar)
            % Returns the nonlinear constraints for the battery
                        
            norm_curRef = (self.curRef - self.gpcommon.norm_y_min)/self.gpcommon.norm_y_delta - 1;    % normalized reference
            Pb_tracking = norm_curRef - optvar.Delta - optvar.P(self.batt_power_tracking);  % Pb during tracking, which takes the tracking error
            % Pb_tracking is the r(k) - mu(k) in the math formulation
            
            % Expected value of Pb for the horizon
            Pb_vals = casadi.MX(self.curHorizon, 1);
            if ~isempty(self.batt_power_nontracking)
                Pb_vals(self.batt_power_nontracking) = ...
                    ((optvar.Pb+1)*self.batt_params.norm_power_delta + self.batt_params.power_min)/self.gpcommon.norm_y_delta;
            end
            Pb_vals(self.batt_power_tracking) = Pb_tracking;
            
            % bounds of battery power, with variance for steps in self.batt_power_tracking
            % Note that Sigma is the std variance, not sigma^2
            Pb_var_tracking = optvar.Sigma(self.batt_power_tracking);
            
            % The chance constraints for Pb, only when we track
            G = {...
                Pb_tracking - self.optpar.PhiEps*Pb_var_tracking <= self.batt_params.power_max/self.gpcommon.norm_y_delta, ...
                Pb_tracking - self.optpar.PhiEpsCompl*Pb_var_tracking >= self.batt_params.power_min/self.gpcommon.norm_y_delta};
            
            % The last inequality sets an upper bound on certain Sigma
            % variables
            SigmaUB = inf(size(optvar.Sigma));
            SigmaUB(self.batt_power_tracking) = (self.batt_params.power_min-self.batt_params.power_max)/self.gpcommon.norm_y_delta/(2*self.optpar.PhiHalf);
            setUb(optvar.Sigma, SigmaUB);
            
            TrackingIdx = zeros(self.curHorizon, 1);
            TrackingIdx(self.batt_power_tracking) = 1;
            
            % Calculate the SoC over the horizon
            soc_mean = self.batt_soc/self.batt_params.timestep/self.gpcommon.norm_y_delta;
            for k = 1:self.curHorizon
                G = [G, {...
                    soc_mean + Pb_vals(k) == optvar.SOC(k), ...  % The Battery dynamic model
                    optvar.tS(k)^2 >= dot(TrackingIdx(1:k), optvar.Sigma(1:k).^2)}];  %#ok<AGROW> % Lower bound constraint on tS => SOC but squared
                soc_mean = optvar.SOC(k);  % The current SOC mean
            end
            
            G = [G, {...
                optvar.SOC - self.optpar.PhiEpsSOC*optvar.tS <= self.batt_params.soc_max_norm, ...
                optvar.SOC - self.optpar.PhiEpsComplSOC*optvar.tS >= self.batt_params.soc_min_norm}];
            
            % The last inequality sets an upper bound on tS
            setUb(optvar.tS, (self.batt_params.soc_min_norm-self.batt_params.soc_max_norm)/(2*self.optpar.PhiHalfSOC));
            
        end
        
        function G = rampConstraints(self, optvar)
            % Add the ramp limit constraints
            if ~isfinite(self.ramplimit) || self.curHorizon < 1
                G = {};
                return;
            end
            
            % Ramp limit w.r.t the last setpoint
            lastSP = self.signals.getSignal('ChwSP', self.curStep-1);
            G = {optvar.ChwSP(1) - lastSP <= self.ramplimit, optvar.ChwSP(1) - lastSP >= -self.ramplimit};
            
            if self.curHorizon < 2, return; end
            
            % Construct the matrix that calculates the difference between
            % consecutive elements of the setpoint vector
            k = self.curHorizon-1;
            A = full(spdiags([-ones(k,1), ones(k,1)], [0 1], k, length(optvar.ChwSP)));
            
            % Construct the inequalities
            G = [G, {A*optvar.ChwSP <= self.ramplimit, A*optvar.ChwSP >= -self.ramplimit}];
        end
        
        function b = isActive(self)
            % Check if the controller is activated, and save the
            % reference values
            % Also prepares index vectors for power outputs and battery
            % power values.
            % self.batt_power_tracking = indices within the horizon that we
            %   track the reference -> battery should only take the residual
            % self.batt_power_nontracking = indices within the horizon that
            %   we don't track the reference -> battery power can be
            %   optimized
            % The union of these should be the entire horizon.
            
            % Beyond simulation length, or before max lag kicks in, disable
            % the controller.
            if self.curStep > self.reflen % || self.curStep <= self.max_lag
                b = false;
                return;
            end
            
            end_step = min(self.curStep + self.active_ahead - 1, self.reflen);
            b = ~all(isnan(self.ref(self.curStep:end_step)));
            if b
                end_step = min(self.curStep + self.horizon - 1, self.reflen);
                nextref = self.ref(self.curStep:end_step);
                self.curOutputIdx = find(~isnan(nextref));
                
                % Calculate the actual horizon, which can be truncated from
                % 'horizon' towards the end of the DR duration.
                self.curHorizon = max(self.curOutputIdx);
                
                self.curRef = nextref(self.curOutputIdx);   % The reference values when tracking
                
                self.batt_power_tracking = self.curOutputIdx;   % When tracking
                self.batt_power_nontracking = find(isnan(nextref(1:self.curHorizon))); % When not tracking, only before curHorizon
            end
        end
        

    end
    
end