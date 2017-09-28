classdef SignalsModel < handle
    % Class to manage the MISO signals for (machine learning) models.
    %
    %   (C) 2017 by Truong X. Nghiem (truong.nghiem@gmail.com)
    
    properties(GetAccess=public,SetAccess=private)
        %m_maxlag;           % Max lag of all auto-regressive
        
        m_inputs;           % Cell array defining the inputs and their lags
        m_output;           % Name of the output signal (string)
        m_totalinputs;      % Total number of inputs
        
        m_Signals;          % the values of the signals, as fields of this structure
    end
        
    methods
        function self = SignalsModel(other)
            % Construct the SignalsModel object.
            % If other is given, which is a SignalsModel object or a
            % SignalsValues object, the Signal value source is shared.
            % Else, create a new SignalsValues object to store the values.

            %self.m_maxlag = 0;
            self.m_inputs = {};
            self.m_output = '';
            self.m_totalinputs = 0;
            
            if nargin > 0
                if isa(other, 'SignalsValues')
                    self.m_Signals = other;
                elseif isa(other, 'SignalsModel')
                    self.m_Signals = other.Signals;
                else
                    error('Invalid argument of class %s.', class(other));
                end
            else
                self.m_Signals = SignalsValues();
            end
        end
        
        function setInputs(self, inputs)
            %	setInputs(self, (cell) inputs)
            %
            % Set the input signals for modelling.
            % Inputs:
            %   inputs  - a cell array that specifies the inputs /
            %               features. Each cell is either:
            %               + a cell {name, lags} where name is the field
            %               of the input in data structure (must be a valid
            %               field name), and lags is a vector of lag values
            %               as in lagmatrix.m. See the notes below.
            %               + a single string which is equivalent to
            %               {name,0}.
            %               The order of the inputs will be kept as is.
            %
            % Example:
            %	obj.setInputs({'t', {'u', [0,2,4]}});
            %	In this case the input vector (of current time-step k):
            %       x(k) = [t(k), u(k), u(k-2) u(k-4)]

            assert(iscell(inputs), 'INPUTS must be a cell array.');
            
            % Make inputs have a standard form {{name, lags}, ...}; also
            % check the names.
            self.m_totalinputs = 0;
            %self.m_maxlag = 0;
            for k = 1:numel(inputs)
                if ischar(inputs{k})
                    if ~self.m_Signals.isSignal(inputs{k})
                        warn('SignalsModel:setInputs', 'The input "%s" is not an existing signal. It will be added.', inputs{k});
                        self.m_Signals.addSignal(inputs{k});
                    end
                    inputs{k} = {inputs{k}, 0};
                    self.m_totalinputs = self.m_totalinputs + 1;
                else
                    assert(iscell(inputs{k}) && numel(inputs{k}) == 2 && ...
                        ischar(inputs{k}{1}) && isnumeric(inputs{k}{2}) && ...
                        isvector(inputs{k}{2}), ...
                        'Invalid element %d of INPUTS.', k);
                    
                    if ~self.m_Signals.isSignal(inputs{k}{1})
                        warn('SignalsModel:setInputs', 'The input "%s" is not an existing signal. It will be added.', inputs{k}{1});
                        self.m_Signals.addSignal(inputs{k}{1});
                    end

                    assert(~isempty(inputs{k}{2}) && all(inputs{k}{2} >= 0), 'Invalid lags for input %d.', k);
                    self.m_totalinputs = self.m_totalinputs + numel(inputs{k}{2});
                    
                    % self.m_maxlag = max(self.m_maxlag, max(inputs{k}{2}));
                end
            end
            
            self.m_inputs = inputs;
        end
        
        function setOutput(self, output)
            %	setOutput(self,(str) output)
            %
            % set the output (target) signal for modelling
            
            assert(ischar(output), 'OUTPUT must be the name of the output signal.');
            
            if ~self.m_Signals.isSignal(output)
                warn('SignalsModel:setOutput', 'The output "%s" is not an existing signal. It will be added.', output);
                self.m_Signals.addSignal(output);
            end
            self.m_output = output;
        end
        
        function [X, Y] = getIOVectors(self, timesteps, varargin)
            % Returns the input and output vectors at given time steps.
            %       [X, Y] = getIOVectors(self, TS, options)
            % Inputs:
            %   TS      - the time-step vector; if ignored or empty then
            %             all time steps.
            %   options - name-value pairs of options. See below.
            %
            % Outputs:
            %   X   - the input vectors as rows of this matrix.
            %   Y   - the output values in this column vector.
            %
            % Options: the following options are supported:
            %   'stepsahead' : the number of steps into the future
            %           of the output. Both input and output vectors are
            %           shifted by this number of steps, that is at time
            %           step k and with N steps ahead, the output is at
            %           time step k + N, and the input values are counted
            %           from time step k + N (so a lag of 0 means time step
            %           (k+N), a lag of 1 means time step (k+N-1)). Default
            %           value of stepsahead is 0. It must be non-negative.
            %           For example if current time is t = 2.0 and with
            %           sampling time of 15 minutes, stepsahead = 2, then
            %           all lag values are counted from 2h30', i.e., lag
            %           value of 0 is at 2h30', lag value of 1 is at 2h15',
            %           and so on.
            %   'inputsubs' : only if k is a scalar, substitutes for
            %           certain input values can be provided in this
            %           structure. Each field of a signal name is a
            %           structure of two fields 'start' and 'values', which
            %           specify substitute values for the corresponding
            %           signal starting at time step 'start'. This is
            %           useful for predictive optimization and control.
            %   'type' : the values from the SignalsValues source will be
            %           transformed into this given type, if provided, by
            %           calling feval('type', values). This can be used to
            %           convert numerical values into symbolic constants
            %           used for symbolic calculations and optimization.
            %   'removeNaNs' : rows which contain NaN values will be
            %           removed from X and Y.
            
            params = inputParser;
            addParameter(params, 'stepsahead', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
            % addParameter(params, 'excepts', {}, @(x) ischar(x) || iscellstr(x));
            addParameter(params, 'inputsubs', struct(), @isstruct);
            addParameter(params, 'type', '', @(x) ischar(x) || isa(x, 'function_handle'));
            addParameter(params, 'removeNaNs', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
            params.parse(varargin{:});
            
            if nargin < 2
                timesteps = [];
            end
            if ~isempty(timesteps)
                assert(isvector(timesteps) && all(timesteps > 0), 'SignalsModel:getIOVectors', 'Invalid TIMESTEPS.');
            end
            
            with_subs = ~isempty(fieldnames(params.Results.inputsubs));
            if with_subs && ~isscalar(timesteps)
                warning('SignalsModel:getIOVectors', 'When TIMESTEPS is not scalar, input substitution is not supported and will be ignored.');
                with_subs = false;
            end
            
            if isempty(timesteps)
                datalen = getSize(self.m_Signals);
            else
                datalen = max(timesteps) + params.Results.stepsahead;
            end
            
            inputidx = 1:datalen;
            
            if with_subs
                assert(all(structfun(@(s) isfield(s, 'values') && isfield(s, 'start') && isscalar(s.start) && s.start > 0, params.Results.inputsubs)),...
                    'SignalsModel:getIOVectors', 'Invalid input substitution structure.');
            end
            
            % Loop through the inputs and construct the input matrix X
            X = nan(datalen, self.m_totalinputs);    % pre-allocate X matrix
            if ~isempty(params.Results.type)
                % Convert type of X
                X = feval(params.Results.type, X);
            end
            curIdx = 1;
            for k = 1:numel(self.m_inputs)
                thelags = self.m_inputs{k}{2} - params.Results.stepsahead;
                numlags = numel(thelags);
                thefield = self.m_inputs{k}{1};
                thesignal = self.m_Signals.getSignal(thefield, inputidx);
                if with_subs && isfield(params.Results.inputsubs, thefield)
                    Xidx = lagmatrix(inputidx, thelags);
                    
                    % start and end indices of the substitute values
                    startidx = params.Results.inputsubs.(thefield).start;
                    endidx = startidx + length(params.Results.inputsubs.(thefield).values) - 1;
                    
                    for kk = 1:numlags
                        % Assign values from Signals source object for those
                        % indices lower than the start index of the substitute
                        % values.
                        Xidx_from = find((Xidx(:,kk) < startidx) | (Xidx(:,kk) > endidx));
                        if ~isempty(Xidx_from)
                            X(Xidx_from, kk+curIdx-1) = thesignal(Xidx(Xidx_from,kk));
                        end
                        
                        % Assign values from substitute values.
                        % don't use ~Xidx_from as it will include NaNs
                        Xidx_from = find((Xidx(:,kk) >= startidx) & (Xidx(:,kk) <= endidx));
                        if ~isempty(Xidx_from)
                            X(Xidx_from, kk+curIdx-1) = ...
                                params.Results.inputsubs.(thefield).values(Xidx(Xidx_from,kk) - startidx + 1);
                        end
                    end
                else
                    X(:, curIdx:curIdx+numlags-1) = lagmatrix(thesignal, thelags);
                end
                curIdx = curIdx + numlags;
            end
            assert(curIdx == self.m_totalinputs + 1);
            
            if ~isempty(timesteps)
                X = X(timesteps, :);
            end
            
            % Construct the output Y
            if nargout > 1
                Y = nan(datalen, 1);    % pre-allocate Y matrix
                if ~isempty(params.Results.type)
                    % Convert type of Y
                    Y = feval(params.Results.type, Y);
                end
                
                thesignal = self.m_Signals.getSignal(self.m_output, inputidx);
                
                if with_subs && isfield(params.Results.inputsubs, self.m_output)
                    Yidx = lagmatrix(inputidx, -params.Results.stepsahead);
                    
                    startidx = params.Results.inputsubs.(self.m_output).start;
                    endidx = startidx + length(params.Results.inputsubs.(self.m_output).values) - 1;
                    
                    % Assign values from Signals source object for those
                    % indices lower than the start index of the substitute
                    % values.
                    Yidx_from = find((Yidx(:) < startidx) | (Yidx(:) > endidx));
                    Y(Yidx_from) = thesignal(Yidx(Yidx_from));
                    
                    % Assign values from substitute values.
                    % don't use ~Yidx_from as it will include NaNs
                    Yidx_from = find((Yidx(:) >= startidx) & (Yidx(:) <= endidx));
                    Y(Yidx_from) = ...
                        params.Results.inputsubs.(self.m_output).values(Yidx(Yidx_from) - startidx + 1);
                else
                    Y = lagmatrix(thesignal, -params.Results.stepsahead);
                end
                
                if ~isempty(timesteps)
                    Y = Y(timesteps, :);
                end
            end
            
            % Remove NaN samples
            if params.Results.removeNaNs
                nanRows = any(isnan(X), 2);    % Xnan = column vector indicating NaN rows in X
                if nargout > 1
                    nanRows = nanRows | isnan(Y);  % NaN rows in X or Y
                    Y(nanRows) = [];
                end
                
                % Remove these rows
                X(nanRows, :) = [];
            end
        end
        
    end
end
