classdef SignalsValues < handle
    %SIGNALSVALUES Class to store signal values for SignalsModel.
    %   Each signal has a unique name, which is a field / property of an
    %   instance of this class. Multiple SignalsModel objects can share the
    %   same SignalsValues storage.
    %
    %   (C) 2017 by Truong X. Nghiem (truong.nghiem@gmail.com)
    
    properties (Access=private)
        m_curSize = 0; % current size of signal vectors
        m_growstep = 1024; % how to grow vectors
        m_signals = struct;
                
        % Optional parameters for normalizing signals.
        % It's a structure of the min and max values used for normalizing
        % the given fields.
        % It's assumed that the values stored in this object are
        % normalized if the corresponding normalization parameters are
        % given..
        m_normalize;
    end
    
    methods
        function self = SignalsValue(growstep)
            self.m_normalize = struct;
            if nargin > 0
                assert(growstep > 0);
                self.m_growstep = growstep;
            end
        end
        
        function grow(self, newsize)
            % Grow all signal vectors to at least the new size.
            if self.m_curSize >= newsize, return; end
            newsize = ceil(newsize/self.m_growstep)*self.m_growstep;
            self.m_curSize = newsize;
            allprops = fieldnames(self.m_signals);
            for k = 1:numel(allprops)
                self.m_signals.(allprops{k})(end+1:newsize) = NaN;
            end
        end
        
        function addSignal(self, name, values, normalization)
            % Add a new signal with given name and optionally initial
            % values.
            % Normalization parameters [min, max] can also be provided,
            % which will be saved and used to normalize the values of this
            % signal going forward. The initial values, if provided, will
            % be normalized before stored.
            assert(~isfield(self.m_signals, name), 'Signal %s already exists.', name);
            
            if nargin > 3
                if ~isempty(normalization)
                    normalization = normalization(:);
                    assert(length(normalization) == 2 && normalization(1) <= normalization(2),...
                        'SignalsValues:addSignal', 'Invalid normalization parameters.');
                    self.m_normalize.(name) = struct('min', normalization(1), 'max', normalization(2));
                end
            else
                normalization = [];
            end
            
            if nargin > 2 && ~isempty(values)
                newsize = numel(values);
                self.grow(newsize);
                self.m_signals.(name) = NaN(self.m_curSize, 1);
                if isempty(normalization)
                    self.m_signals.(name)(1:newsize) = vec(values);
                else 
                    self.m_signals.(name)(1:newsize) = preNorm(vec(values), normalization(1), normalization(2));
                end
            else
                self.m_signals.(name) = NaN(self.m_curSize, 1);
            end
        end
        
        function addMultiSignals(self, S, normparams)
            % Add signals from a structure. Each field of the given
            % structure will become a signal. If a signal already exists,
            % an error will occur.
            % Optional normalization parameters can be provided, which is a
            % structure whose fields are structure of 'min' and 'max' for
            % each signal.
            assert(isstruct(S), 'A structure of signals is expected.');
            if nargin > 2 && ~isempty(normparams)
                assert(isstruct(normparams) && ...
                    all(structfun(@(s) isfield(s, 'min') && isfield(s, 'max'), normparams)), ...
                    'SignalsValues:addMultiSignals', 'Invalid normalization parameter structure.');
            else
                normparams = struct;
            end
            
            flds = fieldnames(S);
            for k = 1:numel(flds)
                f = flds{k};
                if isfield(normparams, f)
                    self.addSignal(f, S.(f), [normparams.(f).min, normparams.(f).max]);
                else
                    self.addSignal(f, S.(f));
                end
            end
        end
        
        function b = isSignal(self, names)
            % Check if a given signal, or set of signals, exists.
            b = isfield(self.m_signals, names);
        end
        
        function s = getSignalNames(self)
            % Returns names of all signals as a cell array of strings.
            s = fieldnames(self.m_signals);
        end
        
        function v = getSignal(self, name, idx)
            % Get values of given signal with given index idx (if omitted
            % then the entire signal vector).
            % The vectors will grow if necessary so that idx is always
            % valid (uninitialized values will be NaN); so beware of
            % overflowing the memory.
            if nargin < 3
                v = self.m_signals.(name);
            else
                % grow if necessary
                if islogical(idx)
                    self.grow(numel(idx));
                else
                    self.grow(max(idx(:)));
                end
                v = self.m_signals.(name)(idx);
            end
        end
        
        function v = getSignalDenorm(self, name, varargin)
            % Similar to getSignal() but will denormalize the values if
            % corresponding normalization parameters are provided.
            v = self.getSignal(name, varargin{:});
            if isfield(self.m_normalize, name)
                % Denormalize the values
                v = postNorm(v, self.m_normalize.(name).min, self.m_normalize.(name).max);
            end
        end
        
        function setSignalUnnormalized(self, name, idx, v)
            % Set values of given signal with given index idx (if empty
            % then the entire signal vector).
            % This function ignores the normalization parameters, if exist
            % for the given signal name. IOW, the values in V should be
            % already normalized.
            if isempty(idx)
                self.m_signals.(name) = v;
            else
                % grow if necessary
                if islogical(idx)
                    self.grow(numel(idx));
                else
                    self.grow(max(idx(:)));
                end
                self.m_signals.(name)(idx) = v;
            end
        end
        
        function setSignal(self, name, idx, v)
            % Set values of given signal with given index idx (if empty
            % then the entire signal vector).
            % If normalization parameters are given for this signal, the
            % provided values are assumed to be unnormalized and will be
            % normalized before being stored in this object.
            if isfield(self.m_normalize, name)
                self.setSignalUnnormalized(name, idx, ...
                    preNorm(v, self.m_normalize.(name).min, self.m_normalize.(name).max));
            else
                self.setSignalUnnormalized(name, idx, v);
            end
        end
        
        function s = getAllSignals(self, idx)
            % Returns all signals in a structure.
            % If idx is given, only values in those indices will be
            % returned.
            if nargin == 1
                s = self.m_signals;
                return;
            end
            
            s = structfun(@(v) v(idx), self.m_signals, 'UniformOutput', false);
        end
        
        function s = getAllSignalsDenorm(self, idx)
            % Similar to getAllSignals() but will denormalize the values if
            % corresponding normalization parameters are provided.
            s = self.getAllSignals(idx);
            flds = fieldnames(s);
            hasnorm = find(isfield(self.m_normalize, flds));
            for kk = hasnorm
                f = flds{kk};
                s.(f) = postNorm(s.(f), self.m_normalize.(f).min, self.m_normalize.(f).max);
            end
        end
        
        function truncate(self, N)
            % Truncate all signals to a given length (N <= current length).
            if N >= self.m_curSize, return; end
            
            % Truncate
            allprops = fieldnames(self.m_signals);
            for k = 1:numel(allprops)
                self.m_signals.(allprops{k})(N+1:end) = [];
            end
            
            self.m_curSize = N;
        end
        
        function s = getSize(self)
            % Returns current size of the signals.
            s = self.m_curSize;
        end
    end
    
end

