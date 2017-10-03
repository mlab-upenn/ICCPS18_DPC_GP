%% Create an mlepProcess instance and configure it

file = 'simBuilding';

ep = mlepProcess;
ep.arguments = {file, 'USA_IL_Chicago-OHare.Intl.AP.725300_TMY3'};
ep.acceptTimeout = 6000; % in milliseconds

VERNUMBER = 2;  % version number of communication protocol (2 as of
                % E+ 8.1.0)

%% Start EnergyPlus cosimulation
[status, msg] = ep.start;

if status ~= 0
    error('Could not start EnergyPlus: %s.', msg);
end

%% The main simulation loop

EPTimeStep = 4;
SimDays = 30+31+31;
deltaT = (60/EPTimeStep)*60;
kStep = 1;  % current simulation step
MAXSTEPS = SimDays*24*EPTimeStep;  % max simulation time = 7 days

% variables for plotting:
outputs = nan(9,MAXSTEPS);
inputs = nan(4,MAXSTEPS);

% Initialize: parse it to obtain building outputs
packet = ep.read;
if isempty(packet)
    error('could not read outputs from E+.');
end
[flag, eptime, outinit] = mlepDecodePacket(packet);
if flag ~= 0, error('check output flag'); end

rng(0);

nConst = 4;
DRS = [];
for ids = 1:SimDays
    DRS = [DRS;lhsdesign(1*96,4)];
end

% data_type = 'unconstrained';
% data_type = 'constrained';
data_type = 'nominal';
% data_type = 'ramped';
rampVal = 2;

if strcmp(data_type, 'unconstrained') || strcmp(data_type, 'ramped')
    SPmin = repmat([23, 0, 11, 3.7], size(DRS,1),1);
    SPmax = repmat([28, 1, 15, 9.7], size(DRS,1),1);
    DRS = SPmin + (SPmax-SPmin).*DRS;
elseif strcmp(data_type, 'constrained')
    SPmin = repmat([23, 0, 12, 5.2], size(DRS,1),1);
    SPmax = repmat([28, 1, 13, 8.2], size(DRS,1),1);
    DRS = SPmin + (SPmax-SPmin).*DRS;
end



while kStep <= MAXSTEPS    

    % BEGIN Compute next set-points
    dayTime = mod(eptime, 86400);  % time in current day    
    
    % Baseline Schedule.
    % 1-sun, 2-mon, ..., 7-sat
    if kStep ==1    % change if this changed in idf
        dayWeek = 'AllOtherDays';
    else
        if outputs(5,kStep-1)>0 || outputs(4,kStep-1)==1
            dayWeek = 'AllOtherDays';
        elseif outputs(4,kStep-1)==7
            dayWeek = 'Saturday';
        else
            dayWeek = 'Weekdays';
        end
    end
    clgstp = schedule('ClgSP', dayTime, dayWeek);
    lgtstp = schedule('LgtSP', dayTime, dayWeek);
    sat = 13;
    cwstp = 6.7;
    
    SP = [clgstp, lgtstp, sat, cwstp];
        
    
    if ~strcmp(data_type, 'nominal')
        SP(1) = DRS(kStep,1);
        SP(3) = DRS(kStep,3);
        SP(4) = DRS(kStep,4);
    end
    
    if strcmp(data_type, 'ramped') && kStep>1
        ChwMin = 3.7;
        ChwMax = 9.7;
        ChwPrev = inputs(4,kStep-1);
        SP(4) = max(ChwMin,ChwPrev-rampVal)+(min(ChwMax,ChwPrev+rampVal)-max(ChwMin,ChwPrev-rampVal))*rand;
    end
    
    inputs(:,kStep) = SP;
    
    % Write to inputs of E+
    ep.write(mlepEncodeRealData(VERNUMBER, 0, (kStep-1)*deltaT, SP));
    
    % Read a data packet from E+
    packet = ep.read;
    if isempty(packet)
        error('Could not read outputs from E+.');
    end
    
    % Initialize: parse it to obtain building outputs
    [flag, eptime, outputs(:,kStep)] = mlepDecodePacket(packet);
    if flag ~= 0, break; end

    kStep = kStep + 1;
    
end

% Stop EnergyPlus
ep.stop;

disp(['Stopped with flag ' num2str(flag)]);
