
clear; close all;
rng(1);

[YY, MM, DD, HH, MINS, ~] = datevec(now);

%% Create an mlepProcess instance and configure it

building = 'LargeOffice';
ctrl_vars = {'ClgSP', 'SupplyAirSP', 'ChwSP'};

cd(['ep_model/' building '/'])
eplusfile = 'Building';

ep = mlepProcess;
ep.arguments = {eplusfile, 'USA_IL_Chicago-OHare.Intl.AP.725300_TMY3'};
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
SimDays = 30+31+31+30+31;   % Jun-Oct
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

% data_type = 'constrained';
% data_type = 'nominal';
% data_type = 'unconstrained';
data_type = 'ramped';
rampVal = 2;

if strcmp(data_type, 'unconstrained') || strcmp(data_type, 'ramped')
    ClgMin = 23;
    ClgMax = 28;
    LgtMin = 0;
    LgtMax = 1;
    SupplyAirMin = 11;
    SupplyAirMax = 15;
    ChwMin = 3.7;
    ChwMax = 9.7;

elseif strcmp(data_type, 'constrained')
    ClgMin = 23;
    ClgMax = 28;
    LgtMin = 0;
    LgtMax = 1;
    SupplyAirMin = 12;
    SupplyAirMax = 13;
    ChwMin = 5.2;
    ChwMax = 8.2;

end

sample_type = 'uniform';
% sample_type = 'prbs';
switch sample_type
    case 'uniform'
        ClgSPrand = ClgMin+(ClgMax-ClgMin)*rand(1,MAXSTEPS);
        LtgSPrand = LgtMin+(LgtMax-LgtMin)*rand(1,MAXSTEPS);
        SupplyAirSPrand = SupplyAirMin+(SupplyAirMax-SupplyAirMin)*rand(1,MAXSTEPS);
        ChwSPrand = ChwMin+(ChwMax-ChwMin)*rand(1,MAXSTEPS);
    case 'prbs'
        ClgSPrand = postNorm(idinput(MAXSTEPS,'prbs')', ClgMin, ClgMax);
        LtgSPrand = postNorm(idinput(MAXSTEPS,'prbs')', LgtMin,LgtMax);
        SupplyAirSPrand = postNorm(idinput(MAXSTEPS,'prbs')', SupplyAirMin, SupplyAirMax);
        ChwSPrand = postNorm(idinput(MAXSTEPS,'prbs')', ChwMin, ChwMax);
        ChwSPprbs = idinput(MAXSTEPS,'prbs')';
end

tic;

while kStep <= MAXSTEPS    

    % BEGIN Compute next set-points
    dayTime = mod(eptime, 86400);  % time in current day    
    
    % Baseline Schedule.
    % 1-sun, 2-mon, ..., 7-sat
    if kStep ==1    % change if this changed in idf
        dayWeek = 'Weekdays';
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
        SP(1) = ClgSPrand(kStep);
        SP(3) = SupplyAirSPrand(kStep);
        SP(4) = ChwSPrand(kStep);
    end
    
    if strcmp(data_type, 'ramped') && kStep>1
        ChwPrev = inputs(4,kStep-1);
        SP(4) = max(ChwMin,ChwPrev-rampVal)+(min(ChwMax,ChwPrev+rampVal)-max(ChwMin,ChwPrev-rampVal))*rand;
        if strcmp(sample_type, 'prbs')
            SP(4) = max(ChwMin,ChwPrev-rampVal)*(ChwSPprbs(kStep)<0)+min(ChwMax,ChwPrev+rampVal)*(ChwSPprbs(kStep)>0);
        end
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
toc;

cd('../../')

disp(['Stopped with flag ' num2str(flag)]);

%% random sampling post processing

data.Ambient = outputs(6,:)';
data.Humidity =  outputs(7,:)';
data.TotalLoad =  outputs(9,:)';
data.Month = outputs(1,:)';
data.DOM = outputs(2,:)';
data.TOD =  outputs(3,:)';
data.DOW = outputs(4,:)';
data.ClgSP = inputs(1,:)';
data.LgtSP = inputs(2,:)';
data.SupplyAirSP = inputs(3,:)';
data.ChwSP = inputs(4,:)';


ctrl_vars_all = {'ClgSP', 'LgtSP', 'SupplyAirSP', 'ChwSP'};
ctrl_idx = [1, 2, 3, 4];
for idn = 1:numel(ctrl_vars)
    figure('Name', 'random sampling'); grid on;
    plot(inputs(ctrl_idx(strcmp(ctrl_vars{idn},ctrl_vars_all)),:), 'LineWidth', 2)
    ylabel(ctrl_vars{idn})
    xlabel('sample number')
end

%% Save results

saveStr = sprintf('random_%s_%s_%dramped_%dinput_%dday_%04d%02d%02d_%02d%02d.mat',...
    building, sample_type, rampVal, numel(ctrl_vars), SimDays, YY, MM, DD, HH, MINS);
save(fullfile('data', saveStr),'-struct','data');
