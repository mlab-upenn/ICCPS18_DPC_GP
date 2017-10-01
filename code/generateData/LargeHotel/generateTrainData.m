%% Create an mlepProcess instance and configure it

file = 'trainBuilding';

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
SimDays = 30;
deltaT = (60/EPTimeStep)*60;
kStep = 1;  % current simulation step
MAXSTEPS = SimDays*24*EPTimeStep;  % max simulation time = 7 days

% variables for plotting:
outputs = nan(9,MAXSTEPS);
inputs = nan(8,MAXSTEPS);

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
    DRS = [DRS;lhsdesign(1*96,8)];
end

% data_type = 'unconstrained';
% data_type = 'constrained';
% data_type = 'nominal';
data_type = 'ramped';
rampVal = 2;

if strcmp(data_type, 'unconstrained') || strcmp(data_type, 'ramped')
    SPmin = repmat([22, 21, 24, 19, 22, 21, 12, 3.7], size(DRS,1),1);
    SPmax = repmat([32, 21, 32, 19, 26, 21, 14, 9.7], size(DRS,1),1);
    DRS = SPmin + (SPmax-SPmin).*DRS;
elseif strcmp(data_type, 'constrained')
    SPmin = repmat([22, 21, 24, 19, 23, 21, 12.5, 5.2], size(DRS,1),1);
    SPmax = repmat([32, 21, 32, 19, 25, 21, 13.5, 8.2], size(DRS,1),1);
    DRS = SPmin + (SPmax-SPmin).*DRS;
end



while kStep <= MAXSTEPS    

    % BEGIN Compute next set-points
    dayTime = mod(eptime, 86400);  % time in current day
    
    % Baseline Schedule.
    if(dayTime <= 7*3600)
        
        clgstp = 30;
        htgstp = 16;
        kitclgstp = 30;
        kithtgstp = 16;
        guestclgstp = 24;
        guesthtgstp = 21;
        sat = 13;
        cwstp = 6.7;
        
        SP = [clgstp, htgstp, kitclgstp, kithtgstp, guestclgstp, guesthtgstp, sat, cwstp];

    else
        
        clgstp = 24;
        htgstp = 21;
        kitclgstp = 26;
        kithtgstp = 19;
        guestclgstp = 24;
        guesthtgstp = 21;
        sat = 13;
        cwstp = 6.7;
        
        SP = [clgstp, htgstp, kitclgstp, kithtgstp, guestclgstp, guesthtgstp, sat, cwstp];
        
    end
    
    if ~strcmp(data_type, 'nominal')
        SP(5) = DRS(kStep,5);
        SP(7) = DRS(kStep,7);
        SP(8) = DRS(kStep,8);
    end
    
    if strcmp(data_type, 'ramped') && kStep>1
        ChwMin = 3.7;
        ChwMax = 9.7;
        ChwPrev = inputs(8,kStep-1);
        SP(8) = max(ChwMin,ChwPrev-rampVal)+(min(ChwMax,ChwPrev+rampVal)-max(ChwMin,ChwPrev-rampVal))*rand;
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


figure
plot(1:MAXSTEPS,outputs(1,:));
figure
plot(1:MAXSTEPS,outputs(2,:));
figure
plot(1:MAXSTEPS,outputs(3,:));
figure
plot(1:MAXSTEPS,outputs(9,:));
