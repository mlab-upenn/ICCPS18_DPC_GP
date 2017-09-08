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
SimDays = 30+31+31;
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
data_type = 'rulebased';

if strcmp(data_type, 'unconstrained')
    SPmin = repmat([22, 21, 24, 19, 22, 21, 12, 3.7], size(DRS,1),1);
    SPmax = repmat([32, 21, 32, 19, 26, 21, 14, 9.7], size(DRS,1),1);
    DRS = SPmin + (SPmax-SPmin).*DRS;
elseif strcmp(data_type, 'constrained')
    SPmin = repmat([22, 21, 24, 19, 22, 21, 12, 3.7], size(DRS,1),1);
    SPmax = repmat([27, 21, 32, 19, 26, 21, 14, 9.7], size(DRS,1),1);
    DRS = SPmin + (SPmax-SPmin).*DRS;
elseif strcmp(data_type, 'rulebased')
    SPbefore7 = repmat([30, 21, 30, 19, 24, 21, 13, 6.7], 7*4,1);
    SPafter7 = repmat([24, 21, 26, 19, 24, 21, 13, 6.7], 17*4,1);
    DRS = repmat([SPbefore7; SPafter7], SimDays, 1);
end



while kStep <= MAXSTEPS    

    % BEGIN Compute next set-points
    dayTime = mod(eptime, 86400);  % time in current day

    SP = DRS(kStep,:);
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

% save_results(data_type, 'train', 1, 1)
