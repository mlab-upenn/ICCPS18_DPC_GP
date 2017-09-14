%% Create an mlepProcess instance and configure it

file = 'testBuilding';

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
SimDays = 31;
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
