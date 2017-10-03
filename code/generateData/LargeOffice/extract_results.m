function RESULTS = extract_results(fname)
%EXTRACT_RESULTS Summary of this function goes here
%   Detailed explanation goes here

% pause to make sure the E+ file was written
ktrial = 0;
while ktrial < 20
    if exist(fname, 'file') > 0
        break;
    else
        pause(1);
    end
    ktrial = ktrial + 1;
end
assert(exist(fname, 'file') > 0, 'File %s does not exist.', fname);

fprintf('Loading results from file %s...\n', fname);

% Load the results
% Results for each simulation
% - *Load: vector of actual loads for the simulation duration
% - Ambient air temperature, Humidity
% - Time
    
% RESULTS = struct('Time', [],...
%     'TotalLoad', [], 'HVACLoad', [], 'OtherLoad', [],...
%     'Ambient', [], 'Humidity', []);

% Read the names of the variables from output file
VARS = mlepLoadEPResults(fname, 'vars');
allObjs = {VARS.object};
allVars = {VARS.name};

% Find the variables we are intersted in
LoadVar = findEPResultVar('Whole Building', 'Facility Total Electric Demand Power', allObjs, allVars);
assert(~isempty(LoadVar), 'Total load data is missing.');

ElecLoadVar = findEPResultVar('Whole Building', 'Facility Total Building Electric Demand Power', allObjs, allVars);
assert(~isempty(ElecLoadVar), 'Electrical load data is missing.');

HVACLoadVar = findEPResultVar('Whole Building', 'Facility Total HVAC Electric Demand Power', allObjs, allVars);
assert(~isempty(HVACLoadVar), 'HVAC load data is missing.');

AMBIENTVar = findEPResultVar('ENVIRONMENT', 'Site Outdoor Air Drybulb Temperature', allObjs, allVars);
assert(~isempty(AMBIENTVar), 'Ambient temp data is missing.');

HumidVar = findEPResultVar('ENVIRONMENT', 'Site Outdoor Air Relative Humidity', allObjs, allVars);
assert(~isempty(HumidVar), 'Humidity data is missing.');

MonthVar = findEPResultVar('EMS', 'currentMonth', allObjs, allVars);
DoMVar = findEPResultVar('EMS', 'currentDayOfMonth', allObjs, allVars);
assert(~isempty(MonthVar) && ~isempty(DoMVar), 'Month data is missing.');

DoWVar = findEPResultVar('EMS', 'currentDayOfWeek', allObjs, allVars);
HolidayVar = findEPResultVar('EMS', 'currentHoliday', allObjs, allVars);
assert(~isempty(DoWVar) && ~isempty(HolidayVar), 'Day data is missing.');

ToDVar = findEPResultVar('EMS', 'currentTimeOfDay', allObjs, allVars);
assert(~isempty(ToDVar), 'Time of day data is missing.');

% control variables
ClgSPVar = findEPResultVar('CLGSETP_SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(ClgSPVar), 'Cooling set point data is missing.');

LgtSPVar = findEPResultVar('BLDG_LIGHT_SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(LgtSPVar), 'Lighting set point data is missing.');

SeasonSATVar = findEPResultVar('Seasonal-Reset-Supply-Air-Temp-Sch', 'Schedule Value', allObjs, allVars);
assert(~isempty(SeasonSATVar), 'SAP data is missing.');

ChWaterVar = findEPResultVar('CW-Loop-Temp-Schedule', 'Schedule Value', allObjs, allVars);
assert(~isempty(ChWaterVar), 'Chilled water data is missing.');

% Load the data for the interested variables
[~, DATA, TS] = mlepLoadEPResults(...
    fname,...
    [LoadVar, ElecLoadVar, HVACLoadVar,...
    AMBIENTVar, HumidVar,...
    ToDVar, DoWVar,...
    MonthVar, DoMVar, ...
    ClgSPVar, LgtSPVar, ...
    SeasonSATVar, ChWaterVar, ...
    ]);

RESULTS.TotalLoad = DATA(:, 1);
RESULTS.OtherLoad = DATA(:, 2);
RESULTS.HVACLoad = DATA(:, 3);

RESULTS.Ambient = DATA(:, 4);
RESULTS.Humidity = DATA(:, 5);

% RESULTS.Time = [TS, DATA(:, 6:7)];
RESULTS.DOW = DATA(:,7);
RESULTS.TOD = DATA(:,6); % floor((60*TS(:,3)+TS(:,4))/15)+1;
RESULTS.Month = DATA(:,8);
RESULTS.DOM = DATA(:,9);

RESULTS.ClgSP = DATA(:, 10);
RESULTS.LgtSP = DATA(:, 11);
RESULTS.SupplyAirSP = DATA(:, 12);
RESULTS.ChwSP = DATA(:, 13);

end
