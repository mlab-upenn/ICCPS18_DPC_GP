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
    
RESULTS = struct('Time', [],...
    'TotalLoad', [], 'HVACLoad', [], 'OtherLoad', [],...
    'Ambient', [], 'Humidity', []);

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

DoWVar = findEPResultVar('EMS', 'currentDayOfWeek', allObjs, allVars);
HolidayVar = findEPResultVar('EMS', 'currentHoliday', allObjs, allVars);
assert(~isempty(DoWVar) && ~isempty(HolidayVar), 'Day data is missing.');

% The PMV variables
PMVVars = findEPResultVar('', 'Zone Thermal Comfort Fanger Model PMV', allObjs, allVars);
nPMVVars = length(PMVVars);

% The setpoints
ChwSPVar = findEPResultVar('CW-LOOP-TEMP-SCHEDULE', 'Schedule Value', allObjs, allVars);
assert(~isempty(ChwSPVar), 'Chilled water SP is missing.');

ClgSPVar = findEPResultVar('CLGSETP_SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(ClgSPVar), 'Cooling SP is missing.');

LgtSPVar = findEPResultVar('BLDG_LIGHT_SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(LgtSPVar), 'BLDG_LIGHT_SCH is missing.');

SupplyAirSPVar = findEPResultVar('SEASONAL-RESET-SUPPLY-AIR-TEMP-SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(SupplyAirSPVar), 'Seasonal-Reset-Supply-Air-Temp-Sch is missing.');

% Load the data for the interested variables
[~, DATA, TS] = mlepLoadEPResults(...
    fname,...
    [LoadVar, ElecLoadVar, HVACLoadVar,...
    AMBIENTVar, HumidVar,...
    DoWVar, HolidayVar,...
    LgtSPVar, ChwSPVar, ClgSPVar, SupplyAirSPVar,...
    PMVVars]);

nData = size(DATA,1);

RESULTS.TotalLoad = DATA(:, 1);
RESULTS.OtherLoad = DATA(:, 2);
RESULTS.HVACLoad = DATA(:, 3);

RESULTS.Ambient = DATA(:, 4);
RESULTS.Humidity = DATA(:, 5);

RESULTS.Time = [TS, DATA(:, 6:7)];

RESULTS.LgtSP = DATA(:, 8);
RESULTS.ChwSP = DATA(:, 9);
RESULTS.ClgSP = DATA(:, 10);
RESULTS.SupplyAirSP = DATA(:, 11);

PMVbegin = 12;

RESULTS.PMV = struct('zone', allObjs(PMVVars),...
    'PMV', mat2cell(DATA(:, PMVbegin:(PMVbegin-1+nPMVVars)), nData, ones(nPMVVars,1)));

end
