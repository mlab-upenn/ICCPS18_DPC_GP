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

DoWVar = findEPResultVar('EMS', 'currentDayOfWeek', allObjs, allVars);
HolidayVar = findEPResultVar('EMS', 'currentHoliday', allObjs, allVars);
assert(~isempty(DoWVar) && ~isempty(HolidayVar), 'Day data is missing.');
ToDVar = findEPResultVar('EMS', 'currentTimeOfDay', allObjs, allVars);
assert(~isempty(ToDVar), 'Time of day data is missing.');

% control variables
ClgSPVar = findEPResultVar('CLGSETP_SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(ClgSPVar), 'Cooling set point data is missing.');

KitchenClgSPVar = findEPResultVar('KITCHEN_CLGSETP_SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(KitchenClgSPVar), 'Kitchen cooling set point data is missing.');

GuestClgSPVar = findEPResultVar('GUEST_CLGSETP_SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(GuestClgSPVar), 'Guest cooling set point data is missing.');

HtgSPVar = findEPResultVar('HTGSETP_SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(HtgSPVar), 'Heating set point data is missing.');

KitchenHtgSPVar = findEPResultVar('KITCHEN_HTGSETP_SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(KitchenHtgSPVar), 'Kitchen heating set point data is missing.');

GuestHtgSPVar = findEPResultVar('GUEST_HTGSETP_SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(GuestHtgSPVar), 'Guest heating set point data is missing.');

SeasonSATVar = findEPResultVar('SEASONAL-RESET-SUPPLY-AIR-TEMP-SCH', 'Schedule Value', allObjs, allVars);
assert(~isempty(SeasonSATVar), 'SAP data is missing.');

ChWaterVar = findEPResultVar('CW-LOOP-TEMP-SCHEDULE', 'Schedule Value', allObjs, allVars);
assert(~isempty(ChWaterVar), 'Chilled water data is missing.');

% Load the data for the interested variables
[~, DATA, TS] = mlepLoadEPResults(...
    fname,...
    [LoadVar, ElecLoadVar, HVACLoadVar,...
    AMBIENTVar, HumidVar,...
    ToDVar, DoWVar,...
    ClgSPVar, KitchenClgSPVar, GuestClgSPVar,...
    HtgSPVar, KitchenHtgSPVar, GuestHtgSPVar,...
    SeasonSATVar, ChWaterVar]);

RESULTS.TotalLoad = DATA(:, 1);
RESULTS.OtherLoad = DATA(:, 2);
RESULTS.HVACLoad = DATA(:, 3);

RESULTS.Ambient = DATA(:, 4);
RESULTS.Humidity = DATA(:, 5);

% RESULTS.Time = [TS, DATA(:, 6:7)];
RESULTS.DOW = DATA(:,7);
RESULTS.TOD = DATA(:,6); % floor((60*TS(:,3)+TS(:,4))/15)+1;

RESULTS.ClgSP = DATA(:, 8);
RESULTS.KitchenClgSP = DATA(:, 9);
RESULTS.GuestClgSP = DATA(:, 10);
RESULTS.HtgSP = DATA(:, 11);
RESULTS.KitchenHtgSP = DATA(:, 12);
RESULTS.GuestHtgSP = DATA(:, 13);
RESULTS.SupplyAirSP = DATA(:, 14);
RESULTS.ChwSP = DATA(:, 15);

end
