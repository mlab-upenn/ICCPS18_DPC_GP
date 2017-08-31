% Setup the toolboxes used by the ICCPS'18 code, for Matlab.

me = mfilename;                                            % what is my filename
disp(mfilename);
disp(me);
mydir = which(me); mydir = mydir(1:end-3-numel(me));        % where am I located

% GPML
if exist('gp', 'file')
    warning('GPML may already exist. You should remove that version so that the correct version will be used.');
end
addpath(fullfile(mydir, 'dependencies', 'gpml-old'))
addpath(fullfile(mydir, 'dependencies', 'gpml-old', 'cov'))
addpath(fullfile(mydir, 'dependencies', 'gpml-old', 'doc'))
addpath(fullfile(mydir, 'dependencies', 'gpml-old', 'inf'))
addpath(fullfile(mydir, 'dependencies', 'gpml-old', 'lik'))
addpath(fullfile(mydir, 'dependencies', 'gpml-old', 'mean'))
addpath(fullfile(mydir, 'dependencies', 'gpml-old', 'util'))
addpath(fullfile(mydir, 'dependencies', 'gpml-old', 'prior'))
addpath(fullfile(mydir, 'dependencies', 'gpml-old', 'util', 'minfunc'))
addpath(fullfile(mydir, 'dependencies', 'gpml-old', 'util', 'minfunc', 'compiled'))

% NEXTGP
addpath(genpath(fullfile(mydir, 'dependencies', 'nextgp')));

% active_gp_hyperlearning
addpath(fullfile(mydir, 'dependencies', 'active_gp_hyperlearning'));

% gpml_extensions
addpath(genpath(fullfile(mydir, 'dependencies', 'gpml_extensions')));

% mgp
addpath(fullfile(mydir, 'dependencies', 'mgp'));

clear me mydir