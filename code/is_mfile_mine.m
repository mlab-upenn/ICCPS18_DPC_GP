function b = is_mfile_mine(fname, mydir)
% Check if a given M-file fname is in a given directory mydir.
% Returns 0 if file doesn't exist, 1 if file exists and is in mydir, 2 if
% file exists but not in mydir.
%
% Example: is_mfile_mine('gp.m', '/path/to/my/dir')

b = 0;
if exist(fname, 'file')
    % Find the path to the existing file
    mytempdir = fileparts(which(fname));
    if strcmp(mytempdir, mydir)
        b = 1;
    else
        b = 2;
    end
end
end