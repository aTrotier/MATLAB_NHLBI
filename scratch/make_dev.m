% make_dev.m 
%
% creates a dev "class" from scratch_functions.m
% functions then can be called by dev.<function_name>(X)
% R Ramasawmy Jan 2019
% This and scratch_functions.m can be combined in to a single file from
% 2016b.

% Get a cell array of function handles
disp('Calling "scratch_functions.m"');
fh = scratch_functions;

% And use the function names as substructures..
scratch_list = cell(1,length(fh));
for i =1:length(fh)
    scratch_list{i} = func2str(fh{i});
end

% ### DEV ###
disp('Making "dev" class');
dev = cell2struct(fh, scratch_list); % One can call this anything they like, but this is nice and short. 

% Tidy up
clear fh i scratch_list; 