% scratch_functions.m
% called by make_dev script

% -- Search for local functions --  
function fh = scratch_functions
    fh = localfunctions;
end

% -------- SCRATCH ---------------
% Put mess here:

function y = add1(x)
y = x + 1;
end

function y = is_even(x)
    y=~rem(x,2);
end

function y = xpi(x)
    y=x*pi;
end