% scratch_functions.m
% called by make_dev script

% -- Search for local functions --  
function fh = scratch_functions
    fh = localfunctions;
end

% -------- SCRATCH ---------------
% Put mess here:
function p = plot_poet_data(poet_data)

scale_pd = max(abs(poet_data),[],1);

yax_des= [1 3.1 5.2 6.3 7.4];

figure, hold on,
    plot(poet_data(:,1)/scale_pd(1)+yax_des(5), 'k-')
    plot(poet_data(:,2)/scale_pd(2)+yax_des(4), 'k-')
    plot(poet_data(:,3)/scale_pd(3)+yax_des(3), 'k-')
    plot(poet_data(:,4)/scale_pd(4)+yax_des(2), 'k-')
    plot(poet_data(:,5)/scale_pd(5)+yax_des(1), 'k-')
    
    ytick('off'); ylim([-0.1 8.5]);

names = {'Z';'Y';'X';'RF';'ADC'};
set(gca,'ytick',yax_des,'yticklabel',names)

end

function MRI_time(time_sec)
min = floor(time_sec/60);
sec = rem(time_sec,60);
disp([num2str(min) 'm ' num2str(sec) 's']);
end

function y = is_even(x)
    y=~rem(x,2);
end

function res = grs2rgb(img, map)
%%Convert grayscale images to RGB using specified colormap.
%	IMG is the grayscale image. Must be specified as a name of the image 
%	including the directory, or the matrix.
%	MAP is the M-by-3 matrix of colors.
%
%	RES = GRS2RGB(IMG) produces the RGB image RES from the grayscale image IMG 
%	using the colormap HOT with 64 colors.
%
%	RES = GRS2RGB(IMG,MAP) produces the RGB image RES from the grayscale image 
%	IMG using the colormap matrix MAP. MAP must contain 3 columns for Red, 
%	Green, and Blue components.  
%
%	Example 1:
%	open 'image.tif';	
%	res = grs2rgb(image);
%
%	Example 2:
%	cmap = colormap(summer);
% 	res = grs2rgb('image.tif',cmap);
%
% 	See also COLORMAP, HOT
%
%	Written by 
%	Valeriy R. Korostyshevskiy, PhD
%	Georgetown University Medical Center
%	Washington, D.C.
%	December 2006
%
% 	vrk@georgetown.edu
% Check the arguments
if nargin<1
	error('grs2rgb:missingImage','Specify the name or the matrix of the image');
end;
if ~exist('map','var') || isempty(map)
	map = hot(64);
end;
[l,w] = size(map);
if w~=3
	error('grs2rgb:wrongColormap','Colormap matrix must contain 3 columns');
end;
if ischar(img)
	a = imread(img);
elseif isnumeric(img)
	a = img;
else
	error('grs2rgb:wrongImageFormat','Image format: must be name or matrix');
end;
% Calculate the indices of the colormap matrix
a = double(a);
a(a==0) = 1; % Needed to produce nonzero index of the colormap matrix
ci = ceil(l*a/max(a(:))); 
% Colors in the new image
[il,iw] = size(a);
r = zeros(il,iw); 
g = zeros(il,iw);
b = zeros(il,iw);
r(:) = map(ci,1);
g(:) = map(ci,2);
b(:) = map(ci,3);
% New image
res = zeros(il,iw,3);
res(:,:,1) = r; 
res(:,:,2) = g; 
res(:,:,3) = b;
end

function grab_study(datPath)
% requires specific data organisation
% pastes all data in to the command window

if nargin < 1
    if ismac
        datPath = '/Volumes/DIRHome/';
    else
        datPath = uigetdir('\\hl-share.nhlbi.nih.gov\DIRHome\RamasawmyR\Scan Data');
    end
    datPath = [datPath filesep];
end

dir_dp = dir(datPath);

study_cell = cell(length(dir_dp), 1);

for i = 1:length(dir_dp)
    temp = dir_dp(i).name;
    if regexp(temp, '.dat')
        study_cell{i} = [['[] = recon_([' 'dirPath ''' temp(1:end-4) '.h5''' ']'] ',' [ '[' 'noisePath ''noise_' temp(1:end-4) '.h5''' ']);' ] ];
    end
end
study_cell{1} = ['dirPath = ''' datPath 'h5' filesep ''';'];
study_cell{2} = ['noisePath = ''' datPath 'noise' filesep ''';'];

fprintf(1, '%s \n', study_cell{:});
fprintf(1, '\n');

end

function z = quickmean(x,y)
if nargin == 1
    y = x;
end
roi = x > mean(x(:));
z = mean(y(roi));
end

function y = nrr(x)
y = x./mean(x(:));
end

% #####################
% Hargreaves...

function [Afp,Bfp]=freeprecess(T,T1,T2,df)
%
%	Function simulates free precession and decay
%	over a time interval T, given relaxation times T1 and T2
%	and off-resonance df.  Times in ms, off-resonance in Hz.

phi = 2*pi*df*T/1000;	% Resonant precession, radians.
E1 = exp(-T/T1);	
E2 = exp(-T/T2);

Afp = [E2 0 0;0 E2 0;0 0 E1]*zrot(phi);
Bfp = [0 0 1-E1]';

end

function Rz=zrot(phi)

Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];

end

function Rx=xrot(phi)

Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];

end

function Ry=yrot(phi)

Ry = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];

end