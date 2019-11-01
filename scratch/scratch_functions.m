% scratch_functions.m
% called by make_dev script

% -- Search for local functions --  
function fh = scratch_functions
    fh = localfunctions;
end

% -------- SCRATCH ---------------
function edit_toolbox
edit(mfilename)
end

% Put mess here:
function [ROI] = rroipoly(IM, numROIs)
if nargin < 2
    numROIs = 1;
end

numSlices = size(IM,3);
IM = IM./mean(IM(:));

ROI = zeros(size(IM));

figure,
for i = 1:numSlices
    imshow(IM(:,:,i),[0 4]);
    
    for j = 1:numROIs
        set(gcf, 'Name', ['ROI ' num2str(j) '/' num2str(numROIs) ', Slice ' num2str(i) '/' num2str(numSlices)]);
        ROI(:,:,i) = ROI(:,:,i) + roipoly;
        
    end
    hold on;
    contour(ROI(:,:,i), 'g');
    hold off;
    pause(1);
end

end

function plot_rois(xi,x)
montage_RR(xi);
hold on; 
contour(x.GM_roi, 'b'); 
contour(x.WM_roi, 'y');
end

function subplot_link(data)
num_wvd = size(data, 2);
sp1 = floor(sqrt(num_wvd));
sp2 = ceil(num_wvd/sp1);
figure,
hl = [];
for i = 1:num_wvd
    h.(matlab.lang.makeValidName(['x' num2str(i)])) = subplot(sp1, sp2, i);
    plot(data(:,i))
    hl = [hl, h.(matlab.lang.makeValidName(['x' num2str(i)]))];
end
linkaxes(hl, 'x');
end

function [output] = get_rois(xi,x)
output = [mean(xi(x.WM_roi)) std(xi(x.WM_roi)) mean(xi(x.GM_roi)) std(xi(x.GM_roi))];
end

function implay_flow(M,P, gif_filename)
if ndims(M) == 4
    M = M(:,:,:,1);
end
M = M./(0.5*max(M(:)));

if (~isempty(P))
P = (min(P(:))<0)*0.5 + P./(max(P(:)) - min(P(:)));
end

if (size(M,2) == size(M,1))
    implay_RR([M P],[0 1]);
else
    implay_RR([M; P],[0 1]);
end

if nargin > 2
    %     gif_filename = 'testAnimated.gif'; % Specify the output file name
    disp(['Writing gif to ' gif_filename]);
    for idx = 1:size(M,3)
        if (~isempty(P))
             if (size(M,2) == size(M,1))
                [A,map] = rgb2ind(repmat([M(:,:,idx) P(:,:,idx)],[1 1 3]),256);
             else
                 [A,map] = rgb2ind(repmat([M(:,:,idx); P(:,:,idx)],[1 1 3]),256);
             end
            
        else
            [A,map] = rgb2ind(repmat([M(:,:,idx)],[1 1 3]),256);
        end
        %         [A] = [M(:,:,idx) P(:,:,idx)];
        if idx == 1
            imwrite(A,map,gif_filename,'gif','LoopCount',Inf,'DelayTime',.12); % 0.25 s
        else
            imwrite(A,map,gif_filename,'gif','WriteMode','append','DelayTime',.25);
        end
    end
end

end

function R_info = RR_rot_m(R)
% https://en.wikipedia.org/wiki/Rotation_matrix
symmetric_flag = R - R';

if sum(symmetric_flag) == 0
    R_info.axis = eig(R); % ? Not sure what to do here
else
    R_info.axis = [R(8) - R(6); R(3) - R(7) ; R(4) - R(2)];
end

R_info.axis = R_info.axis/norm(R_info.axis);
R_info.angle = acosd((sum(diag(R)) - 1)/2);

end

function plot_poet_data(poet_data)

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

function plot_rf_spectrum(freq, coil_label)
samples = size(freq,1);
channels = size(freq,3);

sp1 = floor(sqrt(channels));
sp2 = ceil(channels/sp1);

figure,
for i = 1:channels
    subplot(sp1,sp2, i);
    shadedErrorBar(1:samples, mean(freq(:,:,i),2)', std(freq(:,:,i),[],2)');
    if nargin > 1
        title(coil_label(i,2), 'Interpreter', 'none');
    end
end

end

function G0_2 = maxwell_phase_grad(x, y, interleaves)

grads = zeros(length(x), interleaves, 2);

for i = 1:interleaves
    rot = (i-1)*(2*pi/interleaves);
    grads(:,i,2)  = ( x *cos(rot) + y *sin(rot));
    grads(:,i,1)  = (-x *sin(rot) + y *cos(rot));
end

G0_2 = sum(grads.^2,1);

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

[~, current_dir] = fileparts(datPath(1:end-1));
header_cell{1} = ['% 20' current_dir];
header_cell{2} = [' ' ];
header_cell{3} = ['% save(''2019 studies\d20' current_dir '.mat'')'];
header_cell{4} = ['% load(''2019 studies\d20' current_dir '.mat'')'];
header_cell{5} = ['whos'];
header_cell{6} = [' '];
header_cell{7} = ['%%'];
header_cell{8} = [' '];
fprintf(1, '%s \n', header_cell{:});

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

function grab_dicom_study(dicomPath)
% 
% pastes all data in to the command window

if nargin < 1
    if ismac
        dicomPath = '/Volumes/DIRHome/';
    else
        dicomPath = uigetdir('\\hl-share.nhlbi.nih.gov\DIRHome\RamasawmyR\Scan Data');
    end
    dicomPath = [dicomPath filesep];
else
    if ~(regexp(dicomPath(end), filesep))
        dicomPath = [dicomPath filesep];
    end
end

% [fp, current_dir] = fileparts(dicomPath(1:end-1));
% header_cell{1} = ['% 20' current_dir];
% header_cell{2} = [' ' ];
% header_cell{3} = ['% save(''2019 studies\d20' current_dir '.mat'')'];
% header_cell{4} = ['% load(''2019 studies\d20' current_dir '.mat'')'];
% header_cell{5} = ['whos'];
% header_cell{6} = [' '];
% header_cell{7} = ['%%'];
% header_cell{8} = [' '];
% fprintf(1, '%s \n', header_cell{:});

dir_dp = dir(dicomPath);

% first check if the dicoms have been sorted
dir_check = [dir_dp(3:end).isdir];
if round(mean(dir_check)) == 0
    % dicoms need to be sorted
    dicom_sort_folder(dicomPath);
    
    % now grab sorted folder
    dir_dp = dir(dicomPath);
    
end

%
study_cell = cell(length(dir_dp ), 1);

for i = 1:length(dir_dp)
    temp = dir_dp(i).name;
    study_cell{i} = ['[s' temp ', info_' temp '] = dicom_load_scan([dirPath ''' temp ''']);' ];
    
end
study_cell{1} = ['dirPath = ''' dicomPath ''';'];
study_cell{2} = [''];

fprintf(1, '%s \n', study_cell{:});
fprintf(1, '\n');

end

function grab_mat_study(matPath)
% requires specific data organisation
% pastes all data in to the command window

if nargin < 1
    if ismac
        matPath = '/Volumes/DIRHome/';
    else
        matPath = uigetdir('\\hl-share.nhlbi.nih.gov\DIRHome\RamasawmyR\Scan Data');
    end
    matPath = [matPath filesep];
end

[fp, current_dir] = fileparts(matPath(1:end-1));
header_cell{1} = ['% 20' current_dir];
header_cell{2} = [' ' ];
header_cell{3} = ['% save(''2019 studies\d20' current_dir '.mat'')'];
header_cell{4} = ['% load(''2019 studies\d20' current_dir '.mat'')'];
header_cell{5} = ['whos'];
header_cell{6} = [' '];
header_cell{7} = ['%%'];
header_cell{8} = [' '];
fprintf(1, '%s \n', header_cell{:});

dir_dp = dir(matPath);

study_cell = cell(length(dir_dp )+ 1, 1);

for i = 1:length(dir_dp)
    temp = dir_dp(i).name;
    if regexp(temp, '.mat')
        study_cell{i} = [['load([dirPath ''' temp ''']);'] ' ' [ 'XXX = imrec;' ] ];
    end
end
study_cell{1} = ['dirPath = ''' matPath ''';'];
study_cell{2} = [''];
study_cell{end} = ['clear imrec'];

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