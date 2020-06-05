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
% output = [mean(xi(x.WM_roi)) std(xi(x.WM_roi)) mean(xi(x.GM_roi)) std(xi(x.GM_roi))];
output = [mean(xi(x)) std(xi(x)) ];
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
            imwrite(A,map,gif_filename,'gif','WriteMode','append','DelayTime',.12);
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

function flow_info = rtfm_to_flow_info(RTFM_structure)

flow_info.NominalInterval = RTFM_structure.RB_ECG.NominalInterval; 

if regexpi(RTFM_structure.mrd_header.userParameters.userParameterLong.name, 'venc')
    flow_info.venc = RTFM_structure.mrd_header.userParameters.userParameterLong.value;
end

flow_info.CardiacNumberOfImages = size(RTFM_structure.RB_ECG.PHA,3);

flow_info.PixelSpacing = [RTFM_structure.mrd_header.encoding.encodedSpace.fieldOfView_mm.x/RTFM_structure.mrd_header.encoding.encodedSpace.matrixSize.x ...
                          RTFM_structure.mrd_header.encoding.encodedSpace.fieldOfView_mm.y/RTFM_structure.mrd_header.encoding.encodedSpace.matrixSize.y];

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
% local config?

username=getenv('USERNAME');

if nargin < 1
    if regexp(username, 'ramasawmyr')
        if ismac
            datPath = '/Volumes/DIRHome/';
        else
            datPath = uigetdir('\\hl-share.nhlbi.nih.gov\DIRHome\RamasawmyR\Scan Data');
        end
    else
       datPath = pwd;
    end
    datPath = [datPath filesep];
end

[~, current_dir] = fileparts(datPath(1:end-1));
header_cell{1} = ['% s' current_dir];
header_cell{2} = [' ' ];
if regexp(username, 'ramasawmyr')
    header_cell{3} = ['% save(''H:\Backup\MATLAB R2016a\2020_studies\d' current_dir '.mat'')'];
    header_cell{4} = ['% load(''H:\Backup\MATLAB R2016a\2020_studies\d' current_dir '.mat'')'];
else
    header_cell{3} = ['% save(''d' current_dir '.mat'')'];
    header_cell{4} = ['% load(''d' current_dir '.mat'')'];
end

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
        dicomPath = uigetdir('/Volumes/');
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
    study_cell{i} = ['[s' temp ', info_' temp '] = dicom_load_scan([ddirPath ''' temp ''']);' ];
    
end
study_cell{1} = ['ddirPath = ''' dicomPath ''';'];
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
    roi = x > mean(x(:));
    z = mean(y(roi));
else
    z = mean(x(find(y)));
end
end

function y = nrr(x)
y = x./mean(x(:));
end

function [traj_order] = skip_step(num_int, fib_flag)
if nargin < 2
    fib_flag = 0;
end

traj_order = zeros(1, num_int);
if fib_flag == 0
    x1 = 1:floor(num_int/2);
    x2 = fliplr((floor(num_int/2)+1):num_int);
    
    if mod(num_int,2) == 0
        traj_order(1:2:end) = x1;
        traj_order(2:2:end) = x2;
    else
        traj_order(1:2:end) = x2;
        traj_order(2:2:end) = x1;
    end
else 
    fib_series = [3 5 8 13 21 34 55 89 144];
    fib_step = fib_series(find(num_int < fib_series, 1, 'first')-fib_flag); % 1, normal, 2, tiny step (better?)
    temp = 1 + mod([0:fib_step:(fib_step*num_int)], num_int);
    traj_order = temp(1:num_int);
    % check
    if ~length(unique(traj_order))==num_int
        warning('Fib series error');
    end
end

end

function [x,y] = radial_point(r, theta)
% avoiding tan theta ~= 90 degrees
% some k-space assumptions

if ( ((pi*(1/4)) < abs(theta)) && ( abs(theta) < (pi*(3/4))) )
    theta_orig = theta;
    theta = abs(theta - pi/2);
    sign = int8(theta <= pi/2); if sign == 0; sign = -1; end
    
    x = r*cos(theta)*sign;
    y = r*sin(theta);
    
else
    
    x = r*sin(theta);
    y = r*cos(theta);
end

end

function USC_cmap(ax)
if nargin < 1
    ax = gca;
end
set(ax, 'Colormap', [0 0 0.600000023841858;0.00116421573329717 0.00974264740943909 0.606250047683716;0.00232843146659434 0.0194852948188782 0.612500011920929;0.00349264731630683 0.0292279422283173 0.618750035762787;0.00465686293318868 0.0389705896377563 0.625;0.00582107855007052 0.0487132370471954 0.631250023841858;0.00698529463261366 0.0584558844566345 0.637500047683716;0.00814951024949551 0.0681985318660736 0.643750011920929;0.00931372586637735 0.0779411792755127 0.650000035762787;0.0104779414832592 0.0876838266849518 0.65625;0.011642157100141 0.0974264740943909 0.662500023841858;0.0128063727170229 0.10716912150383 0.668750047683716;0.0139705892652273 0.116911768913269 0.675000011920929;0.0151348048821092 0.126654416322708 0.681250035762787;0.016299020498991 0.136397063732147 0.6875;0.0174632351845503 0.146139711141586 0.693750023841858;0.0186274517327547 0.155882358551025 0.700000047683716;0.0197916682809591 0.165625005960464 0.706250011920929;0.0209558829665184 0.175367653369904 0.712500035762787;0.0221200995147228 0.185110300779343 0.71875;0.0232843142002821 0.194852948188782 0.725000023841858;0.0244485307484865 0.204595595598221 0.731249988079071;0.0256127454340458 0.21433824300766 0.737500011920929;0.0267769619822502 0.224080890417099 0.743750035762787;0.0279411785304546 0.233823537826538 0.75;0.0291053932160139 0.243566185235977 0.756250023841858;0.0302696097642183 0.253308832645416 0.762499988079071;0.0314338244497776 0.263051480054855 0.768750011920929;0.032598040997982 0.272794127464294 0.775000035762787;0.0337622575461864 0.282536774873734 0.78125;0.0349264703691006 0.292279422283173 0.787500023841858;0.036090686917305 0.302022069692612 0.793749988079071;0.0372549034655094 0.311764717102051 0.800000011920929;0.0384191200137138 0.32150736451149 0.806250035762787;0.0395833365619183 0.331250011920929 0.8125;0.0407475493848324 0.340992659330368 0.818750023841858;0.0419117659330368 0.350735306739807 0.824999988079071;0.0430759824812412 0.360477954149246 0.831250011920929;0.0442401990294456 0.370220601558685 0.837500035762787;0.0454044118523598 0.379963248968124 0.84375;0.0465686284005642 0.389705896377563 0.850000023841858;0.0477328449487686 0.399448543787003 0.856249988079071;0.048897061496973 0.409191191196442 0.862500011920929;0.0500612780451775 0.418933838605881 0.868750035762787;0.0512254908680916 0.42867648601532 0.875;0.052389707416296 0.438419133424759 0.881250023841858;0.0535539239645004 0.448161780834198 0.887499988079071;0.0547181405127048 0.457904428243637 0.893750011920929;0.0558823570609093 0.467647075653076 0.899999976158142;0.0570465698838234 0.477389723062515 0.90625;0.0582107864320278 0.487132370471954 0.912500023841858;0.0593750029802322 0.496875017881393 0.918749988079071;0.0605392195284367 0.506617665290833 0.925000011920929;0.0617034323513508 0.516360282897949 0.931249976158142;0.0628676488995552 0.526102960109711 0.9375;0.0640318617224693 0.535845637321472 0.943750023841858;0.0651960819959641 0.545588254928589 0.949999988079071;0.0663602948188782 0.555330872535706 0.956250011920929;0.0675245150923729 0.565073549747467 0.962499976158142;0.068688727915287 0.574816226959229 0.96875;0.0698529407382011 0.584558844566345 0.975000023841858;0.0710171610116959 0.594301462173462 0.981249988079071;0.07218137383461 0.604044139385223 0.987500011920929;0.0733455941081047 0.613786816596985 0.993749976158142;0.0745098069310188 0.623529434204102 1;0.0889705941081047 0.629411816596985 1;0.10343137383461 0.635294139385223 1;0.117892161011696 0.641176462173462 1;0.132352948188782 0.647058844566345 1;0.146813735365868 0.652941226959229 1;0.161274507641792 0.658823549747467 1;0.175735294818878 0.664705872535706 1;0.190196081995964 0.670588254928589 1;0.20465686917305 0.676470637321472 1;0.219117656350136 0.682352960109711 1;0.23357842862606 0.688235282897949 1;0.248039215803146 0.694117665290833 1;0.262499988079071 0.700000047683716 1;0.276960790157318 0.705882370471954 1;0.291421562433243 0.711764693260193 1;0.30588236451149 0.717647075653076 1;0.320343136787415 0.723529458045959 1;0.334803909063339 0.729411780834198 1;0.349264711141586 0.735294103622437 1;0.363725483417511 0.74117648601532 1;0.378186285495758 0.747058868408203 1;0.392647057771683 0.752941191196442 1;0.40710785984993 0.75882351398468 1;0.421568632125854 0.764705896377563 1;0.436029404401779 0.770588278770447 1;0.450490206480026 0.776470601558685 1;0.464950978755951 0.782352924346924 1;0.479411780834198 0.788235306739807 1;0.493872553110123 0.79411768913269 1;0.508333325386047 0.800000011920929 1;0.522794127464294 0.805882334709167 1;0.537254929542542 0.811764717102051 1;0.551715672016144 0.817647099494934 1;0.566176474094391 0.823529422283173 1;0.580637276172638 0.829411745071411 1;0.59509801864624 0.835294127464294 1;0.609558820724487 0.841176509857178 1;0.624019622802734 0.847058832645416 1;0.638480365276337 0.852941155433655 1;0.652941167354584 0.858823537826538 1;0.667401969432831 0.864705920219421 1;0.681862771511078 0.87058824300766 1;0.69632351398468 0.876470565795898 1;0.710784316062927 0.882352948188782 1;0.725245118141174 0.888235330581665 1;0.739705860614777 0.894117653369904 1;0.754166662693024 0.899999976158142 1;0.768627464771271 0.905882358551025 1;0.783088207244873 0.911764740943909 1;0.79754900932312 0.917647063732147 1;0.812009811401367 0.923529386520386 1;0.826470613479614 0.929411768913269 1;0.840931355953217 0.935294151306152 1;0.855392158031464 0.941176474094391 1;0.869852960109711 0.947058796882629 1;0.884313702583313 0.952941179275513 1;0.89877450466156 0.958823561668396 1;0.913235306739807 0.964705884456635 1;0.927696049213409 0.970588207244873 1;0.942156851291656 0.976470589637756 1;0.956617653369904 0.98235297203064 1;0.971078455448151 0.988235294818878 1;0.985539197921753 0.994117617607117 1;1 1 1;1 0.984375 0.984375;1 0.96875 0.96875;1 0.953125 0.953125;1 0.9375 0.9375;1 0.921875 0.921875;1 0.90625 0.90625;1 0.890625 0.890625;1 0.875 0.875;1 0.859375 0.859375;1 0.84375 0.84375;1 0.828125 0.828125;1 0.8125 0.8125;1 0.796875 0.796875;1 0.78125 0.78125;1 0.765625 0.765625;1 0.75 0.75;1 0.734375 0.734375;1 0.71875 0.71875;1 0.703125 0.703125;1 0.6875 0.6875;1 0.671875 0.671875;1 0.65625 0.65625;1 0.640625 0.640625;1 0.625 0.625;1 0.609375 0.609375;1 0.59375 0.59375;1 0.578125 0.578125;1 0.5625 0.5625;1 0.546875 0.546875;1 0.53125 0.53125;1 0.515625 0.515625;1 0.5 0.5;1 0.484375 0.484375;1 0.46875 0.46875;1 0.453125 0.453125;1 0.4375 0.4375;1 0.421875 0.421875;1 0.40625 0.40625;1 0.390625 0.390625;1 0.375 0.375;1 0.359375 0.359375;1 0.34375 0.34375;1 0.328125 0.328125;1 0.3125 0.3125;1 0.296875 0.296875;1 0.28125 0.28125;1 0.265625 0.265625;1 0.25 0.25;1 0.234375 0.234375;1 0.21875 0.21875;1 0.203125 0.203125;1 0.1875 0.1875;1 0.171875 0.171875;1 0.15625 0.15625;1 0.140625 0.140625;1 0.125 0.125;1 0.109375 0.109375;1 0.09375 0.09375;1 0.078125 0.078125;1 0.0625 0.0625;1 0.046875 0.046875;1 0.03125 0.03125;1 0.015625 0.015625;1 0 0;0.99421101808548 0.00124494242481887 0.00292561482638121;0.988422036170959 0.00248988484963775 0.00585122965276241;0.982633054256439 0.00373482750728726 0.00877684447914362;0.976844072341919 0.00497976969927549 0.0117024593055248;0.971055090427399 0.00622471235692501 0.014628074131906;0.965266108512878 0.00746965501457453 0.0175536889582872;0.959477126598358 0.00871459767222404 0.0204793028533459;0.953688144683838 0.00995953939855099 0.0234049186110497;0.947899162769318 0.0112044820562005 0.0263305325061083;0.942110180854797 0.01244942471385 0.0292561482638121;0.936321198940277 0.0136943673714995 0.0321817621588707;0.930532217025757 0.0149393100291491 0.0351073779165745;0.924743235111237 0.016184251755476 0.038032989948988;0.918954253196716 0.0174291953444481 0.0409586057066917;0.913165271282196 0.018674137070775 0.0438842214643955;0.907376289367676 0.019919078797102 0.0468098372220993;0.901587307453156 0.0211640223860741 0.0497354492545128;0.895798325538635 0.022408964112401 0.0526610650122166;0.890009343624115 0.0236539077013731 0.0555866807699203;0.884220361709595 0.0248988494277 0.0585122965276241;0.878431379795074 0.026143791154027 0.0614379085600376;0.872642397880554 0.0273887347429991 0.0643635243177414;0.866853415966034 0.028633676469326 0.0672891363501549;0.861064434051514 0.0298786200582981 0.070214755833149;0.855275452136993 0.0311235617846251 0.0731403678655624;0.849486470222473 0.032368503510952 0.0760659798979759;0.843697488307953 0.0336134470999241 0.07899159938097;0.837908506393433 0.0348583906888962 0.0819172114133835;0.832119524478912 0.036103330552578 0.0848428308963776;0.826330542564392 0.0373482741415501 0.087768442928791;0.820541560649872 0.0385932177305222 0.0906940549612045;0.814752578735352 0.0398381575942039 0.0936196744441986;0.808963596820831 0.041083101183176 0.0965452864766121;0.803174614906311 0.0423280447721481 0.0994708985090256;0.797385632991791 0.0435729846358299 0.10239651799202;0.791596651077271 0.044817928224802 0.105322130024433;0.78580766916275 0.0460628718137741 0.108247749507427;0.78001868724823 0.0473078154027462 0.111173361539841;0.77422970533371 0.048552755266428 0.114098973572254;0.768440723419189 0.0497976988554001 0.117024593055248;0.762651741504669 0.0510426424443722 0.119950205087662;0.756862759590149 0.052287582308054 0.122875817120075;0.751073777675629 0.0535325258970261 0.125801429152489;0.745284795761108 0.0547774694859982 0.128727048635483;0.739495813846588 0.0560224093496799 0.131652668118477;0.733706831932068 0.057267352938652 0.13457827270031;0.727917850017548 0.0585122965276241 0.137503892183304;0.722128868103027 0.0597572401165962 0.140429511666298;0.716339886188507 0.061002179980278 0.143355116248131;0.710550904273987 0.0622471235692501 0.146280735731125;0.704761922359467 0.0634920671582222 0.149206355214119;0.698972940444946 0.064737007021904 0.152131959795952;0.693183958530426 0.0659819543361664 0.155057579278946;0.687394976615906 0.0672268941998482 0.15798319876194;0.681605994701385 0.06847183406353 0.160908818244934;0.675817012786865 0.0697167813777924 0.163834422826767;0.670028030872345 0.0709617212414742 0.166760042309761;0.664239048957825 0.0722066611051559 0.169685661792755;0.658450067043304 0.0734516084194183 0.172611266374588;0.652661085128784 0.0746965482831001 0.175536885857582;0.646872103214264 0.0759414881467819 0.178462505340576;0.641083121299744 0.0771864354610443 0.181388109922409;0.635294139385223 0.0784313753247261 0.184313729405403]);
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

function Rth=throt(phi,theta)
Rz = zrot(-theta);
Rx = xrot(phi);
Rth = inv(Rz)*Rx*Rz;
end