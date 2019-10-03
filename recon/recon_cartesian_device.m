%%-------------------------------%%
%%----MRI reconstruction code----%%
%%-------------------------------%%
% function [img_s,kspace,header] = recon_cartesian_device(file, nfile, user_opts)
% Fully-sampled/inline-GRAPPA Cartesian recon for real-time active device
% channels.
%
% input:    file = ISMRMRD .h5 data file
%           nfile = ISMRMRD .h5 noise dependency file | (optional)
%           user_opts = user_opts structure | [] | (optional)
% default if no argument:
%     user_opts.L7_only = 0;
%     user_opts.SNR_flag = 0;
%     user_opts.Subset_recon = 0;
% Empty input will call a GUI for options
%
% output:   img_s = [structure]
%           img  = reconstructed image
%           header = ISMRMRD .xml header + modified output
%           snr = pseudo-reconstructed SNR
%           timestamp = reconstruction start time
%
% R Ramasawmy & K Yildrim Mar 2019 NHLBI
%
% July overhaul to chunk data for large real-time datasets
% Splits data in to 10 pieces, with some guesswork and assumptions to
% figure out where the repetitions lie.. not sure how it deals with pausing
% etc! Will probably be discontinuous at the chunks for the moment..

%  === test recon ===
% dirPath = '\\hl-share.nhlbi.nih.gov\tmb\Lab-Lederman\Non-Clinical\Yildirim\Low Field Active Guidewire\MRI Scans\07.16.2019 raw data\h5\';
% noisePath = '\\hl-share.nhlbi.nih.gov\tmb\Lab-Lederman\Non-Clinical\Yildirim\Low Field Active Guidewire\MRI Scans\07.16.2019 raw data\noise\';
% test_recon = recon_cartesian_device([dirPath 'meas_MID00051_FID27645_BEAT_interactive_1123_gadgetron_[IA].h5'],[noisePath 'noise_meas_MID00051_FID27645_BEAT_interactive_1123_gadgetron_[IA].h5'],[]);
%  ======

function [img_s] = recon_cartesian_device(file, nfile, user_opts)
%% UI select files and recon options
if nargin < 1
    [fname, dirPath] = uigetfile('*.*', 'Choose data .h5 file');
    [nfname, noisePath] = uigetfile('*.*', 'Choose corresponding noise .h5 file');
    file = [dirPath fname];
    nfile = [noisePath nfname]; clear fname dirPath nfname noisePath;
    SNR_flag = 0;
    user_opts = [];
end

if nargin < 3
    %     user_opts default
    user_opts.L7_only = 0;
    user_opts.SNR_flag = 0;
    user_opts.Subset_recon = 0;
end

if (isempty(user_opts))
    list = {'Recon L7 only','calculate SNR','Recon subset frames'};
    [indx,tf] = listdlg('ListString',list);
    if tf == 0
        user_opts.L7_only = 0;
        user_opts.SNR_flag = 0;
        user_opts.Subset_recon = 0;
    else
        user_opts.L7_only = sum(ismember(indx,1));
        user_opts.SNR_flag = sum(ismember(indx,2));
        user_opts.Subset_recon = sum(ismember(indx,3));
    end
else
    % needs user_opts overwrite option here
end

%% Read data file
df = dir(file);
file_size = round(df.bytes/1e6); % MB
[userview, system] = memory;
sysmem = round(system.PhysicalMemory.Available/1e6); %MB
% sysmem = round(system.SystemMemory.Total/1e6); %MB

if file_size/10 > sysmem
    error(['File-size ' num2str(file_size) 'MB is going to push available memory: ' num2str(sysmem) 'MB.' ]);
    
end
if file_size*3/10 > sysmem
    warning(['File-size ' num2str(file_size) 'MB is going to push available memory: ' num2str(sysmem) 'MB.' ]);
    %error?
    disp('Segmented reconstruction activated');
    % can increase chunk steps here
end
clear userview system file_size df;

% % chunk the data loading
file_info=h5info(file);
data_length = file_info.Groups.Datasets(1).Dataspace.Size; clear file_info;

ismrmrd_s = read_h5_header(file); disp(' ');disp('### Protocol Name ###');disp(ismrmrd_s.measurementInformation.protocolName);disp(' ');
img_s.timestamp = datetime(now, 'ConvertFrom', 'datenum');

%% == scout the number of repetitions ==
scout_steps = 5 * ismrmrd_s.encoding.encodingLimits.kspace_encoding_step_1.maximum * ismrmrd_s.encoding.encodingLimits.slice.maximum;
if scout_steps > data_length; scout_steps = data_length; end;
raw_data_scout = h5read(file, '/dataset/data', 1, scout_steps);

rep_step = mean(diff(find(diff(single(raw_data_scout.head.idx.repetition))))) + 1;
slice_step = mean(diff(find(diff(single(raw_data_scout.head.idx.slice))))) + 1;
steps_per_scan = min([rep_step slice_step]); clear rep_step slice_step
if rem(data_length,steps_per_scan) > 0
    disp('Scan stopped before end of last acquisition, not reconstructing this one');
end
num_scans = floor(data_length/steps_per_scan);

%% Extract header info

samples = double(raw_data_scout.head.number_of_samples(1));
% asymmetric echo?
asym_e = 0; if (samples~=ismrmrd_s.encoding.encodedSpace.matrixSize.x); asym_e=1; end;
echo_vec = (ismrmrd_s.encoding.encodedSpace.matrixSize.x-samples+1):ismrmrd_s.encoding.encodedSpace.matrixSize.x;
channels = double(raw_data_scout.head.active_channels(1));
dt = raw_data_scout.head.sample_time_us(1)*1e-6;

pe1 = ismrmrd_s.encoding.reconSpace.matrixSize.y;
pe2 = ismrmrd_s.encoding.reconSpace.matrixSize.z;

% defaults to max set reps % reps = ismrmrd_s.encoding.encodingLimits.repetition.maximum; % reps = double(max(raw_data.head.idx.repetition)+1);
raw_data_scout = h5read(file, '/dataset/data', steps_per_scan*num_scans, 1);
reps        = double(max(raw_data_scout.head.idx.repetition)+1);
slices      = 1 + ismrmrd_s.encoding.encodingLimits.slice.maximum; % slices = double(max(raw_data.head.idx.slice))+1;
averages    = 1 + ismrmrd_s.encoding.encodingLimits.average.maximum; % averages = double(max(raw_data_scout.head.idx.average))+1;
contrasts   = 1 + ismrmrd_s.encoding.encodingLimits.contrast.maximum; % contrasts = double(max(raw_data.head.idx.contrast))+1;
phases      = 1 + ismrmrd_s.encoding.encodingLimits.phase.maximum; % phases = double(max(raw_data.head.idx.phase))+1;
sets        = 1 + ismrmrd_s.encoding.encodingLimits.set.maximum; % sets = double(max(raw_data.head.idx.set))+1;
clear raw_data_scout;

% im_header.samples = samples;
% im_header.dt = dt;
% im_header.number_aqs = num_scans;
% im_header.averages = averages;
% im_header.channels = channels;
%
% header.im_header = im_header;

%% Print out to window

if (samples < ismrmrd_s.encoding.encodedSpace.matrixSize.x); disp('Asymmetric Echo'); end;
% disp(['BW: ' num2str(dt)])
acc_factor = [ismrmrd_s.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 ismrmrd_s.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2];
disp(' ');disp('### Acquisition Dimensions ###');disp(' ');
header_info = {'Encoded_Res','Encoded_FOV','Recon_Res','Recon_FOV'}';
X_dim = [ismrmrd_s.encoding.encodedSpace.matrixSize.x ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.x ismrmrd_s.encoding.reconSpace.matrixSize.x ismrmrd_s.encoding.reconSpace.fieldOfView_mm.x]';
Y_dim = [ismrmrd_s.encoding.encodedSpace.matrixSize.y ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.y ismrmrd_s.encoding.reconSpace.matrixSize.y ismrmrd_s.encoding.reconSpace.fieldOfView_mm.y]';
Z_dim = [ismrmrd_s.encoding.encodedSpace.matrixSize.z ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.z ismrmrd_s.encoding.reconSpace.matrixSize.z ismrmrd_s.encoding.reconSpace.fieldOfView_mm.z]';
disp(table(header_info, X_dim, Y_dim, Z_dim)); clear header_info X_dim Y_dim Z_dim;

disp(' ');disp('### Experiment Dimensions ###');disp(' ');
Experiment_parameters = {'RO', 'PE1', 'PE2', 'Averages', 'Slices', 'Contrasts', 'Phases', 'Reps', 'Sets', 'Channels'}';
Value = [ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 pe2 averages slices contrasts phases reps sets channels]';
disp(table( Experiment_parameters,Value )); clear Experiment_parameters Value; disp(' ');

%% L7 detection and noise checks

disp(['Dependency ID: ' ismrmrd_s.measurementInformation.measurementDependency.measurementID]);
asi_names = fieldnames(ismrmrd_s.acquisitionSystemInformation);
n_iRD = read_h5_header(nfile);
nasi_names = fieldnames(n_iRD.acquisitionSystemInformation);

solid_int = 0;
coil_label = cell(channels,4);
coil_flag = zeros(channels,2);
for i = 1:length(asi_names)
    if regexp(asi_names{i}, 'coilLabel')
        %         [asi_names{i} ' ' nasi_names{i}]
        solid_int = solid_int + 1;
        coil_label{solid_int,1} = ismrmrd_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilNumber;
        coil_label{solid_int,2} = ismrmrd_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilName;
        
        if regexp(ismrmrd_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilName, 'L7');
            L7_channel = solid_int;
        end
        
        coil_label{solid_int,3} = n_iRD.acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilNumber;
        coil_label{solid_int,4} = n_iRD.acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilName;
        
        %         coil_label{solid_int,5} = coil_label{solid_int,1} == coil_label{solid_int,3};
        %         coil_label{solid_int,6} = regexp(coil_label{solid_int,2}, coil_label{solid_int,4});
        coil_flag(solid_int,1) = coil_label{solid_int,1} == coil_label{solid_int,3};
        coil_flag(solid_int,2) = regexp(coil_label{solid_int,2}, coil_label{solid_int,4});
    end
end

image_channels = 1:channels;
if(exist('L7_channel'))
    image_channels(L7_channel)= []; % requires a error catch for non-L7 scans. KOREL.
else
    warning('No L7 channel, reconstructing image channels')
    
    % switch override:
    user_opts.L7_only = 1;
    L7_channel = image_channels;
end

if sum(coil_flag(:,1)) ~= channels
    disp('Coil order mismatch!');
elseif sum(coil_flag(:,2)) ~= channels
    disp('Coil name mismatch!');
end

%% Calculate Noise Decorrelation Matrix
% Need to separate L7 and imaging coil noise-decorr
% or just dont do it..

if ~exist(nfile) || isempty(nfile)
    disp('No noise adjustment: ');
    dmtx_imaging = diag(ones(1,channels-1));
    dmtx_L7 = 1;
else
    [dmtx_imaging,dmtx_L7]= ismrm_dmtx_device(nfile, dt, L7_channel); % inline function
end

% if nargin > 1
%     if isempty(nfile)
%         disp('No noise adjustment: ');
%         dmtx = diag(ones(1,channels));
%     else
%         disp('Required noise ID: ');
%         disp(ismrmrd_s.measurementInformation.measurementDependency.measurementID);
% dmtx = ismrm_dmtx_RR(nfile, dt);
%     end
% else

% end
disp(' ');

%% User selection/ chunking
if user_opts.Subset_recon == 1
    % User selection
    
    prompt = {'Repetitions start:','Repetitions end:'};
    dlgtitle = 'Input subset recon range (by number of scans)';
    dims = [1 35];
    definput = {'1',num2str(num_scans)};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    range_start = str2double(answer{1});
    range_finish = str2double(answer{2});
    if range_finish > num_scans; range_finish = num_scans; end;
    
    clear prompt dlgtitle dims definput answer
    
    if range_finish-range_start > 25
        chunk_step = round(linspace(range_start, range_finish, 11));
    else
        chunk_step = round(linspace(range_start, range_finish, 2));
    end
    
    chunk_step = chunk_step*steps_per_scan;
    reps = range_finish - range_start;
    
else
    % partition large datasets (>25 reps) in to ten chunks
    
    if num_scans > 25
        chunk_step = round(linspace(0, num_scans, 11));
        chunk_step(2:end) = chunk_step(2:end)*steps_per_scan;
    else
        chunk_step = round(linspace(1, num_scans, 2));
    end
    
end

%% Chunked data loading & recon
% load all the data
% raw_data= h5read(file, '/dataset/data');

% initialise FFT images (for real-time intereactive, i.e. only slice and rep extra dims)
if user_opts.L7_only == 0
    %     cCoil_imgs      = zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 pe2 averages slices contrasts phases reps sets ]);
    cCoil_imgs      = zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 slices reps]);
end

L7_imgs         = zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 slices reps]);

% initialise grappa images
if acc_factor(1) > 1
    if user_opts.L7_only == 0
        cCoil_imgs_g         = zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 slices reps]);
    end
    
    L7_imgs_g         = zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 slices reps]);
end

for iCS = 1:length(chunk_step)-1
    
    range_start  = chunk_step(iCS)+1;
    range_finish = chunk_step(iCS+1);
    
    raw_data = h5read(file, '/dataset/data', range_start, range_finish-range_start+1);
    
    tempreps = 1 + double(max(raw_data.head.idx.repetition)+1) - double(min(raw_data.head.idx.repetition)+1);
    range_start = double(min(raw_data.head.idx.repetition)+1);
    
    kspace = complex(zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 slices tempreps channels],'single'));
    %     whos kspace
    for ii = 1:length(raw_data.data);
        
        d1 = raw_data.data{ii};
        d2 = complex(d1(1:2:end), d1(2:2:end));
        d3 = reshape(d2, samples, channels); %  RE & IM (2)
        
        %         d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
        d3(:,image_channels) = ismrm_apply_noise_decorrelation_mtx(d3(:,image_channels),dmtx_imaging);
        
        if length(L7_channel)>1
            d3(:,L7_channel) = ismrm_apply_noise_decorrelation_mtx(d3(:,L7_channel),dmtx_L7);
            
        else
            d3(:,L7_channel) = d3(:,L7_channel)*dmtx_L7;
        end
        
        
        kspace(uint16(echo_vec),...
            raw_data.head.idx.kspace_encode_step_1(ii)+1 + uint16(ismrmrd_s.encoding.encodedSpace.matrixSize.y/2) - uint16(ismrmrd_s.encoding.encodingLimits.kspace_encoding_step_1.center), ...
            raw_data.head.idx.slice(ii)+1 , ...
            raw_data.head.idx.repetition(ii)+1 - (range_start - 1), ...
            :) = d3;
        
    end
    
    %     whos kspace
    %    size(kspace)
    
    %% fft recon
    
    % check chunk indexing :
    %     [min(raw_data.head.idx.repetition+1) max(raw_data.head.idx.repetition+1)]
    %     range_start:(1+max(raw_data.head.idx.repetition))
    
    if user_opts.L7_only == 0
        temp = ismrm_transform_kspace_to_image(kspace(:,:,:,:,image_channels), [1 2]);
        cCoil_imgs(:,:,:,range_start:(1+max(raw_data.head.idx.repetition))) = ismrm_rss(temp);
    end
    
    temp = ismrm_transform_kspace_to_image(kspace(:,:,:,:,L7_channel), [1 2]);
    if length(L7_channel) > 1 % imaging channel override (i.e. no L7 channel present)
        temp = ismrm_rss(temp);
    end
    
    L7_imgs(:,:,:,range_start:(1+max(raw_data.head.idx.repetition))) = temp;
    
    if acc_factor(1) > 1
        
        %% grappa recon
        
        % should move outside the chunk loop >>
        
        sampling_pattern = squeeze(kspace(:,:,raw_data.head.idx.slice(round(ii/2))+1,raw_data.head.idx.repetition(round(ii/2))+1 - (range_start - 1),1));
        sampling_pattern = abs(sampling_pattern) > 0;
        
        sampling_pattern_1 = sampling_pattern(1,:);
        %     acc = round(ismrmrd_s.encoding.encodingLimits.kspace_encoding_step_1.maximum/sum(sampling_pattern_1));
        acc = acc_factor(1);
        
        temp0 = find(sampling_pattern_1);
        temp1 = diff(temp0);
        temp2 = temp0(find(temp1==1));
        sampling_pattern = single(sampling_pattern);
        
        sampling_pattern(:,temp2) = sampling_pattern(:,temp2)*3; disp('Assuming in-line GRAPPA'); % ASSUMING IN-LINE REFERENCE
        
        % should move outside the chunk loop <<
        
        % figure, imshow(sampling_pattern)
        % figure, plot(sampling_pattern_1); hold on, plot(temp2, sampling_pattern_1(temp2), 'r-'); ylim([-0.1 1.1]);
        
        for slc = 1:slices
            for repc = 1:tempreps
                if user_opts.L7_only == 0
                    cCoil_imgs_g(:,:,slc,repc - 1 + range_start) = abs(ismrm_cartesian_GRAPPA(squeeze(kspace(:,:,slc,repc,image_channels)),sampling_pattern, acc));
                end
                
                temp = ismrm_transform_kspace_to_image(kspace(:,:,slc,repc,L7_channel));
                csm_L7 = ismrm_estimate_csm_mckenzie(squeeze(temp));
                L7_imgs_g(:,:,slc,repc - 1 + range_start) = abs(ismrm_cartesian_GRAPPA(squeeze(kspace(:,:,slc,repc,L7_channel)),sampling_pattern, acc, csm_L7));
                
            end
        end
        
    end
end

%% remove OS

OverSampling = ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.x./ismrmrd_s.encoding.reconSpace.fieldOfView_mm.x;
if OverSampling == 2
    if user_opts.L7_only == 0
        cCoil_imgs = remove_OS(cCoil_imgs);
        if acc_factor(1) > 1
            cCoil_imgs_g = remove_OS(cCoil_imgs_g);
        end
    end
    
    L7_imgs = remove_OS(L7_imgs);
    if acc_factor(1) > 1
        L7_imgs_g = remove_OS(L7_imgs_g);
    end
end


%% ==== VISUALISE DATA ====
% Needs some intelligence to display subset reps only

if user_opts.L7_only == 0
    temp_ci = cCoil_imgs;
end
temp_li = L7_imgs;

if acc_factor(1) > 1;
    if user_opts.L7_only == 0
        temp_ci = cCoil_imgs_g;
    end
    temp_li = L7_imgs_g;
end

if user_opts.L7_only == 0
    figure
    for i = 1:reps
        % Green overlay
        E = abs(squeeze(temp_ci(:,:,:,i)));
        E = E(:,:)./max(E(:));
        I = abs(squeeze(temp_li(:,:,:,i)));
        I = I(:,:)./max(I(:));
        
        % show image
        imshow(E, 'InitialMag', 'fit');
        title(['Rep: ' num2str(i)]);
        % Make a truecolor all-green image.
        green = cat(3, zeros(size(E)),ones(size(E)), zeros(size(E)));
        hold on
        h = imshow(green);
        hold off
        
        set(h, 'AlphaData', I) ;
        
        pause(0.1);
        
    end
else
    implay_RR(reshape(temp_li,[size(temp_li,1) size(temp_li,2)*size(temp_li,3) size(temp_li,4) ]))
end

%% ==== Return data structure img_s ====

if user_opts.L7_only == 0
    img_s.cCoil_imgs    = squeeze(cCoil_imgs);
    
    if acc_factor(1) > 1
        img_s.cCoil_imgs_g  = squeeze(cCoil_imgs_g);
    end
end

if acc_factor(1) > 1
    img_s.L7_imgs_g     = squeeze(L7_imgs_g);
end

img_s.L7_imgs       = squeeze(L7_imgs);
img_s.header        = ismrmrd_s;

end

%% ==== inline functions  ====
function new_imgs = remove_OS(imgs)
%   xdim = ismrmrd_s.encoding.encodedSpace.matrixSize.x;
xdim = size(imgs, 1);
temp = reshape(imgs, [xdim prod(size(imgs))/xdim]);
temp = temp( (xdim/4): (xdim/4 + xdim/2 -1), : );
new_dims = size(imgs); new_dims(1) = xdim/2; %ismrmrd_s.encoding.reconSpace.matrixSize.x;
new_imgs = reshape(temp, new_dims);
end


function [dmtx_imaging, dmtx_L7] = ismrm_dmtx_device(nfile, data_samp_time, L7_channel)
% Process noise data and return decorrelation matrix
%% Load and sort channels
noise_test = h5read(nfile,'/dataset/data');
iRD_s = read_h5_header(nfile);
disp('Noise ID: ');
disp(iRD_s.measurementInformation.measurementID);
Siemens_rBW = iRD_s.acquisitionSystemInformation.relativeReceiverNoiseBandwidth;
% Siemens_rBW = 0.79;

n_samples = double(noise_test.head.number_of_samples(1));
n_channels = double(noise_test.head.active_channels(1));
% assuming Siemens noise cal, using 2 averages:
noise_ind = find(noise_test.head.idx.average==1, 1, 'last');

image_channels = 1:n_channels;
image_channels(L7_channel)= []; % requires a error catch for non-L7 scans. KOREL.


%% Calculate DMTX
nt2 = zeros(n_samples, noise_ind, n_channels);
for i = 1:noise_ind
    nt2(:,i,:)=  double(reshape(complex(noise_test.data{i}(1:2:end), noise_test.data{i}(2:2:end)), [n_samples, 1, n_channels ]));
end
n_scaling = Siemens_rBW * data_samp_time / (noise_test.head.sample_time_us(1)*1e-6);

if isempty(image_channels)
    dmtx_imaging = [];
else
    dmtx_imaging = ismrm_calculate_noise_decorrelation_mtx(nt2(:,:,image_channels), n_scaling ); % figure,imagesc(abs(dmtx));
end
if length(L7_channel) > 1 % image recon override
    dmtx_L7 = ismrm_calculate_noise_decorrelation_mtx(nt2(:,:,L7_channel), n_scaling ); % figure,imagesc(abs(dmtx));
else
    % single channel decorrelation factor
    noise = nt2(:,:,L7_channel);
    noise = reshape(noise, 1, prod(size(noise)));
    dmtx_L7 = (1./(length(noise)-1))*(noise*noise');
    dmtx_L7 = 1./(sqrt(dmtx_L7));
    dmtx_L7 = dmtx_L7*sqrt(2)*sqrt(n_scaling);
    
end

end

function [out] = ismrm_apply_noise_decorrelation_mtx_device(in, dmtx)

end

function snr = cartesian_pseudoreps(crt_k)
% SNR_s = cartesian_psuedoreps(crt_k, slice)
% Expects pre-whitened 2D k-space data [RO,PE,AV,COILS]
% This can be expanded to generic dimensions, but will require header info
% on dimension - better used as a subfunction within recon_cartesian
% R Ramasawmy NHLBI 2019

if length(size(crt_k)) == 3
    crt_k = reshape(crt_k,[size(crt_k,1) size(crt_k,2) 1 size(crt_k,3)]);
end
% size(crt_k)
av_dim = 3;

pseudo_reps = 100;
disp(['Running ' num2str(pseudo_reps) ' pseudo-reps']);

% estimate csm from data
ncoils = size(crt_k,ndims(crt_k));
img = ismrm_transform_kspace_to_image(squeeze(mean(crt_k,av_dim)),[1 2]); % montage_RR(abs(img));figure, imshow(sqrt(sum(img.*conj(img),3)),[])
csm = ismrm_estimate_csm_walsh(img); % montage_RR((csm));
% calculate coil combine
ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(ncoils)); % with pre-whitened

img_pr = zeros(size(crt_k,1), size(crt_k,2), pseudo_reps);
for i = 1:pseudo_reps
    data_pr = crt_k + complex(randn(size(crt_k)),randn(size(crt_k)));
    data_pr = squeeze(mean(data_pr,av_dim));
    x = ismrm_transform_kspace_to_image(data_pr,[1 2]);
    %     img_pr(:,:,i) = sqrt(sum(x.*conj(x),3)); % NO!
    img_pr(:,:,i) = abs(sum(x .* ccm_roemer_optimal, 3));
end

% img_pr_orig = img_pr;
img_pr = img_pr( size(crt_k,1)/4:size(crt_k,1)/4 +size(crt_k,1)/2 -1, :,:);

g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
g(g < eps) = 1;
snr = mean(img_pr,3)./g;
%
% figure,
% subplot(1,2,1); imshow(snr,[0 50]); colorbar; colormap('parula')
% subplot(1,2,2); imshow(std(img_pr,[],3),[]); colorbar; colormap('parula')
%
% SNR_s.img_pr = img_pr;
% SNR_s.snr = snr;
% SNR_s.gr = g;

end