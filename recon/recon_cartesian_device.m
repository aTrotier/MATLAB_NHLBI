%%-------------------------------%%
%%----MRI reconstruction code----%%
%%-------------------------------%%
% function [img_s,kspace,header] = recon_cartesian(file, <nfile>, <SNR_flag>)
% Fully-sampled Cartesian recon (handles partial-fourier asymmetric-echo
% and interpolation).
% input:    file = ISMRMRD .h5 data file 
%           nfile = ISMRMRD .h5 noise dependency file (optional)
% output:   img_s = [structure] 
    %           img  = reconstructed image
    %           header = ISMRMRD .xml header + modified output
    %           snr = pseudo-reconstructed SNR
%           kspace = raw k-space (prewhitened if noise-file is present)
%
% R Ramasawmy Dec 2018 NHLBI 

function [img_s,kspace] = recon_cartesian_device(file, nfile, SNR_flag)
%%
if nargin < 1
     [fname, dirPath] = uigetfile('*.*');
     [nfname, noisePath] = uigetfile('*.*');
     file = [dirPath fname];
     nfile = [noisePath nfname]; clear fname dirPath nfname noisePath;
     SNR_flag = 0;
end

%% Read data file
raw_data= h5read(file, '/dataset/data');
ismrmrd_s = read_h5_header(file); disp(' ');disp('### Protocol Name ###');disp(ismrmrd_s.measurementInformation.protocolName);disp(' ');
header = ismrmrd_s;

samples = double(raw_data.head.number_of_samples(1));
% asymmetric echo?
asym_e = 0; if (samples~=ismrmrd_s.encoding.encodedSpace.matrixSize.x); asym_e=1; end;
echo_vec = (ismrmrd_s.encoding.encodedSpace.matrixSize.x-samples+1):ismrmrd_s.encoding.encodedSpace.matrixSize.x;
channels = double(raw_data.head.active_channels(1));

pe1 = ismrmrd_s.encoding.reconSpace.matrixSize.y; % pe1 = 1+double(max(raw_data.head.idx.kspace_encode_step_1));
pe2 = ismrmrd_s.encoding.reconSpace.matrixSize.z; % pe2 = 1+double(max(raw_data.head.idx.kspace_encode_step_2));

averages = double(max(raw_data.head.idx.average))+1;
slices = double(max(raw_data.head.idx.slice))+1;
contrasts = double(max(raw_data.head.idx.contrast))+1;
phases = double(max(raw_data.head.idx.phase))+1;
reps = double(max(raw_data.head.idx.repetition)+1);
sets = double(max(raw_data.head.idx.set))+1;
dt = raw_data.head.sample_time_us(1)*1e-6;

im_header.samples = samples;
im_header.dt = dt;
im_header.number_aqs = length(raw_data.data);
im_header.averages = averages;
im_header.channels = channels;

header.im_header = im_header;

%% Print out to window 

if (samples < ismrmrd_s.encoding.encodedSpace.matrixSize.x); disp('Asymmetric Echo'); end;
% disp(['BW: ' num2str(dt)])

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

%% Noise checks
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
image_channels(L7_channel)= [];

if sum(coil_flag(:,1)) ~= channels
    disp('Coil order mismatch!');
elseif sum(coil_flag(:,2)) ~= channels
    disp('Coil name mismatch!');
end




%% Calculate Noise Decorrelation Matrix
if nargin > 1
    if isempty(nfile)
        disp('No noise adjustment: ');
        dmtx = diag(ones(1,channels));
    else
        disp('Required noise ID: ');
        disp(ismrmrd_s.measurementInformation.measurementDependency.measurementID);
        dmtx = ismrm_dmtx_RR(nfile, dt);
    end
else
    disp('No noise adjustment: ');
    dmtx = diag(ones(1,channels));
end
disp(' ');




%%
if(length(unique(double(raw_data.head.number_of_samples))) > 1)
    % non-uniform sample number (grappa reference lines?)
    error('GRAPPA recon needed')
    
else

    kspace = zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 pe2 averages slices contrasts phases reps sets channels]);
    % disp(['Kspace dims: ' num2str(size(kspace))])

    for ii = 1:length(raw_data.data)
                
        d1 = raw_data.data{ii};
        d2 = complex(d1(1:2:end), d1(2:2:end));
        d3 = reshape(d2, samples, channels); %  RE & IM (2)
        
        if nargin > 1
            d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
        end
        
        kspace(echo_vec,...
            raw_data.head.idx.kspace_encode_step_1(ii)+1 + ismrmrd_s.encoding.encodedSpace.matrixSize.y/2 - ismrmrd_s.encoding.encodingLimits.kspace_encoding_step_1.center, ...
            raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
            raw_data.head.idx.average(ii)+1, ....
            raw_data.head.idx.slice(ii)+1 , ...
            raw_data.head.idx.contrast(ii)+1, ...
            raw_data.head.idx.phase(ii)+1, ...
            raw_data.head.idx.repetition(ii)+1, ...
            raw_data.head.idx.set(ii)+1, ...
            :) = d3;
     
    end
end

%% Image Recon
% detect grappa

sampling_pattern = squeeze(kspace(:,:,1,1,1,1,1,1,1,1));
sampling_pattern = abs(sampling_pattern) > 0; 

sampling_pattern_1 = sampling_pattern(1,:);
acc = round(ismrmrd_s.encoding.encodingLimits.kspace_encoding_step_1.maximum/sum(sampling_pattern_1));

if acc > 1
    temp0 = find(sampling_pattern_1);
    temp1 = diff(temp0);
    temp2 = temp0(find(temp1==1));
    sampling_pattern = single(sampling_pattern);
    
    sampling_pattern(:,temp2) = sampling_pattern(:,temp2)*3; disp('Assuming in-line GRAPPA'); % ASSUMING IN-LINE REFERENCE
end

% figure, plot(sampling_pattern_1); hold on, plot(temp2, sampling_pattern_1(temp2), 'r-'); ylim([-0.1 1.1]);

% figure, 
% subplot(1,2,1); imshow(abs(image_grappa),[0 200])
% subplot(1,2,2); imshow(abs(L7_alias)./max(abs(L7_alias(:))),[0 1])

%% transform to image-space

kspace = mean(kspace,4);% kspace = squeeze(kspace);cCoil_imgs = ismrm_transform_kspace_to_image(kspace);

image_dims = size(kspace);
image_dims = image_dims(1:9);

% Just FFT images
cCoil_imgs      = zeros(image_dims);
L7_imgs         = zeros(image_dims);

% Grappa reconstructed images
cCoil_imgs_g    = zeros(image_dims);
L7_imgs_g       = zeros(image_dims);

tic;
for slc = 1:slices
    for coc = 1:contrasts
        for phc = 1:phases
            for repc = 1:reps
                for setc = 1:sets
                    % FFT & GRAPPA recon
                    
                    cCoil_imgs(:,:,:,1,slc,coc,phc,repc,setc) = ismrm_rss(ismrm_transform_kspace_to_image(kspace(:,:,:,1,slc,coc,phc,repc,setc,image_channels)),ndims(kspace));
                    
                    cCoil_imgs_g(:,:,:,1,slc,coc,phc,repc,setc) = ismrm_cartesian_GRAPPA(squeeze(kspace(:,:,1,1,slc,coc,phc,repc,setc,image_channels)),sampling_pattern, acc);
                    
                    temp = ismrm_transform_kspace_to_image(kspace(:,:,:,1,slc,coc,phc,repc,setc,L7_channel));
                    L7_imgs(:,:,:,1,slc,coc,phc,repc,setc) = ismrm_rss(temp, ndims(kspace));
                    csm_L7 = ismrm_estimate_csm_mckenzie(temp);
                    L7_imgs_g(:,:,:,1,slc,coc,phc,repc,setc) = ismrm_cartesian_GRAPPA(squeeze(kspace(:,:,1,1,slc,coc,phc,repc,setc,L7_channel)),sampling_pattern, acc, csm_L7);
                    
                end
            end
        end
    end
end
toc;

figure, 
subplot(2,2,1); imshow(abs(cCoil_imgs(:,:,:,1,slc,coc,phc,repc,setc)),[0 200])
subplot(2,2,2); imshow(abs(L7_imgs(:,:,:,1,slc,coc,phc,repc,setc)),[0 50])
subplot(2,2,3); imshow(abs(cCoil_imgs_g(:,:,:,1,slc,coc,phc,repc,setc)),[0 200])
subplot(2,2,4); imshow(abs(L7_imgs_g(:,:,:,1,slc,coc,phc,repc,setc)),[0 50])

%% remove OS
OverSampling = ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.x./ismrmrd_s.encoding.reconSpace.fieldOfView_mm.x;
if OverSampling == 2
    cCoil_imgs = remove_OS(cCoil_imgs);
    cCoil_imgs_g = remove_OS(cCoil_imgs_g);
    L7_imgs = remove_OS(L7_imgs);
    L7_imgs_g = remove_OS(L7_imgs_g);
end

%% PLAY DATA

% implay_RR([squeeze(abs(cCoil_imgs(:,:,:,1,1,1,1,:,1))) squeeze(abs(cCoil_imgs_g(:,:,:,1,1,1,1,:,1)))])
% implay_RR([squeeze(abs(L7_imgs(:,:,:,1,1,1,1,:,1))) squeeze(abs(L7_imgs_g(:,:,:,1,1,1,1,:,1)))])

% non -grappa recon
% figure,
% for i = 1:reps
%     E = abs(squeeze(cCoil_imgs(:,:,:,1,:,coc,phc,i,setc)));
%     E = E(:,:)./max(E(:));
%     I = abs(squeeze(L7_imgs(:,:,:,1,:,coc,phc,i,setc)));
%     I = I(:,:)./max(I(:));
%     
%     
%     imshow(E, 'InitialMag', 'fit')
%     % Make a truecolor all-green image.
%     green = cat(3, zeros(size(E)),ones(size(E)), zeros(size(E)));
%     hold on
%     h = imshow(green);
%     hold off
%     
%     set(h, 'AlphaData', I) ;
%         
%     pause(0.1);
% end

figure,
for i = 1:reps
    E = abs(squeeze(cCoil_imgs_g(:,:,:,1,:,coc,phc,i,setc)));
    E = E(:,:)./max(E(:));
    I = abs(squeeze(L7_imgs_g(:,:,:,1,:,coc,phc,i,setc)));
    I = I(:,:)./max(I(:));
    
    
    imshow(E, 'InitialMag', 'fit')
    % Make a truecolor all-green image.
    green = cat(3, zeros(size(E)),ones(size(E)), zeros(size(E)));
    hold on
    h = imshow(green);
    hold off
    
    set(h, 'AlphaData', I) ;
        
    pause(0.1);
end

% overlay movies - needs some Dan H. 
% E = abs(squeeze(cCoil_imgs_g(:,:,:,1,:,coc,phc,:,setc)));
% E = reshape(E, [size(E,1), size(E,2)*slices, size(E,4)]);
% E = E./max(E(:));
% 
% I = abs(squeeze(L7_imgs_g(:,:,:,1,:,coc,phc,:,setc)));
% I = reshape(I, [size(I,1), size(I,2)*slices, size(I,4)]);
% I = I./max(I(:));
% 
% Rough thresholder 
% implay_DAH( E , round(I/(15*mean(I(:)))) );
% 

%% Return variables

img_s.cCoil_imgs    = squeeze(cCoil_imgs);
img_s.cCoil_imgs_g  = squeeze(cCoil_imgs_g);
img_s.L7_imgs       = squeeze(L7_imgs);
img_s.L7_imgs_g     = squeeze(L7_imgs_g);
img_s.header        = header;

end

function new_imgs = remove_OS(imgs)
%   xdim = ismrmrd_s.encoding.encodedSpace.matrixSize.x;
xdim = size(imgs, 1);
    temp = reshape(imgs, [xdim prod(size(imgs))/xdim]);
    temp = temp( (xdim/4): (xdim/4 + xdim/2 -1), : );
    new_dims = size(imgs); new_dims(1) = xdim/2; %ismrmrd_s.encoding.reconSpace.matrixSize.x;
    new_imgs = reshape(temp, new_dims);
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