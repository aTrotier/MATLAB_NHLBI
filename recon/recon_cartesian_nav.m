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

function [img_s,kspace] = recon_cartesian_nav(file, nfile, SNR_flag)
%% Read data file
if nargin == 1
    nfile = [];
end

% file = RR_run_on_mac(file); % incorporate with NHLBI toolbox
% nfile = RR_run_on_mac(nfile);
%%
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

IM_R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ];
 
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

%% Calculate Noise Decorrelation Matrix

if isempty(nfile)
    disp('No noise adjustment: ');
    dmtx = diag(ones(1,channels));
else
    disp('Required noise ID: ');
    disp(ismrmrd_s.measurementInformation.measurementDependency.measurementID);
    
    %% Noise checks
    
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
            
            coil_label{solid_int,3} = n_iRD.acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilNumber;
            coil_label{solid_int,4} = n_iRD.acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilName;
            
            %         coil_label{solid_int,5} = coil_label{solid_int,1} == coil_label{solid_int,3};
            %         coil_label{solid_int,6} = regexp(coil_label{solid_int,2}, coil_label{solid_int,4});
            coil_flag(solid_int,1) = coil_label{solid_int,1} == coil_label{solid_int,3};
            coil_flag(solid_int,2) = regexp(coil_label{solid_int,2}, coil_label{solid_int,4});
        end
    end
    
    if sum(coil_flag(:,1)) ~= channels
        disp('Coil order mismatch!');
    elseif sum(coil_flag(:,2)) ~= channels
        disp('Coil name mismatch!');
    end
    
    dmtx = ismrm_dmtx_RR(nfile, dt);
end

disp(' ');

%% Recon
% navigator extract
% 
nav_samples = single(raw_data.head.number_of_samples(1));
data_samples = single(raw_data.head.number_of_samples(3));

acq_frames = find(single(raw_data.head.number_of_samples) == data_samples);
nav_frames = find(single(raw_data.head.number_of_samples) == nav_samples);
      
    echo_vec = 1:data_samples;
    kspace = complex(zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 pe2 averages slices contrasts phases reps sets channels],'single'));
    % disp(['Kspace dims: ' num2str(size(kspace))])
    
    for ii = acq_frames'%2:length(raw_data.data) % %
        
        d1 = raw_data.data{ii};
        d2 = complex(d1(1:2:end), d1(2:2:end));
        d3 = reshape(d2, data_samples, channels); %  RE & IM (2)
        
        if nargin > 1
            d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
        end
        
        kspace(echo_vec,...
            raw_data.head.idx.kspace_encode_step_1(ii)+1, ...
            raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
            raw_data.head.idx.average(ii)+1, ....
            raw_data.head.idx.slice(ii)+1 , ...
            raw_data.head.idx.contrast(ii)+1, ...
            raw_data.head.idx.phase(ii)+1, ...
            raw_data.head.idx.repetition(ii)+1, ...
            raw_data.head.idx.set(ii)+1, ...
            :) = d3;
        
    end
    
    echo_vec = 1:nav_samples;
    nav_kspace = complex(zeros([nav_samples pe1 pe2 averages slices contrasts phases reps sets channels],'single'));
    
    for ii = nav_frames'%2:length(raw_data.data) % %
        
        d1 = raw_data.data{ii};
        d2 = complex(d1(1:2:end), d1(2:2:end));
        d3 = reshape(d2, nav_samples, channels); %  RE & IM (2)
        
        if nargin > 1
            d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
        end
        
        nav_kspace(:,...
            raw_data.head.idx.kspace_encode_step_1(ii)+1, ...
            raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
            raw_data.head.idx.average(ii)+1, ....
            raw_data.head.idx.slice(ii)+1 , ...
            raw_data.head.idx.contrast(ii)+1, ...
            raw_data.head.idx.phase(ii)+1, ...
            raw_data.head.idx.repetition(ii)+1, ...
            raw_data.head.idx.set(ii)+1, ...
            :) = d3;
        
    end

    nav_kspace_store = complex(zeros([nav_samples length(nav_frames) channels],'single')); %128*64

    for ii = 1:length(nav_frames)%2:length(raw_data.data) % %
        
        d1 = raw_data.data{nav_frames(ii)};
        d2 = complex(d1(1:2:end), d1(2:2:end));
        d3 = reshape(d2, nav_samples, channels); %  RE & IM (2)
        
        if nargin > 1
            d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
        end
        
        nav_kspace_store(:,ii,:) = d3;
        
    end

%% SNR
if nargin < 3
    SNR_flag = 0;
end

if SNR_flag
    
    if isempty(SNR_flag)
        pseudoRep_slices = 1:slices;
    end
    if SNR_flag > 0 % specify slice(s);
        pseudoRep_slices = SNR_flag;
        
        if pseudoRep_slices > slices
            pseudoRep_slices = RR_slice_order(round(slices/2));
            warning(['PR slice > #slices, using slice ' num2str(pseudoRep_slices)])
        end
    end
    
    %     snr = zeros(ismrmrd_s.encoding.reconSpace.matrixSize.x, ismrmrd_s.encoding.reconSpace.matrixSize.y, slices);
    snr = zeros(ismrmrd_s.encoding.reconSpace.matrixSize.x, size(kspace,2), slices);
    
    for iSlice = pseudoRep_slices
        %         snr(:,:,iSlice) = cartesian_pseudoreps(squeeze(kspace(:,:,:,:,iSlice,:,:,:,:,:)));
        snr(:,:,iSlice) = cartesian_pseudoreps(squeeze(kspace(:,:,:,:,iSlice,1,1,1,1,:)));
    end
end


%% transform to image-space

kspace = mean(kspace,4);% kspace = squeeze(kspace);cCoil_imgs = ismrm_transform_kspace_to_image(kspace);

if pe2 == 1
    % 2D
    cCoil_imgs = ismrm_transform_kspace_to_image(kspace, [1 2]);    
else
    % 3D
    cCoil_imgs = ismrm_transform_kspace_to_image(kspace, [1 2 3]);
end
 
%% DC nav..

temp = squeeze(nav_kspace);
size(temp)
figure, colorful_plots(abs(temp(:,:,1,1,1)))

figure, colorful_plots(abs(nav_kspace_store(:,1:128,1)))
nks2 = reshape(nav_kspace_store, [size(nav_kspace_store,1), 2, size(nav_kspace_store,2)/2, size(nav_kspace_store,3)]);
nks2 = squeeze(cat(1, nks2(:,1,:,:), nks2(:,2,:,:)));

figure, colorful_plots(abs(nks2(:,:,1)))
figure, colorful_plots(abs(nks2(:,1:4096,1)))

temp2 = [squeeze(temp(:,32,:,:,1))];
size(temp2)
figure, 
subplot(1,2,1); colorful_plots(abs(temp2(:,:,1)))
subplot(1,2,2); colorful_plots(abs(temp2(:,:,2)))

%% remove OS
OverSampling = ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.x./ismrmrd_s.encoding.reconSpace.fieldOfView_mm.x;
if OverSampling == 2
    
    xdim = ismrmrd_s.encoding.encodedSpace.matrixSize.x;
    temp = reshape(cCoil_imgs, [xdim prod(size(cCoil_imgs))/xdim]);
    temp = temp( (xdim/4): (xdim/4 + xdim/2 -1), : );
    new_dims = size(cCoil_imgs); new_dims(1) = ismrmrd_s.encoding.reconSpace.matrixSize.x;
    cCoil_imgs = reshape(temp, new_dims);
    
end

%% Roemer coil combination "slice-by-slice"
% if ismrmrd_s.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 == 1
% dims = size(cCoil_imgs);
% CCM_img =  zeros(dims(1:end-1));
% for par = 1:pe2
%     for slc = 1:slices
%         for coc = 1:contrasts
%             for phc = 1:phases
%                 for repc = 1:reps
%                     for setc = 1:sets
%                         
%                         temp = cCoil_imgs(:,:,par,1,slc,coc,phc,repc,setc,:);
%                         csm = ismrm_estimate_csm_walsh( squeeze( temp ) );
%                         ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
%                         CCM_img(:,:,par,1,slc,coc,phc,repc,setc) = abs( sum( squeeze( temp ) .* ccm_roemer_optimal, 3) );
%                         
%                     end
%                 end
%             end
%         end
%     end
% end
% else
%     CCM_img = cCoil_imgs;
% end

% Sum-of-Squares coil combine
 SOS_img = squeeze(sqrt(sum(cCoil_imgs.*conj(cCoil_imgs),length(size(cCoil_imgs)) )));
 CCM_img = squeeze(SOS_img);

%% Pre-scan normalise
prescan_normalise = 0;
if prescan_normalise
    %% example prescan normalise..
    make_nhlbi_toolbox; make_dev;
    % nhlbi_toolbox.edit_toolbox()
    
    mid_slice = RR_slice_order(slices, round(slices/2));
    %     mid_slice = 9; % axial LF
    % %     mid_slice = 24; % axial Rad
    
    ii = find(raw_data.head.idx.slice+1 == mid_slice, 1, 'first');
    
    [psn, n_ismrmrd_s] = nhlbi_toolbox.prescan_normalize(nfile, ismrmrd_s, nhlbi_toolbox);
    im_fov =  [ismrmrd_s.encoding.reconSpace.fieldOfView_mm.x ismrmrd_s.encoding.reconSpace.fieldOfView_mm.y ismrmrd_s.encoding.reconSpace.fieldOfView_mm.z];
    im_res = im_fov./[ismrmrd_s.encoding.reconSpace.matrixSize.x ismrmrd_s.encoding.reconSpace.matrixSize.y ismrmrd_s.encoding.reconSpace.matrixSize.z];
    
    %     load('d20181113_NV1_brain.mat', 'sag_s24');
    %     img_temp = sag_s24(:,:,8);
    
    % load('d20181116_NV_brain.mat','sag_s24_ss'); whos sag_s24_ss
    % img_temp = sag_s24_ss;
    
    img = squeeze(CCM_img);
    img_temp = img(:,:,mid_slice);
    
    img_temp = rot90(img_temp,2); % ax LF
    % %     img_temp = rot90(img_temp,1);
    
    figure, imshow([dev.nrr(img_temp)],[0 3])
    
    % == sample to image FOV ==
    %     implay_RR(psn.psn_masked, [0 100])
    %     implay_RR(permute(psn.psn_masked,[2 3 1]), [0 100])
    %%  2D pre-scan normalise correction..
    
    % match Pre-scan normalise slice to imaging slice..
    % fudge the offsets until it looks right >>
    sample_x = round( size(psn.psn_image,1)*([psn.fov(1)/2 + [-im_fov(1)/2:im_fov(1)/2]]/psn.fov(1)));
    sample_y = round( size(psn.psn_image,2)*([psn.fov(2)/2 + [-im_fov(2)/2:im_fov(2)/2]]/psn.fov(2)));
    sample_z = size(psn.psn_image,3)/2;
    
    sample_shift = round(raw_data.head.position(:,ii)'./(psn.res/2)); % sagittal?
    %     sample_x = sample_x - sample_shift(1);
    %     sample_y = sample_y - sample_shift(1);
    %     sample_z = sample_z + sample_shift(2);
    sample_x = sample_x + 0;
    sample_y = sample_y + 1;
    sample_z = sample_z + 2;
    
    psn_slice = psn.psn_masked(sample_x, sample_y, sample_z); % psn_slice = psn_slice./max(psn_slice(:));
    psn_slice = rot90(psn_slice,2);
    figure, imshow(psn_slice,[])
    
    % == resample to image dimensions ==
    psn_slice_rs = zeros(size(img_temp,1),size(img_temp,2));
    
    psn_slice_rs(size(psn_slice_rs,1)/2 + ([1:length(psn_slice)]-round(length(psn_slice)/2)), size(psn_slice_rs,2)/2 + ([1:length(psn_slice)]-round(length(psn_slice)/2))) = ...
        ismrm_transform_image_to_kspace(psn_slice);
    psn_slice_rs = abs(ismrm_transform_kspace_to_image(psn_slice_rs));
    
    %     figure, imshow(psn_slice_rs,[])
    mask = imclose(bwareaopen(img_temp > mean(img_temp(:)) - 0.5*std(img_temp(:)),11),strel('disk',3)); %figure, imshow(mask)
    
    % == polynomial fit ==
    matrixSize = size(psn_slice_rs);
    x = 1:matrixSize(1); y = 1:matrixSize(2);
    [X, Y] = meshgrid(y,x); clear XY a;
    XY(:,:,1) = X; XY(:,:,2) = Y;
    temp=XY(:,:,1); a(:,1) = temp(find(mask));
    temp=XY(:,:,2); a(:,2) = temp(find(mask));
    
    opts= optimset('Display', 'none');
    surfit_v = @(B,a) B(1)*a(:,1).^2 + B(2)*a(:,2).^2 + B(3)*a(:,1).*a(:,2) + B(4)*a(:,1) + B(5)*a(:,2) + B(6);
    pn_factor = max(psn_slice(:));
    %    lsqcurvefit(func handle,   initial guess,      ROI,    original data,      fit minimum,                        fit maximum)
    B =  lsqcurvefit(surfit_v,      [0 0 0 0 0 pn_factor],   a,      psn_slice_rs(find(mask)),  [-10 -10 -10 -10 -10 0]*pn_factor, [10 10 10 10 10 1]*pn_factor, opts);
    
    surfit_s = @(B,XY) B(1)*XY(:,:,1).^2 + B(2)*XY(:,:,2).^2 + B(3)*XY(:,:,1).*XY(:,:,2) + B(4)*XY(:,:,1) + B(5)*XY(:,:,2) + B(6);
    %     coil_sens = surfit_s(B,XY);
    coil_sens = surfit_s(B,XY) + mean(img_temp(:)); % mean to offset low values.. needs optimisation!
    
    figure, imshow([psn_slice_rs coil_sens],[]); hold on; contour([mask mask], 'w')
    figure, imshow(dev.nrr([img_temp./coil_sens]),[0 3])
    
    %% ex. apply to image space
    
    % image rotation matrix
    IM_R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ];
    
    % image in physical coordinates
    Xinit = linspace(-(im_fov(1)/2 - im_res(1)/2), (im_fov(1)/2 - im_res(1)/2), ismrmrd_s.encoding.reconSpace.matrixSize.x);
    Yinit = linspace(-(im_fov(2)/2 - im_res(2)/2), (im_fov(2)/2 - im_res(2)/2), ismrmrd_s.encoding.reconSpace.matrixSize.y);
    Zinit = 0; % ??
    
    % Set up meshgrid for rotation
    [Xq,Yq,Zq] = meshgrid(Xinit,Yinit,Zinit);
    XYZ = [Xq(:) Yq(:) Zq(:)];
    
    % Rotate to image plan
    rotXYZ=XYZ*IM_R';
    % Add imaging offsets
    rotXYZ=rotXYZ + repmat(raw_data.head.position(:,mid_slice)',[size(rotXYZ,1) 1]);
    % Reshape for coil model
    Xqr = reshape(rotXYZ(:,1), size(Xq,1), []);
    Yqr = reshape(rotXYZ(:,2), size(Yq,1), []);
    Zqr = reshape(rotXYZ(:,3), size(Zq,1), []);
    XYZ = cat(4, Xqr, Yqr, Zqr); % figure, imshow([Xqr Yqr Zqr],[]); colormap('parula')
    
    % Generate coil model for imaging slice
    coil_sens = psn.model_3D( psn.model_3D_coefs , XYZ);
    
    figure('Name', 'Slice Coil Map'), imshow(coil_sens,[]);
    figure('Name', 'Normalised to model'), imshow(dev.nrr([img_temp./(coil_sens)]),[0 3])
    figure('Name', 'Normalised + mean(img) offset'), imshow(dev.nrr([img_temp./(coil_sens + mean(img_temp(:)) )]),[0 3])
    
    coil_sens3 = psn.model3_3D( psn.model3_3D_coefs , XYZ);
    
    figure('Name', 'Slice Coil Map'), imshow(coil_sens3,[]);
    figure('Name', 'Normalised to model'), imshow(dev.nrr([img_temp./(coil_sens3)]),[0 3])
    figure('Name', 'Normalised + mean(img) offset'), imshow(dev.nrr([img_temp./(coil_sens3 + mean(img_temp(:)) )]),[0 3])
    
    % %     % Check :: what is the corresponding PSN slice
    % %     Zqrp = round(Zqr./psn.res(1)) + size(psn.psn_image,1)/2;
    % %     Yqrp = round(Yqr./psn.res(2)) + size(psn.psn_image,2)/2;
    % %     Xqrp = round(Xqr./psn.res(3)) + size(psn.psn_image,3)/2; % figure, imshow([Xqrp Yqrp Zqrp],[]); colormap('parula')
    % % % %
    % % % %     length(unique(Zqr(1,:))),
    % % % %     length(unique(Zqr(:,1)))
    % %
    % % %     round(Zqr(1,:)./psn.psn_image(1)) + size(psn.psn_image,1)/2;
    % % %     round(Yqr(:,1)'./psn.psn_image(2)) + size(psn.psn_image,2)/2;
    % %
    % %    psni = psn.psn_image(Zqrp(1,:),Yqrp(:,1),Xqrp(1));
    % %    figure, imshow(psni,[])
    
end

%% Return variables
% kspace = squeeze(kspace);
img = squeeze(CCM_img);
montage_RR(img(:,:,:,1));
% new_dims = size(img);
% montage_RR(reshape(img, [new_dims(1) new_dims(2) prod(new_dims)/(new_dims(1)*new_dims(2))]));

img_s.img = img;
img_s.header = header;

if SNR_flag
    img_s.snr = snr;
end

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

ncoils = size(crt_k,ndims(crt_k));
img_pr = zeros(size(crt_k,1), size(crt_k,2), pseudo_reps);

for i = 1:pseudo_reps
    RR_loop_count(i,pseudo_reps);
    data_pr = crt_k + complex(randn(size(crt_k)),randn(size(crt_k)));
    data_pr = squeeze(mean(data_pr,av_dim));
    x = ismrm_transform_kspace_to_image(data_pr,[1 2]);
    
    csm = ismrm_estimate_csm_walsh(x);
    ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(ncoils));
    
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
% subplot(1,2,2); imshow(std(img_pr,[1 2],3),[]); colorbar; colormap('parula')
%
% SNR_s.img_pr = img_pr;
% SNR_s.snr = snr;
% SNR_s.gr = g;

end