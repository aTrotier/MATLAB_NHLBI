function [RTFM_output] = recon_spiral_DNF(dfile, nfile, user_opts)
% recon_spiral_RTFM(dfile, nfile, <user_opts>)
% 
% Calls nhlbi_toolbox & scratch_functions
%

%% Make NHLBI tools and setup paths
make_nhlbi_toolbox;
make_dev;

dfile = nhlbi_toolbox.run_path_on_sys(dfile);
nfile = nhlbi_toolbox.run_path_on_sys(nfile);

addpath(nhlbi_toolbox.run_path_on_sys('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Ramasawmy\local_MATLAB\ismrm_sunrise_matlab-master\irt\mex\v7'));
addpath(nhlbi_toolbox.run_path_on_sys('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Ramasawmy\local_MATLAB\ismrm_sunrise_matlab-master\irt\mri'));
addpath(nhlbi_toolbox.run_path_on_sys('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\test_data\2D_spiral')); % RR
addpath(nhlbi_toolbox.run_path_on_sys('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\test_data\3D_spiral')); % RR
addpath(nhlbi_toolbox.run_path_on_sys('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\trajectories')); % RR
addpath(nhlbi_toolbox.run_path_on_sys('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\spiral_GRAPPA')); % RR


%% Grab XML header
iRD_s = nhlbi_toolbox.h5read_xml(dfile);
RTFM_output.mrd_header = iRD_s;

RTFM_output.timestamp = datetime(now, 'ConvertFrom', 'datenum'); disp(' ');
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
disp(['Reconstructing: ' iRD_s.measurementInformation.protocolName]); disp(' ');
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');

%% Load data
raw_data = h5read(dfile,'/dataset/data');

%% user_opts default:
% need to configure to avoid overwriting input selection
if exist('user_opts','var')
    user_opts_in = user_opts;
else
    user_opts_in = struct();
end
user_opts_in_list = fieldnames(user_opts_in);

% Future work: toggle recons
user_opts.bin = 0;
user_opts.SNR = 0;
user_opts.GRAPPA = 0;
user_opts.GRAPPA_ref_dfile = ''; % uigetfile?
user_opts.arm_blocks = 1;
user_opts.import_weights = [];
user_opts.export_weights = 0;

%%% Recipe selection
% user_opts.traj_file = 'meas_MID00221_FID16805_VDSFLMEASVE_24x48_192m_20190531.mat'; 
% user_opts.gridReadHigh = 1231;
user_opts.traj_file = 'meas_MID00279_FID37048_VDSFLMEASVE_m256_fov300_s48.mat'; 
user_opts.gridReadLow = 20;
user_opts.gridReadHigh = 2155; % end-of-spiral
% user_opts.gridReadHigh = 2632; % rewinder-included
user_opts.GRAPPA_Xall_file = []; 
% user_opts.GRAPPA_Xall_file = '\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\trajectories\Xall_fov300_matrix256_AF2_VD24x48.mat';
% user_opts.GRAPPA_Xall_file = '\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\trajectories\meas_MID00359_FID42470_VDSFLMEASVE_fov320_s48_dt.mat';


%% Grab sequence-specific info ### Requires user attention! ### 
% === update user_opts with input data === 
for i = 1:length(user_opts_in_list)
    if isfield(user_opts, user_opts_in_list{i})
        user_opts.(matlab.lang.makeValidName(user_opts_in_list{i})) = user_opts_in.(matlab.lang.makeValidName(user_opts_in_list{i}));
    end
end


disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
disp('User options:');
disp(user_opts);
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');

clear user_opts_in user_opts_in_list

%% NOISE PREWHITEN
if ~exist('nfile', 'var')
    nfile = [];
end
dmtx = nhlbi_toolbox.noise_adjust(nfile, iRD_s, raw_data.head.sample_time_us(1)*1e-6, nhlbi_toolbox);

%% GRAPPA CALIBRATION DATA
if user_opts.GRAPPA == 1 && isempty(user_opts.import_weights)
cal_data = h5read(nhlbi_toolbox.run_path_on_sys(user_opts.GRAPPA_ref_dfile),'/dataset/data');
% cal_dmtx = ismrm_dmtx_RR(user_opts.GRAPPA_ref_nfile, cal_data.head.sample_time_us(1)*1e-6); % NOISE FILE MUST BE THE SAME!

samples     = double(cal_data.head.number_of_samples(1));
interleaves_FS = (1 + double(max(cal_data.head.idx.kspace_encode_step_1)));
channels    = double(cal_data.head.active_channels(1));
reps_CAL        = (1 + double(max(cal_data.head.idx.repetition)));

kspaceCalibration = complex(zeros([samples interleaves_FS reps_CAL channels],'single'));

for ii = 1:length(cal_data.data)
    
    d1 = cal_data.data{ii};
    d2 = complex(d1(1:2:end), d1(2:2:end));
    d3 = reshape(d2, samples, channels); %  RE & IM (2)
    
    d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
    
    
    kspaceCalibration(:,...
        cal_data.head.idx.kspace_encode_step_1(ii)+1, ...
        cal_data.head.idx.repetition(ii)+1, ...
        :) = d3;
    
end

kspaceCalibration = permute(kspaceCalibration, [1 2 4 3]);
% size(kspaceCalibration)
% [samples,interleaves_FS,channels,reps_CAL]=size(kspaceCalibration);
end

%%  Grab general info
% Future: switch to generic script set-up
interleaves     = (1 + double(max(raw_data.head.idx.kspace_encode_step_1)));
reps      = (1 + double(max(raw_data.head.idx.repetition)));
samples     = double(raw_data.head.number_of_samples(1));
channels    = double(raw_data.head.active_channels(1));
% pe2         = 1+double(max(raw_data.head.idx.kspace_encode_step_2));
% contrasts   = (1 + double(max(raw_data.head.idx.contrast)));
% phases      = (1 + double(max(raw_data.head.idx.phase)));
% sets        = (1 + double(max(raw_data.head.idx.set)));
% averages    = (1 + double(max(raw_data.head.idx.average)));
% slices      = (1 + double(max(raw_data.head.idx.slice)));

% rep_vec  = (1 + double(raw_data.head.idx.repetition));
% sets_vec = (1 + double(raw_data.head.idx.set));
% ksp1_vec = (1 + double(raw_data.head.idx.kspace_encode_step_1));

% visual debugging;
%     nhlbi_toolbox.plot_experiment(raw_data)

matrix = iRD_s.encoding.reconSpace.matrixSize.x;
dt = raw_data.head.sample_time_us(1)*1e-6;
matrix_size = [matrix matrix]; user_opts.matrix_size = matrix_size; 
   
fov = iRD_s.encoding.reconSpace.fieldOfView_mm.x;
zf_inplane = 0;
    
%% Grab data 

kspace = complex(zeros([samples interleaves reps channels],'single'));

for ii = 1:length(raw_data.data)
    
    d1 = raw_data.data{ii};
    d2 = complex(d1(1:2:end), d1(2:2:end));
    d3 = reshape(d2, samples, channels); %  RE & IM (2)
    
    d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
    
    kspace(:,...
        raw_data.head.idx.kspace_encode_step_1(ii)+1, ...
        raw_data.head.idx.repetition(ii)+1, ...
        :) = d3;
    
end

if user_opts.arm_blocks > 1
    
    kspace_temp = complex(zeros([samples length(raw_data.data) channels],'single'));
    kspace_b = complex(zeros([samples interleaves ((reps+1)/2 -2) channels],'single')); % size(kspace_b)
    sample_indices = 1:length(raw_data.data);
    
    for ii = 1:length(raw_data.data)
        
        d1 = raw_data.data{ii};
        d2 = complex(d1(1:2:end), d1(2:2:end));
        d3 = reshape(d2, samples, channels); %  RE & IM (2)
        
        d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
        
        kspace_temp(:,ii,:) = d3;
        
    end
    
    kspace_temp = kspace_temp(:,(2*interleaves + 1):end,:);
    sample_indices = sample_indices(:,(2*interleaves + 1):end,:);
    
    averages = ((reps+1)/2 -2)/user_opts.arm_blocks;
%     kspace_temp = reshape(kspace_temp, [samples, interleaves/user_opts.arm_blocks user_opts.arm_blocks*((reps+1)/2 -2)/averages averages channels]);
    kspace_temp = reshape(kspace_temp, [samples, interleaves/user_opts.arm_blocks user_opts.arm_blocks*user_opts.arm_blocks averages channels]);
%     size(kspace_temp)
    sample_indices = reshape(sample_indices, [interleaves/user_opts.arm_blocks user_opts.arm_blocks*user_opts.arm_blocks averages ]);
    
    repc = 0;
    for jj = 1:averages
        for ii = 1:user_opts.arm_blocks
            repc = repc + 1;
            temp_vec = ii:(user_opts.arm_blocks):(user_opts.arm_blocks*user_opts.arm_blocks);
             temp_ints= vec(sample_indices(:,temp_vec,jj));
% %              [temp_ints raw_data.head.idx.kspace_encode_step_1(temp_ints)+1]
            kspace_b(:,raw_data.head.idx.kspace_encode_step_1(temp_ints)+1,repc,:) = reshape(squeeze(kspace_temp(:,:,temp_vec,jj,:)), [samples interleaves channels]);
            
        end
    end

    kspace_b = permute(kspace_b, [1 2 4 3]); clear kspace_temp sample_indices 
    
end

kspace = permute(kspace, [1 2 4 3]);

%% Load traj and crop 
 
gridReadLow = user_opts.gridReadLow;
gridReadHigh = user_opts.gridReadHigh;

disp(['Loading ' user_opts.traj_file])
whos('-file',user_opts.traj_file);
load(user_opts.traj_file);
if gridReadHigh > size(kxall,1) %#ok<*NODEF>
    warning('gridReadHigh is too high! Using length traj');
    gridReadHigh = size(kxall,1);
end

kxall = kxall(gridReadLow+1:gridReadHigh,:);
kyall = kyall(gridReadLow+1:gridReadHigh,:);

% Crop kspace & calibration data
kspace = kspace(gridReadLow+1:gridReadHigh,:,:,:);

if user_opts.arm_blocks > 1
    kspace_b = kspace_b(gridReadLow+1:gridReadHigh,:,:,:);
end

if user_opts.GRAPPA == 1 && isempty(user_opts.import_weights)
kspaceCalibration = kspaceCalibration(gridReadLow+1:gridReadHigh,:,:,:);
end

% Scale trajectory
[nr,ns] = size(kxall);
traj = zeros(nr,ns,2);
traj(:,:,1) = kxall;
traj(:,:,2) = kyall;
traj = reshape(traj,nr*ns,2);
traj = traj*(pi/max(traj(:)));

% Prep density compensation weights
kx = kxall./max(kxall(:)).*matrix_size(1)/2;
ky = kyall./max(kyall(:)).*matrix_size(2)/2;
ksp = [kx(:) ky(:)];
ksp = ksp./fov;
mask = true(matrix_size(1),matrix_size(2));
sizeMask = size(mask);
nufft_args = {sizeMask, [6 6], 2*sizeMask, sizeMask/2, 'table', 2^12, 'minmax:kb'};
G = Gmri(ksp, mask, 'fov', fov, 'basis', {'dirac'}, 'nufft', nufft_args);
wi = abs(mri_density_comp(ksp, 'pipe','G',G.arg.Gnufft));

% Prepare NUFFT operator
st = nufft_init(traj, matrix_size, [6 6], matrix_size.*2, matrix_size./2);
    
    
%% Perform 2D recon
if user_opts.GRAPPA == 0
    % === fully sampled data ===
    
    [nr,ns,nc,nframes] = size(kspace);
    
    % Grid to image
    kspace = reshape(kspace,nr*ns,nc,nframes);
    % im_recon = zeros(matrix_size(1),matrix_size(2),nc,nframes);
    imrec = zeros(matrix_size(1), matrix_size(2),nframes);
    for f = 1:nframes
        
        x = nufft_adj(kspace(:,:,f).*repmat(wi,[1,nc]),st);
        imrec(:,:,f) = sqrt(sum(x.*conj(x),3));
    end
    
    if user_opts.bin == 1
        [bin_data, RTFM_output.NominalInterval] = physio_Binning(raw_data.head.physiology_time_stamp, 30); % figure, plot(raw_data.head.physiology_time_stamp(1,:), 'r-');
        
        if RTFM_output.NominalInterval == 0
            warning('No ECG trace!'); % and continue
            
        else
            kspace = reshape(kspace,[nr ns nc nframes]);
            kspace = permute(kspace, [1 2 4 3]);
            bin_temp = zeros([matrix_size length(bin_data)]);
            
            % ==
            
            csm_data = squeeze(mean(kspace, 3));
            csm_data = reshape(csm_data, [nr*ns, nc]);
            x = nufft_adj(csm_data.*repmat(wi,[1,nc]),st);
            csm  = ismrm_estimate_csm_mckenzie(squeeze(x)); % montage_RR(abs(csm));
            
            kspace = reshape(kspace, [nr ns*nframes nc]);
            
            
            % ==
            
            for i = 1:length(bin_data)
                RR_loop_count(i, length(bin_data));
                sample_window = bin_data{i};
                
                % Set-up plan
                
                xi = raw_data.head.idx.kspace_encode_step_1(sample_window)+1;
                data = kspace(:,sample_window,:);
                
                data2 = zeros(nr, ns, nc);
                
                for j = 1:ns
                    nnn = length(find(xi==j));
                    xitemp(j) = nnn;
                    data2(:,j,:) = mean(data(:,xi==j,:),2);
                end
                xi = find(xitemp);
                data = data2(:,xi,:);
                
                data = reshape(data,length(xi)*nr,channels);
                
                kx = kxall(:,xi)./max(kxall(:)).*matrix_size(1)/2;
                ky = kyall(:,xi)./max(kyall(:)).*matrix_size(2)/2;
                ksp = [kx(:) ky(:)];
                ksp = ksp./fov;
                mask = true(matrix_size(1),matrix_size(2));
                sizeMask = size(mask);
                nufft_args = {sizeMask, [6 6], 2*sizeMask, sizeMask/2, 'table', 2^12, 'minmax:kb'};
                G = Gmri(ksp, mask, 'fov', fov, 'basis', {'dirac'}, 'nufft', nufft_args);
                wi = abs(mri_density_comp(ksp, 'pipe','G',G.arg.Gnufft));
                
                % Prepare NUFFT operator
                clear traj;
                traj(:,:,1) = kxall(:,xi);
                traj(:,:,2) = kyall(:,xi);
                traj = reshape(traj,nr*length(xi),2);
                traj = traj*(pi/max(traj(:)));
                st = nufft_init(traj, matrix_size, [6 6], matrix_size.*2, matrix_size./2);
               
%                 xg = nufft_adj(data.*repmat(wi,[1,channels]), st);
                
                I = ones(matrix_size(1)); D = repmat(wi, [1 channels]);
                x = cg_RR(data(:), st, I, D, csm, wi, 3);
                
                [bin_temp(:,:,i)] = sqrt(sum(x.*conj(x),3));
            end
            
            RTFM_output.bin_imrec = bin_temp;
        end
        
    end
    
else % GRAPPA == 1
    % === under-sampled data ===
    
    
    if ~isempty(user_opts.import_weights)
        disp('Importing calibration');
        
        g_weights       = user_opts.import_weights.g_weights;
        source_indices  = user_opts.import_weights.source_indices;
        Xall            = user_opts.import_weights.Xall;
        interleaves_FS  = user_opts.import_weights.interleaves_FS;
        af = interleaves_FS/interleaves;
    
    else
        
        rseg = 8; %read segment size. default=8
        pseg = 1; %projection segment size. default=4
        af = interleaves_FS/interleaves;
        
        if isempty(user_opts.GRAPPA_Xall_file)
            Xall = prepXall(kx, ky, af); size(Xall)
            %save(['\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\trajectories\Xall_fov' num2str(fov) '_matrix' num2str(matrix) '_AF' num2str(af) '_VD24x48.mat'], 'Xall');
        else
            load(user_opts.GRAPPA_Xall_file);
        end
        
        disp('performing calibration');
        mode='calibration';
        
        tic
        [~,g_weights,source_indices] = throughtime_spiral_grappa_v4_crop(Xall,kx,ky,[],kspaceCalibration,rseg,pseg,af,mode,[],[]);
        toc
    end
    
    disp('performing reconstruction')
    mode='reconstruction';
    
    if user_opts.SNR
        
        disp('SNR calculation')
        
        pseudo_reps = 100; disp(['Running ' num2str(pseudo_reps) ' pseudo-reps']);
        
        img_pr = zeros([matrix_size pseudo_reps]);
        
        for iPR = 1:pseudo_reps
            RR_loop_count(iPR,pseudo_reps);
            
            % == add white noise ==
            kspaceUnderSampled_wn = kspace(:,:,:,reps); % run on last repetition. 
            kspaceUnderSampled_wn = kspaceUnderSampled_wn + complex(randn(size(kspaceUnderSampled_wn)),randn(size(kspaceUnderSampled_wn)));
            
            % == perform recon ==
            sigrecMeas = throughtime_spiral_grappa_v4_crop(Xall,kx,ky,kspaceUnderSampled_wn,[],[],[],af,mode,g_weights,source_indices);
            kfilled = reshape(sigrecMeas,size(kx,1)*interleaves_FS,channels);
            x = nufft_adj(kfilled.*repmat(wi,[1,channels]),st); 
            %    coil_imrec(:,:,:,m) = i;
            
            % == Use Roemer’s linear coil combine ==
            csm = ismrm_estimate_csm_walsh(x); montage_RR(x)
            ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened data
            
            img_pr(:,:,iPR) = abs( sum( x .* ccm_roemer_optimal , 3) );
        end
        
        kfilled = reshape(mean(kspaceCalibration,4),size(kx,1)*interleaves_FS,channels);
        x = nufft_adj(kfilled.*repmat(wi,[1,channels]),st); 
        cal_scaling = sqrt(sum(x.*conj(x),3)); montage_RR(cal_scaling)
        
        % == calculate g-factor and SNR ==
        g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
        g(g < eps) = 1;
        SNR = mean(img_pr,3)./g;
                
        csm_sq = sum(csm .* conj(csm),3); csm_sq(csm_sq < eps) = 1;
        g = g .* sqrt(csm_sq);
    
        figure, subplot(1,2,1); imagesc(SNR,[0 50]); title('SNR map');
        subplot(1,2,2);         imagesc(g,[1 2]);    title('g-factor map');
    end
    
    coil_imrec = zeros([matrix_size channels reps]);
    imrec = zeros([matrix_size reps]);
    
    for m=1:reps
        RR_loop_count(m, reps);
        
        sigrecMeas = throughtime_spiral_grappa_v4_crop(Xall,kx,ky,kspace(:,:,:,m),[],[],[],af,mode,g_weights,source_indices);
        kfilled = reshape(sigrecMeas,size(kx,1)*interleaves_FS,channels);
        x = nufft_adj(kfilled.*repmat(wi,[1,channels]),st);
        coil_imrec(:,:,:,m) = x;
        imrec(:,:,m) = sqrt(sum(x.*conj(x),3));
        
    end
    
    if user_opts.arm_blocks > 1
        coil_imrec_b = zeros([matrix_size channels size(kspace_b,4)]);
        imrec_b = zeros([matrix_size size(kspace_b,4)]);
        
        for m=1:size(kspace_b,4)
            RR_loop_count(m, size(kspace_b,4));
            
            sigrecMeas = throughtime_spiral_grappa_v4_crop(Xall,kx,ky,kspace_b(:,:,:,m),[],[],[],af,mode,g_weights,source_indices);
            kfilled = reshape(sigrecMeas,size(kx,1)*interleaves_FS,channels);
            x = nufft_adj(kfilled.*repmat(wi,[1,channels]),st);
            coil_imrec_b(:,:,:,m) = x;
            imrec_b(:,:,m) = sqrt(sum(x.*conj(x),3));
            
        end
        
        RTFM_output.img_b = imrec_b;
    end
end

%%
prescannorm = 0;
if prescannorm
    
    [psn, n_ismrmrd_s] = nhlbi_toolbox.prescan_normalize(nfile, iRD_s, nhlbi_toolbox);
    im_fov =  [iRD_s.encoding.reconSpace.fieldOfView_mm.x iRD_s.encoding.reconSpace.fieldOfView_mm.y iRD_s.encoding.reconSpace.fieldOfView_mm.z];
    im_res = im_fov./[iRD_s.encoding.reconSpace.matrixSize.x iRD_s.encoding.reconSpace.matrixSize.y iRD_s.encoding.reconSpace.matrixSize.z];
    
    img_temp = imrec(:,:,end);
    
    img_temp = rot90(img_temp,0); % ax LF
    % %     img_temp = rot90(img_temp,1);
    
    figure, imshow([dev.nrr(img_temp)],[0 3])

    %%  
    % image rotation matrix
    IM_R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ];
    
    % image in physical coordinates
    Xinit = linspace(-(im_fov(1)/2 - im_res(1)/2), (im_fov(1)/2 - im_res(1)/2), iRD_s.encoding.reconSpace.matrixSize.x);
    Yinit = linspace(-(im_fov(2)/2 - im_res(2)/2), (im_fov(2)/2 - im_res(2)/2), iRD_s.encoding.reconSpace.matrixSize.y);
    Zinit = 0; % ??
    
    % Set up meshgrid for rotation
    [Xq,Yq,Zq] = meshgrid(Xinit,Yinit,Zinit);
    XYZ = [Xq(:) Yq(:) Zq(:)];
    
    % Rotate to image plan
    rotXYZ=XYZ*IM_R';
    % Add imaging offsets
    rotXYZ=rotXYZ + repmat(raw_data.head.position(:,1)',[size(rotXYZ,1) 1]);
    % Reshape for coil model
    Xqr = reshape(rotXYZ(:,1), size(Xq,1), []);
    Yqr = reshape(rotXYZ(:,2), size(Yq,1), []);
    Zqr = reshape(rotXYZ(:,3), size(Zq,1), []);
    XYZ = cat(4, Xqr, Yqr, Zqr); % figure, imshow([Xqr Yqr Zqr],[]); colormap('parula')
    
    % Generate coil model for imaging slice
    coil_sens = psn.model_3D( psn.model_3D_coefs , XYZ);
    
    figure, 
    subplot(2,4,1); hold on, title('Slice Coil Map'), imshow(coil_sens,[]);
    subplot(2,4,2); hold on, title('Normalised to model'), imshow(dev.nrr([img_temp./(coil_sens)]),[0 3])
    subplot(2,4,3); hold on, title('Normalised + mean(img) offset'), imshow(dev.nrr([img_temp./(coil_sens + mean(img_temp(:)) )]),[0 3])
    subplot(2,4,4); hold on, title('Nada'), imshow(dev.nrr([img_temp]),[0 3])
    
    coil_sens3 = psn.model3_3D( psn.model3_3D_coefs , XYZ);
    
    subplot(2,4,5); hold on, title('Slice Coil Map'), imshow(coil_sens3,[]);
    subplot(2,4,6); hold on, title('Normalised to model'), imshow(dev.nrr([img_temp./(coil_sens3)]),[0 3])
    subplot(2,4,7); hold on, title('Normalised + mean(img) offset'), imshow(dev.nrr([img_temp./(coil_sens3 + mean(img_temp(:)) )]),[0 3])
end
    
    
%%
% whos imrec
if (ispc) || (ismac) % dont draw in linux
    montage_RR(imrec);
end
   
RTFM_output.img = imrec;

if user_opts.export_weights == 1
    user_opts.import_weights.g_weights      = g_weights;
    user_opts.import_weights.source_indices = source_indices;
    user_opts.import_weights.Xall           = Xall;
    user_opts.import_weights.interleaves_FS = interleaves_FS;
end

if (~isempty(user_opts.import_weights) && (user_opts.export_weights == 0))
    user_opts.import_weights = []; % save space on data export
end

RTFM_output.user_opts = user_opts;

end

function Xall = prepXall(kx,ky, af)
[nr,np] =size(kx);
npacc=np/af;
pidxacq = 1:af:np;
pidxtarg = [1:1:np];
pidxtarg(pidxacq) = 0;
pidxtarg(pidxtarg==0) = [];
Xall = zeros(nr,np,6,2); % [read position, projection position, point 1-6, readout value/arm number]

for projind = 1:length(pidxtarg)
    proj = pidxtarg(projind);
    for read = 1:nr
        % Location of one target point
        locX = kx(read,proj);
        locY = ky(read,proj);
        
            % Find distance between target point and all acquired points
            diffX = kx - locX;
            diffY = ky - locY;
            difftot = diffX.^2 + diffY.^2;
            mask = zeros(size(difftot));
            mask(:,pidxacq) = 1;
            difftot = difftot.*mask;
            difftot(difftot==0) = max(difftot(:));
            source_indices = zeros(6,2);
            
            % Find first source point, and then use two surrounding readout
            % points
            [~,idx] = min(difftot(:));
            [x,y] = ind2sub(size(difftot),idx);
            source_indices(1,1) = x;
            source_indices(1,2) = y;
            difftot(:,y) = max(difftot(:));
            source_indices(2,1) = x-2;
            if((x-2)<1) % avoid selecting points outside of collected ones; current recon implementation can't handle this case
                source_indices(2,1) = x+1;
            end
            source_indices(3,1) = x+2;
            if((x+2)>nr)
                source_indices(3,1) = x-1;
            end
            source_indices(2,2) = y;
            source_indices(3,2) = y;
            
            % Zero out all arms except those immediately before and after y
            % identified above (find closest readout point on neighboring
            % arms)
            arms = 1:size(difftot,2);
            if(y<proj)
                difftot(:,find(arms<y)) = max(difftot(:))+1;
            else
                difftot(:,find(arms>y)) = max(difftot(:))+1;
            end
            
            % Zero out quadrants surrounding first target point
            if(diffX(x,y)<0) % source kx is less than target kx, zero out points to the left
                difftot(diffX<0) = max(difftot(:))+1;
            else
                difftot(diffX>0) = max(difftot(:))+1;
            end
            
            if(diffY(x,y)<0) % source ky is less than target kx, zero out points to the left
                difftot(diffY<0) = max(difftot(:))+1;
            else
                difftot(diffY>0) = max(difftot(:))+1;
            end
            
            % Find fourth source point, and then use two surrounding readout
            % points
            [~,idx] = min(difftot(:));
            [x,y] = ind2sub(size(difftot),idx);
            source_indices(4,1) = x;
            source_indices(4,2) = y;

            difftot(x,y) = max(difftot(:));
            source_indices(5,1) = x-2;
            if((x-2)<1) % avoid selecting points outside of collected ones; current recon implementation can't handle this case
                source_indices(5,1) = x+1;
            end
            source_indices(6,1) = x+2;
            if((x+2)>nr)
                source_indices(6,1) = x-1;
            end
            source_indices(5,2) = y;
            source_indices(6,2) = y;     
            
            if(y == proj)
               y = source_indices(1,2);
               x = source_indices(1,1);
               source_indices(4,1) = x-1;
               source_indices(5,1) = x+1;
               source_indices(6,1) = x-3;
               if((x+1)>nr)
                   source_indices(5,1) = x-4;
               end
               source_indices(4,2) = y;
               source_indices(5,2) = y;
               source_indices(6,2) = y;
            end         

        Xall(read,proj,:,:) = source_indices;
    end
end
end