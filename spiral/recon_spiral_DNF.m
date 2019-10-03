function [RTFM_output] = recon_spiral_DNF(dfile, nfile, user_opts)
% recon_spiral_RTFM(dfile, nfile, <user_opts>)
% 
% Calls nhlbi_toolbox & scratch_functions
%

%% Make NHLBI tools and setup paths
make_nhlbi_toolbox;
make_dev;

addpath('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Ramasawmy\local_MATLAB\ismrm_sunrise_matlab-master\irt\mex\v7');
addpath('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Ramasawmy\local_MATLAB\ismrm_sunrise_matlab-master\irt\mri');
addpath('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\test_data\2D_spiral'); % RR
addpath('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\test_data\3D_spiral'); % RR
addpath('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\trajectories'); % RR
addpath('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\spiral_GRAPPA'); % RR


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
user_opts.SNR = 0;
user_opts.GRAPPA = 0;
user_opts.GRAPPA_ref_dfile = ''; % uigetfile?

%%% Recipe selection
% user_opts.traj_file = 'meas_MID00221_FID16805_VDSFLMEASVE_24x48_192m_20190531.mat'; 
% user_opts.gridReadHigh = 1231;
user_opts.traj_file = 'meas_MID00279_FID37048_VDSFLMEASVE_m256_fov300_s48.mat'; 
user_opts.gridReadHigh = 2155; % end-of-spiral
% user_opts.gridReadHigh = 2632; % rewinder-included
% user_opts.GRAPPA_Xall_file = []; 
user_opts.GRAPPA_Xall_file = '\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\trajectories\Xall_fov300_matrix256_AF2_VD24x48.mat';


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
if user_opts.GRAPPA == 1
cal_data = h5read(user_opts.GRAPPA_ref_dfile,'/dataset/data');
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

kspace = permute(kspace, [1 2 4 3]);

%% Load traj and crop 
 
gridReadLow = 20;
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

if user_opts.GRAPPA == 1
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
    
else % GRAPPA == 1
    % === under-sampled data ===
    
    rseg = 8; %read segment size. default=8
    pseg = 1; %projection segment size. default=4
    af = interleaves_FS/interleaves;
    if isempty(user_opts.GRAPPA_Xall_file)
        Xall = prepXall(kx, ky, af); size(Xall)
        save(['\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Franson\Recon_code\trajectories\Xall_fov' num2str(fov) '_matrix' num2str(matrix) '_AF' num2str(af) '_VD24x48.mat'], 'Xall');
    else
        load(user_opts.GRAPPA_Xall_file);
    end
    
    disp('performing calibration');
    mode='calibration';
    
    tic
    [~,g_weights,source_indices] = throughtime_spiral_grappa_v4_crop(Xall,kx,ky,[],kspaceCalibration,rseg,pseg,af,mode,[],[]);
    toc
    
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
    
end

% whos imrec
montage_RR(imrec)
RTFM_output.img = imrec; 

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