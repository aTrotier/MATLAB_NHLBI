function [RTFM_output] = recon_spiral_RTFM(dfile, nfile, user_opts)
% recon_spiral_RTFM(dfile, nfile, <user_opts>)
% 
% Calls nhlbi_toolbox & scratch_functions
%

%% Make NHLBI tool "class" (dev)
make_nhlbi_toolbox;
make_dev;

dfile = nhlbi_toolbox.run_path_on_sys(dfile);
nfile = nhlbi_toolbox.run_path_on_sys(nfile);

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

user_opts.delayFactor = 0;

if ~isfield(user_opts,'vds')
    user_opts.vds = 100;
end

user_opts.iter = 3;
user_opts.number_cardiac_frames = 30;
user_opts.ave_for_csm = 4;
user_opts.IO = 0;
user_opts.SNR = 0;
% user_opts.sliding_window_width = 0;

% Future work: toggle recons
% user_opts.recon_CG_SENSE = 1;
% user_opts.recon_GRAPPA = 0;
% user_opts.recon_CS_TV = 0;

%%% Recipe selection
user_opts.recon_RT_frame    = 0;
user_opts.recon_RT_SW       = 0;
user_opts.recon_RB_ECG      = 1;
user_opts.recon_RB_DCSG     = 0;
user_opts.recon_RB_FSG      = 0;
user_opts.gridAll           = 0;

%% Grab sequence-specific info ### Requires user attention! ### 
% GASFib configuration?
disp(' ');disp('Sequence-specific configuration:'); disp(' ');

iRD_s.encoding.trajectoryDescription.userParameterLong
traj_setup.interleaves = iRD_s.encoding.trajectoryDescription.userParameterLong.value;
user_opts.sliding_window_width = traj_setup.interleaves;

traj_setup.gMax = iRD_s.encoding.userParameterDouble.value;
traj_setup.sMax = iRD_s.encoding.userParameterDouble_1.value;

if isfield(iRD_s.userParameters, 'userParameterLong')
    iRD_s.userParameters.userParameterLong
    venc = iRD_s.userParameters.userParameterLong.value;
else
    venc = 0;
end

if max(double(raw_data.head.idx.repetition)) == 0
    repvec = 0:(length(raw_data.head.idx.repetition)/2 - 1);
    raw_data.head.idx.repetition = reshape(repmat(repvec, [2 1]), [length((raw_data.head.idx.repetition)) 1]);
end

% === update user_opts with input data === 
for i = 1:length(user_opts_in_list)
    if isfield(user_opts, user_opts_in_list{i})
        user_opts.(matlab.lang.makeValidName(user_opts_in_list{i})) = user_opts_in.(matlab.lang.makeValidName(user_opts_in_list{i}));
    end
end

% pseudoRep = 0, no PR
%             n, perform PR on slice n
%             [], perform PR on all slices
if isempty(user_opts.SNR )
    pseudoRep = 1;
    user_opts.SNR  = 1:user_opts.number_cardiac_frames;
end
if user_opts.SNR > 0 % specify slice(s);
    if user_opts.SNR > user_opts.number_cardiac_frames
        error(['PR slice > #slices, using slice '])
    end
end


disp(' ');disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
disp('User options:');
disp(user_opts);  disp(' ');
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');

clear user_opts_in user_opts_in_list

%%  Grab general info
% Future: switch to generic script set-up
gangles     = (1 + double(max(raw_data.head.idx.kspace_encode_step_1)));
pe2         = 1+double(max(raw_data.head.idx.kspace_encode_step_2));
contrasts   = (1 + double(max(raw_data.head.idx.contrast)));
phases      = (1 + double(max(raw_data.head.idx.phase)));
samples     = double(raw_data.head.number_of_samples(1));
channels    = double(raw_data.head.active_channels(1));
sets        = (1 + double(max(raw_data.head.idx.set)));
gareps      = (1 + double(max(raw_data.head.idx.repetition)));
averages    = (1 + double(max(raw_data.head.idx.average)));
slices      = (1 + double(max(raw_data.head.idx.slice)));

% rep_vec  = (1 + double(raw_data.head.idx.repetition));
% sets_vec = (1 + double(raw_data.head.idx.set));
% ksp1_vec = (1 + double(raw_data.head.idx.kspace_encode_step_1));

% visual debugging;
%     nhlbi_toolbox.plot_experiment(raw_data)

matrix = iRD_s.encoding.reconSpace.matrixSize.x;
dt = raw_data.head.sample_time_us(1)*1e-6;
matrix_size = [matrix matrix]; user_opts.matrix_size = matrix_size; 

disp(' ');disp('### Experiment Dimensions ###');disp(' ');
Experiment_parameters = {'Samples','Interleaves','Rotations', 'Matrix', 'PE2', 'Averages', 'Slices', 'Contrasts', 'Phases', 'Reps', 'Sets', 'Channels'}';
Value = [samples traj_setup.interleaves gangles iRD_s.encoding.encodedSpace.matrixSize.x pe2 averages slices contrasts phases gareps sets channels]';
disp(table( Experiment_parameters,Value )); clear Experiment_parameters Value; disp(' ');

%% Noise Adjust
% iRD_s_noise = nhlbi_toolbox.h5read_xml(nfile);
% nhlbi_toolbox.check_noise_dependency(iRD_s, iRD_s_noise)
dmtx = nhlbi_toolbox.noise_adjust(nfile, iRD_s, raw_data.head.sample_time_us(1)*1e-6, nhlbi_toolbox);

%% RT data grab : queue all en-masse
data_store = complex(zeros([samples, length(raw_data.data), channels], 'single'));
for i = 1:length(raw_data.data)
%     RR_loop_count(i,gareps, 'Loading data')
    d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
    d = reshape(d, samples, 1, channels);
    d = ismrm_apply_noise_decorrelation_mtx(d, dmtx);
    data_store(:,i,:) = squeeze(d);
end

% % this method may be more memory-friendly for large datasets: (take only what is needed at the time)
% temp = cast(cat(1,raw_data.data{1:5}),'double');
% temp2 = complex(temp(1:2:end), temp(2:2:end));
% temp2 = reshape(temp2, [samples channels 5]);
% temp2 = permute(temp2,[1 3 2]);

%% Set-up spiral trajectories

FOV = iRD_s.encoding.reconSpace.fieldOfView_mm.x/10;
FOV = [FOV -1*FOV*(1 - user_opts.vds/100)]; disp(['FOV: ' num2str(FOV)])
krmax = 1/(2*(FOV(1)/matrix_size(1)));

[k,g] = vds(traj_setup.sMax, traj_setup.gMax, dt, traj_setup.interleaves, FOV, krmax); close;

% Check if trajectory plan "k" is same length as number of samples
if samples > length(k)
    samples2 = length(k);
else
    samples2 = samples;
end

trajectory_nominal = zeros(samples2,gangles,2);
gradients_nominal =  zeros(samples2,gangles,2);

neg = -1;
for solid_int= 1:gangles
    rot = (solid_int-1)*(2*pi/gangles);
    trajectory_nominal(:,solid_int,2) = neg*-( real(k(1:samples2)) *cos(rot) + imag(k(1:samples2)) *sin(rot));
    trajectory_nominal(:,solid_int,1) = neg*-(-real(k(1:samples2)) *sin(rot) + imag(k(1:samples2)) *cos(rot));
    gradients_nominal(:,solid_int,2)  = neg*-( real(g(1:samples2)) *cos(rot) + imag(g(1:samples2)) *sin(rot));
    gradients_nominal(:,solid_int,1)  = neg*-(-real(g(1:samples2)) *sin(rot) + imag(g(1:samples2)) *cos(rot));
end

if user_opts.IO == 0
    
    %%% GIRF corrections
    gradients_store = gradients_nominal;
    if user_opts.delayFactor > 0 % Option for delayed acquisition to avoid ADC ringdown
        spiral_start = floor(user_opts.delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
        gradients_nominal = cat(1,zeros([spiral_start gangles 2]), gradients_store(1:(samples2-spiral_start),:,:));
    end
    
    % set-up GIRF object
    % R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; % Rotation matrix
    R = [raw_data.head.phase_dir(:,1), raw_data.head.read_dir(:,1), raw_data.head.slice_dir(:,1)  ]; % Rotation matrix the right order..
    
    tRR = 0; % legacy object for sub-dwelltime delay correction
    sR.R = R; % pack rotation matrix and field info in to one object for backwards compatability..
    sR.T = iRD_s.acquisitionSystemInformation.systemFieldStrength_T;
    
    trajectory = apply_GIRF(gradients_nominal, dt, sR, tRR ); % trajectory_nominal = apply_GIRF(gradients_nominal, dt, R, tRR );
    trajectory = trajectory(:,:,1:2);
else
%     figure, hold on, plot(gradients_nominal(:,1,1),'b-');plot(gradients_nominal(:,1,2),'r-');
    
    trajectory_nominalr(:,:,1) = flipud(-1*trajectory_nominal(:,:,1));
    trajectory_nominalr(:,:,2) = flipud(-1*trajectory_nominal(:,:,2));
    trajectory = [trajectory_nominalr;trajectory_nominal];
    
    gradients_nominalr(:,:,1)= flipud(gradients_nominal(:,:,1));
    gradients_nominalr(:,:,2)= flipud(gradients_nominal(:,:,2));
    gradients_nominal = [gradients_nominalr;gradients_nominal];
    
%     plot(trajectory_nominalr(:,:,1), trajectory_nominalr(:,:,2), 'k-')
%     plot(trajectory_nominal(:,:,1), trajectory_nominal(:,:,2), 'r-')

    %%
    rr = 1;
    if rr > 1
        
        [x_m0, y_m0, traj_setup.grad_info] = vds_M0(traj_setup.sMax, traj_setup.gMax, dt, traj_setup.interleaves, FOV, krmax, 1); close;
        
        %   mismatch with vds and vds_m0..
        %     rewinder_grt(:,1) = x_m0((traj_setup.grad_info.n_spiral+1):end);
        %     rewinder_grt(:,2) = y_m0((traj_setup.grad_info.n_spiral+1):end);
        %
        %     interp_factor = samples/traj_setup.grad_info.n_spiral;
        %     l_rew_grt = traj_setup.grad_info.n_rampdown + traj_setup.grad_info.n_M0;
        %     l_rew_dt = round(l_rew_grt*interp_factor);
        %     rewinder_dt = interp1(linspace(0,1,l_rew_grt), rewinder_grt, linspace(0,1,l_rew_dt));
        %
        %     x_m0 = [gradients_nominal(:,1,2);rewinder_dt(:,1)];
        %     y_m0 = [gradients_nominal(:,1,1);rewinder_dt(:,2)];
        %
        %     figure, subplot(1,2,1); plot(x_m0),subplot(1,2,2); plot(y_m0)
        
        interp_factor = 1e-5/dt; %samples/traj_setup.grad_info.n_spiral;
        l_gM0_grt = length(x_m0);
        l_gM0_dt = round(l_gM0_grt*interp_factor);
        x_m0 = interp1(linspace(0,1,l_gM0_grt), x_m0, linspace(0,1,l_gM0_dt));
        y_m0 = interp1(linspace(0,1,l_gM0_grt), y_m0, linspace(0,1,l_gM0_dt));
        
        gradients_M0 =  zeros(l_gM0_dt,gangles,2);
        for i = 1:gangles
            rot = (i-1)*(2*pi/gangles);
            gradients_M0(:,i,1)  = ( x_m0 *cos(rot) + y_m0 *sin(rot));
            gradients_M0(:,i,2)  = (-x_m0 *sin(rot) + y_m0 *cos(rot));
        end
        %     figure, plot(squeeze(gradients_M0(:,1,:)))
        
        gradients_M0r = zeros(size(gradients_M0));
        gradients_M0r(:,:,1)= flipud(gradients_M0(:,:,1));
        gradients_M0r(:,:,2)= flipud(gradients_M0(:,:,2));
        gradients_M0 = [gradients_M0r;gradients_M0];
        %     figure, plot(squeeze(gradients_M0(:,1,:)))
        
        R = [raw_data.head.phase_dir(:,1), raw_data.head.read_dir(:,1), raw_data.head.slice_dir(:,1)  ]; % Rotation matrix the right order..
        
        tRR = 0; % legacy object for sub-dwelltime delay correction
        sR.R = R; % pack rotation matrix and field info in to one object for backwards compatability..
        sR.T = iRD_s.acquisitionSystemInformation.systemFieldStrength_T;
        
        %     figure, colorful_plots(squeeze(gradients_M0(:,1,:)))
        %
        [traj_POET_girf, grad_POET_girf] = apply_GIRF(gradients_M0(:,:,[2 1]) , dt, sR, tRR);
        
        trajectory = traj_POET_girf((l_gM0_dt-samples) + [1:(2*samples2)],:,1:2);         % mod(1-(1:interleaves),interleaves)+1
        gradients_nominal = grad_POET_girf((l_gM0_dt-samples) + [1:(2*samples2)],:,1:2);  % mod(1-(1:interleaves),interleaves)+1
        
        %     figure, colorful_plots(trajectory(:,:,1),trajectory(:,:,2))
        %     figure,
        %     subplot(2,2,1), colorful_plots(squeeze(trajectory(:,1,:))); title('io 1');
        %     subplot(2,2,2), colorful_plots(squeeze([trajectory_nominalr(:,1,:);trajectory_nominal(:,1,:)])); title('nom 1');
        %     subplot(2,2,3), colorful_plots(squeeze(trajectory(:,2,:))); title('io 2');
        %     subplot(2,2,4), colorful_plots(squeeze([trajectory_nominalr(:,2,:);trajectory_nominal(:,2,:)])); title('nom 2');
        %
        %     figure,
        %     subplot(2,2,1), colorful_plots(squeeze(trajectory(:,1,:))); title('io 1');
        %     subplot(2,2,2), colorful_plots(squeeze([trajectory_nominalr(:,1,:);trajectory_nominal(:,1,:)])); title('nom 1');
        %     subplot(2,2,3), colorful_plots(squeeze(trajectory(:,end,:))); title('io 2');
        %     subplot(2,2,4), colorful_plots(squeeze([trajectory_nominalr(:,2,:);trajectory_nominal(:,2,:)])); title('nom 2');
        %
        %
        %     figure, colorful_plots(squeeze(gradients_nominal(:,1,:)))
        %     figure, colorful_plots(squeeze([trajectory_nominalr(:,1,:);trajectory_nominal(:,1,:)]))
    end

%%
    data_store = reshape(data_store, [size(data_store,1)*2 size(data_store,2)/2, channels]);
    
end

% figure, colorful_plots(trajectory(:,:,1), trajectory(:,:,2))

%% Estimate CSM

sample_window = user_opts.ave_for_csm*2*gangles; % both flow sets
if gareps < sample_window
    sample_window = gareps;
    warning(['Using ' num2str(gareps) ' acquisitions for CSM.']);
end

sample_window = 1:sample_window;
user_opts.csm = csm_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);

%% Flow reconstruction cook-book
%% % Real-time reconstruction:
    %% === Frame reconstruction ===
    
    if user_opts.recon_RT_frame
        disp('Real-time reconstruction: Frame-by-frame');
        % RTF and "proper" implementation:
        % reps = gareps;
        % RT_F.MAG = zeros(matrix, matrix, reps, sets);
        % RT_F.PHA = zeros(matrix, matrix, reps);
        %
        % for i = 1:reps
        %      RR_loop_count(i, reps);
        %     sample_window = find(ismember(raw_data.head.idx.repetition, i));
        %     ...
        %     ...
        %     ...
        % end
        
        % Manual GASFib spiral config:
%         gawindow = 5; % acceleration = traj_setup.interleaves/(gawindow/2) - with flow sets;
        gawindow = user_opts.sliding_window_width;
        
        if ~isfield(user_opts,'reps')
            user_opts.reps= (gareps- user_opts.sliding_window_width+1);
        end
        
        reps = user_opts.reps;
        
        RT_F.MAG = zeros(matrix, matrix, reps, sets);
        
%         RT_F.PHA = zeros(matrix, matrix, reps);
        
        for i = 1:reps
%             RR_loop_count(i, reps);
            sample_window = (i-1)*gawindow  + (1:gawindow);
            [RT_F.MAG(:,:,i,:)] = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
        end
        
    end

    %% === Sliding-window reconstruction ===    
    
     if user_opts.recon_RT_SW
        
        disp('Real-time reconstruction: Sliding-window');
        
        if ~isfield(user_opts,'reps')
            user_opts.reps= (gareps- user_opts.sliding_window_width+1);
        end
        
        reps = user_opts.reps;
%         reps = 300;
        RT_SW.MAG = zeros(matrix, matrix, reps, sets);
        
        if sets > 1
        RT_SW.PHA = zeros(matrix, matrix, reps);
        end
        
        for i = 1:reps % 1:gareps
            RR_loop_count(i, reps);
            sample_window = (i*sets - (sets-1)) + (0:(sets*user_opts.sliding_window_width-1));
            
            if sets > 1
                [RT_SW.MAG(:,:,i,:), RT_SW.PHA(:,:,i)] = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
            else
                [RT_SW.MAG(:,:,i)] = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
            end
        end
     end
    
     %% Visualise RT data
     
%      if user_opts.recon_RT_frame || user_opts.recon_RT_SW
%          M_play = [];
%          P_play = [];
%          max_frame = reps;
%          
%          if user_opts.recon_RT_frame && user_opts.recon_RT_SW
%              max_frame= min(size(RT_SW.PHA,3),size(RT_F.PHA,3));
%              if max_frame > 100
%                  max_frame = 100;
%              end
%          end
%          
%          if user_opts.recon_RT_frame
%              M_play = cat(1, M_play, RT_F.MAG(:,:,1:max_frame,1));
%              P_play = cat(1, P_play, RT_F.PHA(:,:,1:max_frame));
%          end
%          if user_opts.recon_RT_SW
%              M_play = cat(1, M_play, RT_SW.MAG(:,:,1:max_frame,1));
%              P_play = cat(1, P_play, RT_SW.PHA(:,:,1:max_frame));
%          end
%          
%          dev.implay_flow(M_play, P_play);
%          
%          clear M_play P_play
%          
%      end

%% % Retrospective-binning reconstruction
% === Limit reconstruction period ===
% sampling_vec = round(linspace(1000,gareps, 2));
% % sampling_vec = gareps; iSV = 1;

    %% === Reference #1: ECG ===
    
    if user_opts.recon_RB_ECG
        disp('Binning reconstruction: ECG');
        %%% Calculate frames from the physiology time stamp
        % === Limit reconstruction period === %         [bin_data, RB_ECG.NominalInterval] = physio_Binning(raw_data.head.physiology_time_stamp(:,1:sampling_vec(iSV)), user_opts.number_cardiac_frames); % figure, plot(raw_data.head.physiology_time_stamp(1,:), 'r-');
        [bin_data, RB_ECG.NominalInterval] = physio_Binning(raw_data.head.physiology_time_stamp, user_opts.number_cardiac_frames); % figure, plot(raw_data.head.physiology_time_stamp(1,:), 'r-');
        
        if user_opts.gridAll == 1
            disp('Gridding all frames together');
            RB_ECG.NominalInterval = 0;
        end
        
        if RB_ECG.NominalInterval == 0
            warning('No ECG trace!'); % and continue
            
            RB_ECG.MAG = zeros([matrix_size 1 sets]);
            if sets > 1
                RB_ECG.PHA = zeros([matrix_size 1]);
            end
            
            sample_window = 1:length(raw_data.data);
            if user_opts.IO == 0
            else
                sample_window = ceil(0.5*sample_window);
            end
            
            if sets > 1
                [RB_ECG.MAG, RB_ECG.PHA] = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
            else
                [RB_ECG.MAG] = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
            end
            
            
        else
            %% Binning Analysis
            RB_ECG.bin_stats = binning_traj_stats(raw_data, bin_data, user_opts);
            
            %% recon
            if user_opts.SNR > 0
                pseudo_reps = 100; disp(['Running ' num2str(pseudo_reps) ' pseudo-reps']);
                
                pseudoRep_phases = user_opts.SNR;
                
                RB_ECG.SNR = zeros( [matrix_size length(pseudoRep_phases)] );
                
                for iPhase = 1:length(pseudoRep_phases)
                    RR_loop_count(iPhase,length(pseudoRep_phases));
                    sample_window = bin_data{ pseudoRep_phases(iPhase) };
                    
                    %                 disp('--normal loop--');
                    %                 tic
                    %                 img_pr = zeros([matrix_size pseudo_reps sets]);
                    %                 %
                    %                 for i = 1:pseudo_reps
                    % %                     RR_loop_count(i,pseudo_reps);
                    %                     img_pr(:,:,i,:) = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
                    %                 end
                    %                 toc
                    %
                    
                    %                 img_pr = abs(img_pr(:,:,:,1));
                    %
                    %                 g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
                    %                 g(g < eps) = 1;
                    %                 RB_ECG.SNR = mean(img_pr,3)./g;
                    % %                 figure, imagesc(RB_ECG.SNR,[0 50]);
                    %                 user_opts.snr_1=RB_ECG.SNR ;
                    
                    img_pr2= cell(1,pseudo_reps);
                    nhlbi_toolbox.parpool_setup(24);
                    
                    parfor i = 1:pseudo_reps
                        img_pr2{i} = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
                        
                    end
                    
                    img_pr = zeros([matrix_size pseudo_reps sets]);
                    for i = 1:pseudo_reps
                        img_pr(:,:,i,:) = img_pr2{i};
                    end
                    %                 toc
                    
                    img_pr = abs(img_pr(:,:,:,1));
                    
                    g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
                    g(g < eps) = 1;
                    RB_ECG.SNR(:,:,iPhase) = mean(img_pr,3)./g;
                    
                    %                 figure, imagesc(RB_ECG.SNR,[0 50]);
                    %                 user_opts.snr_2=RB_ECG.SNR ;
                end
            end
            
            % Standard recon
            RB_ECG.MAG = zeros([matrix_size user_opts.number_cardiac_frames sets]);
            if sets > 1
                RB_ECG.PHA = zeros([matrix_size user_opts.number_cardiac_frames]);
            end
            
            for i = 1:length(bin_data)
%                  RR_loop_count(i, length(bin_data));
                
                if user_opts.IO == 0
                    sample_window = bin_data{i};
                else
                    sample_window = ceil(0.5*bin_data{i});
                end
                
                if sets > 1
                    [RB_ECG.MAG(:,:,i,:), RB_ECG.PHA(:,:,i)] = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
                else
                    [RB_ECG.MAG(:,:,i)] = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
                end
            end
            
            
        end
        
        %% Flow Analysis
%         if sets > 1 && user_opts.SNR == 0
%             [im_size] = [size(RB_ECG.MAG,1) size(RB_ECG.MAG,2)];
%             ao = circ_roi_RR(10,round(im_size/2), im_size); bg = [];
%             
%             % ao = auto_flow_ao_beta(RB_ECG.MAG(:,:,:,1), RB_ECG.PHA);
%             % [~,~,bg] = auto_GREflow_ao_beta(RB_ECG.MAG(:,:,:,1), RR_pi2dicom(RB_ECG.PHA));
%             
%             info.PixelSpacing = [iRD_s.encoding.reconSpace.fieldOfView_mm.x iRD_s.encoding.reconSpace.fieldOfView_mm.y]./iRD_s.encoding.reconSpace.matrixSize.x;
%             info.CardiacNumberOfImages = user_opts.number_cardiac_frames;
%             info.NominalInterval = RB_ECG.NominalInterval;
%             
%             [RB_ECG.flow_data] = GRE_flow_RR2(RB_ECG.PHA, ao, bg, info, venc);
%             RB_ECG.roi_ao = ao; RB_ECG.roi_bg = bg;
%         end
        
    end
    
    %% === Reference #2: DC signal-based self-gating  ===
    
    if user_opts.recon_RB_DCSG
        
%         resp_frames = 8; 
        resp_frames = 1;
        [dcsg] = self_gate_DC(data_store, iRD_s, [user_opts.number_cardiac_frames resp_frames]);
        bin_data = cell(1,user_opts.number_cardiac_frames);
        for i = 1:user_opts.number_cardiac_frames
            bin_data{i} = find(dcsg.ecg_binned == i);
        end
        
%         bin_data = cell(1,resp_frames);
%         for i = 1:resp_frames
%             bin_data{i} = find(dcsg.ecg_binned == i);
%         end
%         
%         RB_DCSG.MAG = zeros([matrix_size resp_frames ]);
%         
%         for i = 1:resp_frames
%             RR_loop_count(i, length(bin_data));
%             sample_window = bin_data{i};
%             [RB_DCSG.MAG(:,:,i,:)] = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
%         end
%         
%         implay_RR(RB_DCSG.MAG(:,:,2:end))
%         
        
        %% Binning Analysis
        RB_DCSG.bin_stats = binning_traj_stats(raw_data, bin_data, user_opts);
            
        %% recon
        
        RB_DCSG.MAG = zeros([matrix_size user_opts.number_cardiac_frames sets]);
        RB_DCSG.PHA = zeros([matrix_size user_opts.number_cardiac_frames]);
        
        for i = 1:user_opts.number_cardiac_frames
            RR_loop_count(i, length(bin_data));
            sample_window = bin_data{i};
            [RB_DCSG.MAG(:,:,i,:), RB_DCSG.PHA(:,:,i)] = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
        end
        
        %% Flow Analysis
        
        ao = auto_flow_ao_beta(RB_DCSG.MAG(:,:,:,1), RB_DCSG.PHA); close;
        [~,~,bg] = auto_GREflow_ao_beta(RB_DCSG.MAG(:,:,:,1), RR_pi2dicom(RB_DCSG.PHA)); close;
        
        info.PixelSpacing = [iRD_s.encoding.reconSpace.fieldOfView_mm.x iRD_s.encoding.reconSpace.fieldOfView_mm.y]./iRD_s.encoding.reconSpace.matrixSize.x;
        info.CardiacNumberOfImages = user_opts.number_cardiac_frames;
        info.NominalInterval = 1000./dcsg.ecg_freq; % ms;
        
        [RB_DCSG.flow_data] = GRE_flow_RR2(RB_DCSG.PHA, ao, bg, info, venc);
        RB_DCSG.roi_ao = ao; RB_DCSG.roi_bg = bg;
        RB_DCSG.dcsg = dcsg;
        
    end

    %% === Flow-based self-gating  ===

    if user_opts.recon_RB_FSG
        
        % === approximate ROI to generate flow-ish curves:
        
        dA = (iRD_s.encoding.reconSpace.fieldOfView_mm.x/10)./iRD_s.encoding.reconSpace.matrixSize.x; dA = dA.^2;
        temp_roi = circ_roi_RR(round(2.0/sqrt(dA)),matrix_size.*0.5, matrix_size); % radius ~ 2cm
%         figure, imshow(mean(RT_SW.PHA,3),[]); hold on; contour(temp_roi,'g');
        [~, temp_sF] = dicomFlowExtract_RR(RT_SW.PHA, temp_roi, dA, venc); close;close;
        
        % peak-to-peak flow detection
        [~,heart_rate,~,beat_vec] = beat2beat_CO_RR2(temp_sF,iRD_s.sequenceParameters.TR*2, 1); close all;  % figure, plot(COdata.sF, 'k-');
        rtsg.heart_rate = heart_rate; rtsg.rt_flow = temp_sF;
        
        % Sliding window indexing (this codecan be majorly compressed..)
        % ===
        temp_d2=cell(1,user_opts.number_cardiac_frames);
        for i = 1:length(beat_vec)
            
            rep_ind = beat_vec{i};
            cardiac_ind_beat = rep_ind -rep_ind(1);
            ind_d1 = floor((1:length(rep_ind))/((length(rep_ind))/(user_opts.number_cardiac_frames-0.01)));
            
            for i = 0:(user_opts.number_cardiac_frames-1)
                temp_d2{i+1} = [temp_d2{i+1} rep_ind(1) + find((ind_d1 == i))-1]; % pad for reps
            end
        end
        
        bin_data = cell(1,user_opts.number_cardiac_frames);
        for i = 1:user_opts.number_cardiac_frames
            frames = temp_d2{i};
            temp_d3 = [];
            for j = 1:length(frames)
       
                SW_frames = frames(j) + (0:user_opts.sliding_window_width-1);
                
                % choose middle of the frames:
                window_acc = 2; rec_win = (1:window_acc) - round(window_acc/2); % <== HARD CODE
                temp_d3 = [temp_d3 SW_frames(floor(length(SW_frames)/2)+ rec_win)];
                
            end
            frames = temp_d3;
            frames = find(ismember(1 + raw_data.head.idx.repetition, frames)); % this line is reliant on the repetition number changing for each spiral arm (every two arms with flow)
            bin_data{i} = frames;
        end
        % ===
        
        %% Binning Analysis
        RB_FSG.bin_stats = binning_traj_stats(raw_data, bin_data, user_opts);
            
        %% recon
                 
        RB_FSG.MAG = zeros([matrix_size user_opts.number_cardiac_frames sets]);
        RB_FSG.PHA = zeros([matrix_size user_opts.number_cardiac_frames]);
        
        for i = 1:user_opts.number_cardiac_frames
            RR_loop_count(i, length(bin_data));
            sample_window = bin_data{i};
            [RB_FSG.MAG(:,:,i,:), RB_FSG.PHA(:,:,i)] = sample_window_recon(raw_data, data_store, trajectory, gradients_nominal, sample_window, user_opts);
        end
        
        %% Flow Analysis
        
        ao = auto_flow_ao_beta(RB_FSG.MAG(:,:,:,1), RB_FSG.PHA);
        [~,~,bg] = auto_GREflow_ao_beta(RB_FSG.MAG(:,:,:,1), RR_pi2dicom(RB_FSG.PHA));
        
        info.PixelSpacing = [iRD_s.encoding.reconSpace.fieldOfView_mm.x iRD_s.encoding.reconSpace.fieldOfView_mm.y]./iRD_s.encoding.reconSpace.matrixSize.x;
        info.CardiacNumberOfImages = user_opts.number_cardiac_frames;
        info.NominalInterval = 60000/heart_rate;
        
        [RB_FSG.flow_data] = GRE_flow_RR2(RB_FSG.PHA, ao, bg, info, venc);
        RB_FSG.roi_ao = ao; RB_FSG.roi_bg = bg;
        RB_FSG.rtsg = rtsg;
        
    end    
    
    %% Visualise binned data
if sets > 2
    if user_opts.recon_RB_ECG || user_opts.recon_RB_DCSG || user_opts.recon_RB_FSG
       
        M_play = [];
        P_play = [];
        
        if user_opts.recon_RB_ECG
            M_play = cat(1, M_play, RB_ECG.MAG(:,:,:,1));
            P_play = cat(1, P_play, RB_ECG.PHA);
        end
        if user_opts.recon_RB_DCSG
            M_play = cat(1, M_play, RB_DCSG.MAG(:,:,:,1));
            P_play = cat(1, P_play, RB_DCSG.PHA);
        end
        if user_opts.recon_RB_FSG
            M_play = cat(1, M_play, RB_FSG.MAG(:,:,:,1));
            P_play = cat(1, P_play, RB_FSG.PHA);
        end
        
        dev.implay_flow(M_play, P_play);
        
        clear M_play P_play
    end
end

%% Export all

RTFM_output.user_opts = user_opts;

if user_opts.recon_RT_frame
    RTFM_output.RT_F = RT_F;
end
if user_opts.recon_RT_SW
    RTFM_output.RT_SW = RT_SW;
end
if user_opts.recon_RB_ECG
    RTFM_output.RB_ECG = RB_ECG;
end
if user_opts.recon_RB_DCSG
    RTFM_output.RB_DCSG = RB_DCSG;
end
if user_opts.recon_RB_FSG
    RTFM_output.RB_FSG = RB_FSG;
end

end

function [csm] = csm_recon(raw_data, data_store, traj, grad, indices, user_opts)

gangles = 1 + double(max(raw_data.head.idx.kspace_encode_step_1));
if user_opts.IO == 0
    xi = raw_data.head.idx.kspace_encode_step_1(indices)+1;
else
    xi = raw_data.head.idx.kspace_encode_step_1(indices*2)+1;
end

[samples,~,channels] = size(data_store);
data = data_store(:,indices,:);
    
    % ### >>> comment out to enable FULL GRID (slower) ###
    data2 = zeros(samples, gangles, channels);
    for j = 1:gangles
        nnn = length(find(xi==j));
        xitemp(j) = nnn;
        data2(:,j,:) = mean(data(:,xi==j,:),2);
    end
    xi = find(xitemp);
    data = data2(:,xi,:);
    % ### <<< comment out to enable FULL GRID (slower) ###
    
    data = reshape(data,length(xi)*samples,channels);

    % Set-up plan
    traj = traj(:,xi,:)*(pi/(max(traj(:))));
    traj = reshape(traj,length(xi)*size(traj,1),2);    
    st = nufft_init(traj, user_opts.matrix_size, [6 6],user_opts.matrix_size.*2, user_opts.matrix_size./2);
    
    %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
    grad = reshape(grad(:,xi,:),length(xi)*size(grad,1),2);
    grad = complex(grad(:,1),grad(:,2));
    kk = complex(traj(:,1),traj(:,2));
    weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; 
    
    x = nufft_adj(data.*repmat(weights,[1,channels]), st);
    img_coil = squeeze(x);
    
    csm  = ismrm_estimate_csm_mckenzie(squeeze(img_coil)); % montage_RR(abs(csm));

end

function [mag, pha] = sample_window_recon(raw_data, data_store, traj, grad, indices, user_opts)

    sets = 1 + double(max(raw_data.head.idx.set));
    gangles = 1 + double(max(raw_data.head.idx.kspace_encode_step_1));
    
    if user_opts.IO == 0
        traj_list = raw_data.head.idx.kspace_encode_step_1(indices)+1;
    else
        traj_list = raw_data.head.idx.kspace_encode_step_1(2*indices)+1;
    end

    x = zeros([user_opts.matrix_size sets]);
    [samples,~,channels] = size(data_store);
    
    for iSet = 1:sets
        if user_opts.IO == 0
            si = find(1+raw_data.head.idx.set(indices) == iSet);
        else
            si = find(1+raw_data.head.idx.set(2*indices) == iSet);
        end
        
        xi = traj_list(si);
        data = data_store(:,indices(si),:);
         
        if user_opts.SNR > 0
            data = data + complex(randn(size(data)),randn(size(data)));
        end
        
        % ### >>> comment out to enable FULL GRID (slower) ###
        data2 = zeros(samples, gangles, channels);
        for j = 1:gangles
            nnn = length(find(xi==j));
            xitemp(j) = nnn;
            data2(:,j,:) = mean(data(:,xi==j,:),2);
        end
        xi = find(xitemp);
        data = data2(:,xi,:);
        % ### <<< comment out to enable FULL GRID (slower) ###
        
        data = reshape(data,length(xi)*samples,channels);
        
        traj2 = traj(:,xi,:)*(pi/(max(traj(:))));
        traj2 = reshape(traj2,length(xi)*size(traj,1),2);
            
        grad2 = reshape(grad(:,xi,:),length(xi)*size(grad,1),2);
        
       
        
        [x(:,:,iSet)] = recon_cg_RR(data, traj2, grad2,  user_opts);
        
    end
%     mag = abs(x);
    mag = x;
    if sets > 1
        pha = atan2(imag(x(:,:,2).*conj(x(:,:,1))), real(x(:,:,2).*conj(x(:,:,1))));
    end
end

function [complexImage] = recon_cg_RR(data, traj, grad, user_opts)
%         [samples,~,channels] = size(data);
%         data = reshape(data,length(indices)*samples,channels);
%         
        grad = complex(grad(:,1),grad(:,2));
        kk = complex(traj(:,1),traj(:,2));
        weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
        
        st = nufft_init(traj, user_opts.matrix_size, [6 6],user_opts.matrix_size.*2, user_opts.matrix_size./2);
        
        % Custom SENSE implementation 
        I = ones(user_opts.matrix_size(1)); D = repmat(weights, [1 size(user_opts.csm,3)]);
       
        complexImage = cg_RR(data(:), st, I, D, user_opts.csm, weights, user_opts.iter);
end

function [ccc] = binning_traj_stats(raw_data, bin_data, user_opts)

%% % Output stats: Count number of arms per frame
debug = 0;
gangles = 1 + double(max(raw_data.head.idx.kspace_encode_step_1));

if debug
    plot_bin_fig = figure;
    plot_bin_tabgp = uitabgroup(plot_bin_fig); clear a tab;
end


for i = 1:user_opts.number_cardiac_frames
    a = bin_data{i};
    t_arms = raw_data.head.idx.kspace_encode_step_1(a)+1;
    s_ind = find(raw_data.head.idx.set(a) == 0);
    t_arms0 = t_arms(s_ind);
    
    if debug
        bbb(1,i) = length(unique( t_arms0 ));
        s_ind = find(raw_data.head.idx.set(a) == 1);
        t_arms1 = t_arms(s_ind);
        bbb(1,2) = length(unique( t_arms1 ));
    end
    
    [counts] = hist(t_arms0, (1:gangles)); % 0.5 + to form bins?
    ccc(:,i) = counts;
end

if debug
    for iSV = 1:user_opts.number_cardiac_frames
    temp_str = ['Min #unique arms: ' num2str(min(bbb(:))) ' Gangle-x: ' num2str(gangles/min(bbb(:))) ]; % interleaves/min(bbb(:))
    plot_bin_tab{i} = uitab(plot_bin_tabgp,'Title',['iSV ' num2str(iSV)]);
    plot_bin_axes{i} = axes('parent', plot_bin_tab{i});
    imagesc(ccc, 'Parent',plot_bin_axes{i});
    axis image; xlabel('Frames'); ylabel('Gangles'); title(temp_str); shg;
    end
end

end