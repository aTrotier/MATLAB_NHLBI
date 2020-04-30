function [mag_RT] = recon_spiral_SE_beta(dfile,  nfile, user_opts)
% [mag_RT, seheader, complex_out] = recon_spiral_SE_beta(dfile,  nfile, SpiDes)
% Spin-echo sequence reconstruction (uses user_float field to determine TSE
% factor, FOV, matrix etc).
% R Ramasawmy NHLBI Aug 2018

make_nhlbi_toolbox;
make_dev;

addpath(nhlbi_toolbox.run_path_on_sys('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Ramasawmy\local_MATLAB\ismrm_sunrise_matlab-master\irt\mex\v7'));
addpath(nhlbi_toolbox.run_path_on_sys('\\hl-share.nhlbi.nih.gov\tmb\Lab-Campbell\Ramasawmy\local_MATLAB\ismrm_sunrise_matlab-master\irt\mri'));

dfile = nhlbi_toolbox.run_path_on_sys(dfile);
nfile = nhlbi_toolbox.run_path_on_sys(nfile);

%% Grab XML header
iRD_s = nhlbi_toolbox.h5read_xml(dfile);
mag_RT.mrd_header = iRD_s;

%% Load data
raw_data = h5read(dfile,'/dataset/data');
iRD_s = read_h5_header(dfile);

disp(['Reconstructing: ' iRD_s.measurementInformation.protocolName]);

%% user_opts default:
% need to configure to avoid overwriting input selection
if exist('user_opts','var')
    user_opts_in = user_opts;
else
    user_opts_in = struct();
end
user_opts_in_list = fieldnames(user_opts_in);

% default parameters
user_opts.delayFactor   = 0;
user_opts.vds           = 100;
user_opts.iter          = 3;
user_opts.IO            = 0;
user_opts.SNR           = 0;
user_opts.import_traj   = [];
user_opts.field_map     = 0;
user_opts.sep_tse       = 0;
user_opts.TSEf          = 1;
user_opts.ramp_samples  = 0;
user_opts.gmax          = 2.4;       % traj_setup.gMax = iRD_s.encoding.userParameterDouble.value;
user_opts.smax          = 14414.4;   % traj_setup.sMax = iRD_s.encoding.userParameterDouble_1.value;

% === update user_opts with input data === 
for i = 1:length(user_opts_in_list)
    if isfield(user_opts, user_opts_in_list{i})
        user_opts.(matlab.lang.makeValidName(user_opts_in_list{i})) = user_opts_in.(matlab.lang.makeValidName(user_opts_in_list{i}));
    end
end

disp(' ');disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
disp('User options:');
disp(user_opts);  disp(' ');
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');

%%
% nhlbi_toolbox.plot_experiment(raw_data);

samples     = double(raw_data.head.number_of_samples(1));
interleaves = double(max(raw_data.head.idx.kspace_encode_step_1))+1;
pe2         = 1+double(max(raw_data.head.idx.kspace_encode_step_2));
averages    = (1 + double(max(raw_data.head.idx.average)));
slices      = (1 + double(max(raw_data.head.idx.slice)));
contrasts   = (1 + double(max(raw_data.head.idx.contrast)));
phases      = (1 + double(max(raw_data.head.idx.phase)));
reps        = (1 + double(max(raw_data.head.idx.repetition)));
sets        = (1 + double(max(raw_data.head.idx.set)));
channels    = double(raw_data.head.active_channels(1));

matrix = iRD_s.encoding.reconSpace.matrixSize.x;
matrix_size = [matrix matrix];

dt = raw_data.head.sample_time_us(1)*1e-6;
FOV = iRD_s.encoding.reconSpace.fieldOfView_mm.x/10;

temp = raw_data.head.user_float(:,2);
user_opts.TSEf = double(temp(8)); if user_opts.TSEf==0; user_opts.TSEf = 1; end; 
% user_opts.vds = double(temp(6));
clear temp;

TSE_acq_index = 1:length(raw_data.data);
if user_opts.field_map == 1
    samples     = double(raw_data.head.number_of_samples(2));
    samples_map = double(raw_data.head.number_of_samples(1));
       
%     FM_acq_index  = 1:(user_opts.TSEf+1):length(raw_data.data);
    FM_acq_index  = find(raw_data.head.number_of_samples==raw_data.head.number_of_samples(1));
    TSE_acq_index(FM_acq_index) = [];
end

% user_opts.vds = 100;%temp(6);

disp(' ');disp('### Experiment Dimensions ###');disp(' ');
Experiment_parameters = {'Samples','Interleaves','Matrix', 'PE2', 'Averages', 'Slices', 'Contrasts', 'Phases', 'Reps', 'Sets', 'Channels'}';
Value = [samples interleaves iRD_s.encoding.encodedSpace.matrixSize.x pe2 averages slices contrasts phases reps sets channels]';
disp(table( Experiment_parameters,Value )); clear Experiment_parameters Value; disp(' ');

seheader.samples = samples;
seheader.dt = dt;
seheader.number_aqs = length(raw_data.data);
seheader.averages = averages;
seheader.channels = channels;

mag_RT.mrd_header.seheader = seheader; 

%% Noise Adjust
% iRD_s_noise = nhlbi_toolbox.h5read_xml(nfile);
% nhlbi_toolbox.check_noise_dependency(iRD_s, iRD_s_noise)
dmtx = nhlbi_toolbox.noise_adjust(nfile, iRD_s, raw_data.head.sample_time_us(1)*1e-6, nhlbi_toolbox);

%% Trajectory

if isempty(user_opts.import_traj)
    %% Build Nominal Fully Sampled traj and gradients
    
    disp(['FOV: ' num2str(FOV)])
    disp(['TSE factor: ' num2str(user_opts.TSEf)])
    FOV = [FOV -1*FOV*(1 - user_opts.vds/100)];
    
    smax = 14414.4; %
    % smax = 3669.72;
    krmax = 1/(2*(FOV(1)/matrix_size(1)));
    [k,g] = vds(user_opts.smax,user_opts.gmax, dt, interleaves, FOV, krmax); close;
    
    %%
    
    % poet_data = POET_read_textfile('C:\Users\ramasawmyr\Desktop\temp\tses_s40_m320_fov350_vds100_bw980.txt');
    
    %% Rotate Fibonacci steps
    if samples/2 > length(k) % single ADC for IO
        % if samples > length(k) % normal ADC for out
        %     samples2 = length(k);
        samples2 = length(k);
    else
        samples2 = samples/2;
    end
    trajectory_nominal = zeros(samples2,interleaves,2);
    gradients_nominal =  zeros(samples2,interleaves,2);
    
    neg = -1;
    for solid_int= 1:interleaves
        rot = (solid_int-1)*(2*pi/interleaves);
        trajectory_nominal(:,solid_int,1) = neg*-( real(k(1:samples2)) *cos(rot) + imag(k(1:samples2)) *sin(rot));
        trajectory_nominal(:,solid_int,2) = neg*-(-real(k(1:samples2)) *sin(rot) + imag(k(1:samples2)) *cos(rot));
        gradients_nominal(:,solid_int,1)  = neg*-( real(g(1:samples2)) *cos(rot) + imag(g(1:samples2)) *sin(rot));
        gradients_nominal(:,solid_int,2)  = neg*-(-real(g(1:samples2)) *sin(rot) + imag(g(1:samples2)) *cos(rot));
    end
    
    trajectory_nominalr = zeros(samples2, interleaves, 2);
    gradients_nominalr = zeros(samples2, interleaves, 2);
    for solid_int= 1:interleaves
        rot = (solid_int-1)*(2*pi/interleaves);
        trajectory_nominalr(:,solid_int,1) = neg*-( real(fliplr(k(1:samples2))) *cos(rot) + imag(fliplr(k(1:samples2))) *sin(rot));
        trajectory_nominalr(:,solid_int,2) = neg*-(-real(fliplr(k(1:samples2))) *sin(rot) + imag(fliplr(k(1:samples2))) *cos(rot));
        gradients_nominalr(:,solid_int,1)  = neg*-( real(fliplr(g(1:samples2))) *cos(rot) + imag(fliplr(g(1:samples2))) *sin(rot));
        gradients_nominalr(:,solid_int,2)  = neg*-(-real(fliplr(g(1:samples2))) *sin(rot) + imag(fliplr(g(1:samples2))) *cos(rot));
    end
    
    % unified trajectory
    
    trajectory_io = zeros(samples2*2, interleaves, 2);
    gradients_io = zeros(samples2*2, interleaves, 2);
    
    for i = 1:interleaves
        xi = mod(i-1+(interleaves/2),interleaves)+1; xo = i;
        trajectory_io(:,i,:) = [trajectory_nominalr(:,xi,:);trajectory_nominal(:,xo,:)];
        trajectory_io(:,i,:) = [trajectory_nominalr(:,xi,:);trajectory_nominal(:,xo,:)];
        gradients_io(:,i,:) = [gradients_nominalr(:,xo,:);gradients_nominal(:,xo,:)];
        gradients_io(:,i,:) = [gradients_nominalr(:,xo,:);gradients_nominal(:,xo,:)];
    end
    
    trajectory_nominal = trajectory_io;
    gradients_M0 = gradients_io;
   
    
else
    load(user_opts.import_traj);
    samples2 = 0.5*samples;
    
%     R = [raw_data.head.phase_dir(:,1), raw_data.head.read_dir(:,1), raw_data.head.slice_dir(:,1)  ]; % Rotation matrix the right order..
    R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; % Rotation matrix the right order..
    
    tRR = 0; % legacy object for sub-dwelltime delay correction
    sR.R = R; % pack rotation matrix and field info in to one object for backwards compatability..
    sR.T = iRD_s.acquisitionSystemInformation.systemFieldStrength_T;
    
    % v1
    [traj_POET_girf] = apply_GIRF(gradients_M0(:,:,[1 2]) , 1e-6, sR, tRR);
    
    if user_opts.ramp_samples == 1
        adc_start = 0;
        trajectory_nominal = traj_POET_girf( round(linspace(1,size(traj_POET_girf,1),samples)) ,:,1:2);
        
        FOV = [FOV -1*FOV*(1 - user_opts.vds/100)];
        smax = 14414.4; %
        krmax = 1/(2*(FOV(1)/matrix_size(1)));
        [k,g] = vds(user_opts.smax,user_opts.gmax, dt, interleaves, FOV, krmax); close;
        samples2 = length(k);
        
        trajectory_nominal= trajectory_nominal((samples/2 + [-(samples2-1):(samples2)]),:,:);
    
    else
        
        trajectory_nominal = traj_POET_girf( (adc_start+2 + round(dt*1e6*0.5)) + round([0:(samples-1)]*(dt*1e6)) ,:,1:2);
    end
    
    % v2. Resample gradients first. (looks better)
        %     gradients_M0_dt = gradients_M0(round(1:(dt/1e-6):end) ,:,:);
        %     [traj_POET_girf] = apply_GIRF(gradients_M0_dt(:,:,[1 2]) , dt, sR, tRR);
        % %     figure,colorful_plots(traj_POET_girf(:,:,1), traj_POET_girf(:,:,2))
        %     
        %     trajectory_nominal = traj_POET_girf( (round(adc_start*(1e-6/dt))+1) + round([0:(samples-1)]) ,:,1:2); 
        % % %     trajectory_nominal = traj_POET_girf( (round(adc_start*(1e-6/dt))+3) + round([0:(samples-1)]) ,:,1:2); 
        % %     figure,colorful_plots(trajectory_nominal(:,:,1), trajectory_nominal(:,:,2));
end

%% Field map <option>


if user_opts.field_map == 1
    %%
    
    FOV = [FOV(1) -1*FOV(1)*(1 - user_opts.vds/100)];
    
    smax = 14414.4; %
    % smax = 3669.72;
    krmax = 1/(2*(FOV(1)/matrix_size(1)));
    [k,g] = vds(user_opts.smax,user_opts.gmax, dt, interleaves, FOV, krmax); close;
    
    if samples_map > length(k) % 
        samples2_map = length(k);
    else
        samples2_map = samples_map;
    end
    trajectory_FM = zeros(samples2_map,interleaves,2);
    gradients_FM =  zeros(samples2_map,interleaves,2);
    
    neg = 1;
    for solid_int= 1:interleaves
        rot = (solid_int-1)*(2*pi/interleaves);
        trajectory_FM(:,solid_int,1) = neg*-( real(k(1:samples2_map)) *cos(rot) + imag(k(1:samples2_map)) *sin(rot));
        trajectory_FM(:,solid_int,2) = neg*-(-real(k(1:samples2_map)) *sin(rot) + imag(k(1:samples2_map)) *cos(rot));
        gradients_FM(:,solid_int,1)  = neg*-( real(g(1:samples2_map)) *cos(rot) + imag(g(1:samples2_map)) *sin(rot));
        gradients_FM(:,solid_int,2)  = neg*-(-real(g(1:samples2_map)) *sin(rot) + imag(g(1:samples2_map)) *cos(rot));
    end
    
    % =======================
    % apply GIRF
    % =======================
    
    
    R = [ raw_data.head.read_dir(:,1),raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix
    tRR = 0;
    sR.R = R;
    sR.T = iRD_s.acquisitionSystemInformation.systemFieldStrength_T;
    
    trajectory_FM = apply_GIRF(gradients_FM, dt, sR, tRR );
    trajectory_FM = trajectory_FM(:,:,1:2);
    
    %%
    
    kspace_FM = complex(zeros([samples_map interleaves pe2 1 slices 1 phases reps 2 channels],'single'));
    % size(kspace_FM)
    for i = 1:length(FM_acq_index)
        ii = FM_acq_index(i); 
        d1 = raw_data.data{ii};
        d2 = complex(d1(1:2:end), d1(2:2:end));
        d3 = reshape(d2, samples_map, channels); %  RE & IM (2)
        
%         if nargin > 1
            d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
%         end
        
        kspace_FM(:,...
            raw_data.head.idx.kspace_encode_step_1(ii)+1, ...
            raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
            1, ....
            raw_data.head.idx.slice(ii)+1 , ...
            1, ...
            raw_data.head.idx.phase(ii)+1, ...
            raw_data.head.idx.repetition(ii)+1, ...
            raw_data.head.idx.set(ii), ...
            :) = d3;
    end
    
    kspace_FM = squeeze(kspace_FM);  size(kspace_FM)
    if ndims(kspace_FM) == 4 % single slice
       kspace_FM = reshape(kspace_FM, [size(kspace_FM,1) size(kspace_FM,2) 1 size(kspace_FM,3) size(kspace_FM,4)]); 
    end
    
    % >> PCA
%     pca_channels = 8;
%     channels_orig = channels;
%     channels = pca_channels;
%     kspace_FM = nhlbi_toolbox.coil_pca(kspace_FM, pca_channels);
%     
    % 
    
    %%
    fov = FOV(1)*10; % mm
    mask = true(matrix_size(1),matrix_size(2));
    sizeMask = size(mask);
    nufft_args = {sizeMask, [6 6], 2*sizeMask, sizeMask/2, 'table', 2^12, 'minmax:kb'};

    omega = trajectory_FM*(pi/max(max(max(trajectory_FM))));
    omega = reshape(omega,interleaves*samples2_map,2);
    csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
    
    kx = (omega(:,1)./pi).*matrix_size(1)/2;
    ky = (omega(:,2)./pi).*matrix_size(2)/2;
    ksp = [kx(:) ky(:)]./fov;
    G = Gmri(ksp, mask, 'fov', fov, 'basis', {'dirac'}, 'nufft', nufft_args);
    csm_weights = abs(mri_density_comp(ksp, 'pipe','G',G.arg.Gnufft)); % figure, plot(csm_weights(1:samples2_map))

    csm = zeros([matrix_size channels slices]);
    for iSlice= 1:slices
        
        data_temp = zeros([samples_map interleaves channels],'single');
        data_temp(:,1:2:end,:) = kspace_FM(:,1:2:end,iSlice,1,:);
        data_temp(:,2:2:end,:) = kspace_FM(:,2:2:end,iSlice, 2,:);
        data_temp = reshape(data_temp, [samples_map*interleaves channels]);
        
        % -- filter --
        hann = @(n) 0.5*(1 - cos(2*pi*[0:n-1]/n));
        fwin = hann(2*samples_map); 
        fwin = fwin(samples_map+1:end);
% %         figure, plot(fwin)
%         blackman = @(n) 0.42 - 0.5*cos(2*pi*([0:(0.5*n-1)])./(0.5*n-1))+ 0.08*cos(4*pi*([0:(0.5*n-1)])./(0.5*n-1));
%         fwin = blackman(4*samples_map); 
%         fwin = fwin(samples_map+1:end);
% %         figure, plot(fwin)
        data_temp = data_temp.*repmat(reshape(fwin,[samples_map 1]), [interleaves channels]);
        
        x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st)/numel(repmat(csm_weights,[1,channels]));
        csm(:,:,:,iSlice) = ismrm_estimate_csm_walsh( x );
        ccm_roemer_optimal = ismrm_compute_ccm(csm(:,:,:,iSlice), eye(channels)); % with pre-whitened
        
        %     img_tse(:,:,iTSE) = sqrt(sum(x.*conj(x),3));
        img_field_map(:,:,iSlice) = abs( sum( squeeze( x ) .* ccm_roemer_optimal, 3) );
    end
    
    montage_RR(img_field_map);
    
    mag_RT.img_field_map       = img_field_map;
    
    %% FIELD MAP ESTIMATION
    
    for i = 1:2
        
        % omega = trajectory_FM(:,double(raw_data.head.idx.kspace_encode_step_1(FM_acq_index(i:2:end)))+1,:)*(pi/max(max(max(trajectory_FM))));
        omega = trajectory_FM(:,i:2:end,:)*(pi/max(max(max(trajectory_FM))));
        omega = reshape(omega,interleaves/2*samples2_map,2);
        csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
        
        kx = (omega(:,1)./pi).*matrix_size(1)/2;
        ky = (omega(:,2)./pi).*matrix_size(2)/2;
        ksp = [kx(:) ky(:)]./fov;
        G = Gmri(ksp, mask, 'fov', fov, 'basis', {'dirac'}, 'nufft', nufft_args);
        csm_weights = abs(mri_density_comp(ksp, 'pipe','G',G.arg.Gnufft)); % figure, plot(csm_weights(1:samples2_map))
        
        for iSlice= 1:slices
            %         data_temp = zeros([samples_map interleaves/2 channels],'single');
            data_temp = squeeze(kspace_FM(:,(i:2:end),iSlice,i,:));
            data_temp = data_temp.*repmat(reshape(fwin,[samples_map 1]), [1 interleaves/2 channels]);
      
            data_temp = reshape(data_temp, [samples_map*interleaves/2 channels]);
            
            I = ones(matrix_size(1)); D = repmat(csm_weights, [1 channels]);
            
%             complexImage(:,:,i,iSlice) = cg_RR(data_temp(:), csm_st, I, D, csm(:,:,:,iSlice), csm_weights, 5);
            complexImage(:,:,i,iSlice) = cg_RR(data_temp(:), csm_st, I, D, csm_tse(:,:,:,iSlice,1), csm_weights, 5);
            
        end
    end
    
    TE_diff = 0.001; % s
 
    for iSlice= 1:slices
        pha(:,:,iSlice) = atan2(imag(complexImage(:,:,2,iSlice).*conj(complexImage(:,:,1,iSlice))), real(complexImage(:,:,2,iSlice).*conj(complexImage(:,:,1,iSlice))));
        B0_mask(:,:,iSlice) = abs( img_field_map(:,:,iSlice)) > mean(abs(img_field_map(:))) ;
        B0_mask(:,:,iSlice) = bwareaopen(B0_mask(:,:,iSlice),30);
        
        B0_map(:,:,iSlice)= pha(:,:,iSlice) ./(2*pi*TE_diff);
        B0_map_filt(:,:,iSlice) = imguidedfilter(B0_map(:,:,iSlice).*B0_mask(:,:,iSlice), 'DegreeOfSmoothing',150);
    end
    
    montage_RR(abs(complexImage));
    
      montage_RR([B0_map imguidedfilter(B0_map.*B0_mask, 'DegreeOfSmoothing',150)],[-80 80]); colormap('parula');
%     hold on, contour([B0_mask B0_mask],'w');
    
    mag_RT.B0_map       = B0_map;
    mag_RT.B0_map_filt  = B0_map_filt;
    mag_RT.B0_mask      = B0_mask;
    mag_RT.B0_img       = img_field_map;

end

%% RECON MEAN ("brute force")

kspace = complex(zeros([samples interleaves pe2 averages slices 1 phases reps sets channels],'single'));
% disp(['Kspace dims: ' num2str(size(kspace))])

for i = 1:length(TSE_acq_index)
    ii = TSE_acq_index(i);
    
    d1 = raw_data.data{ii};
    d2 = complex(d1(1:2:end), d1(2:2:end));
    d3 = reshape(d2, samples, channels); %  RE & IM (2)
    
%     if nargin > 1
        d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
%     end
    
    kspace(:,...
        raw_data.head.idx.kspace_encode_step_1(ii)+1, ...
        raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
        raw_data.head.idx.average(ii)+1, ....
        raw_data.head.idx.slice(ii)+1 , ...
        1, ...
        raw_data.head.idx.phase(ii)+1, ...
        raw_data.head.idx.repetition(ii)+1, ...
        raw_data.head.idx.set(ii)+1, ...
        :) = d3;
    
end

kspace = mean(kspace,4);

    fov = FOV(1)*10; % mm
    mask = true(matrix_size(1),matrix_size(2));
    sizeMask = size(mask);
    nufft_args = {sizeMask, [6 6], 2*sizeMask, sizeMask/2, 'table', 2^12, 'minmax:kb'};

    omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
    
    omega = reshape(omega,interleaves*samples,2);
    csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);

    kx = (omega(:,1)./pi).*matrix_size(1)/2;
    ky = (omega(:,2)./pi).*matrix_size(2)/2;
    ksp = [kx(:) ky(:)]./fov;
    G = Gmri(ksp, mask, 'fov', fov, 'basis', {'dirac'}, 'nufft', nufft_args);
    csm_weights = abs(mri_density_comp(ksp, 'pipe','G',G.arg.Gnufft)); % figure, plot(csm_weights(1:samples))

img = zeros([matrix_size slices]);
csm = zeros([matrix_size channels slices]);

for iSlice = 1:slices
    data_temp = squeeze(kspace(:,:,1,1,iSlice,1,1,1,1,:));
    data_temp = reshape(data_temp, [samples*interleaves channels]);
    
    x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st)/numel(repmat(csm_weights,[1,channels]));
        csm(:,:,:,iSlice) = ismrm_estimate_csm_walsh( x );
        ccm_roemer_optimal = ismrm_compute_ccm(csm(:,:,:,iSlice), eye(channels)); % with pre-whitened
        
        %     img_tse(:,:,iTSE) = sqrt(sum(x.*conj(x),3));
    img(:,:,iSlice) = abs( sum( squeeze( x ) .* ccm_roemer_optimal, 3) );
    
end
montage_RR(img);
mag_RT.img = img;

%% TSE- Separate echo RECON
if (user_opts.sep_tse)
kspace = complex(zeros([samples interleaves pe2 1 slices contrasts phases reps sets channels],'single'));
% disp(['Kspace dims: ' num2str(size(kspace))])

for i = 1:length(TSE_acq_index)
    ii = TSE_acq_index(i);
    
    d1 = raw_data.data{ii};
    d2 = complex(d1(1:2:end), d1(2:2:end));
    d3 = reshape(d2, samples, channels); %  RE & IM (2)
    
%     if nargin > 1
        d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
%     end
    
    kspace(:,...
        raw_data.head.idx.kspace_encode_step_1(ii)+1, ...
        raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
        1, ....
        raw_data.head.idx.slice(ii)+1 , ...
        raw_data.head.idx.contrast(ii)+1, ...
        raw_data.head.idx.phase(ii)+1, ...
        raw_data.head.idx.repetition(ii)+1, ...
        raw_data.head.idx.set(ii)+1, ...
        :) = d3;
    
end

kspace = mean(kspace,4);

    fov = FOV(1)*10; % mm
    mask = true(matrix_size(1),matrix_size(2));
    sizeMask = size(mask);
    nufft_args = {sizeMask, [6 6], 2*sizeMask, sizeMask/2, 'table', 2^12, 'minmax:kb'};

    omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
    
    omega = reshape(omega,interleaves*samples,2);
    csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);

    kx = (omega(:,1)./pi).*matrix_size(1)/2;
    ky = (omega(:,2)./pi).*matrix_size(2)/2;
    ksp = [kx(:) ky(:)]./fov;
    G = Gmri(ksp, mask, 'fov', fov, 'basis', {'dirac'}, 'nufft', nufft_args);
    csm_weights = abs(mri_density_comp(ksp, 'pipe','G',G.arg.Gnufft)); % figure, plot(csm_weights(1:samples))

    csm_tse = zeros([matrix_size channels slices contrasts]);
    for iSlice = 1:slices
        img_tse = zeros([matrix_size user_opts.TSEf]);
        img_tse_phase = img_tse;
        
        for iCon = 1:contrasts
            data_temp = squeeze(kspace(:,:,1,1,iSlice,iCon,1,1,1,:));
            data_temp = reshape(data_temp, [samples*interleaves channels]);
            
            x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st)/numel(repmat(csm_weights,[1,channels]));
            csm = ismrm_estimate_csm_walsh( x );
            ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
            
            %     img_tse(:,:,iTSE) = sqrt(sum(x.*conj(x),3));
            csm_tse(:,:,:,iSlice,iCon)  = csm;
            img_tse(:,:,iCon)           = abs( sum( squeeze( x ) .* ccm_roemer_optimal, 3) );
%             img_tse_phase(:,:,iCon)     = angle( sum( x , 3) );
            
        end
        
        mag_RT.img_tse(:,:,:,iSlice) = img_tse;
        mag_RT.img_tse_phase(:,:,:,iSlice) = img_tse_phase;
    
    end
    
    % reorder slices
if slices == 1
    montage_RR(img_tse);
else
    montage_RR(RR_slice_order(squeeze(mag_RT.img_tse(:,:,1,:))));
end
end

%% TSE- Separate echo RECON
% if (user_opts.sep_tse)
% % samples     = double(raw_data.head.number_of_samples(1));
% % SLICE LOOP
% 
% % being lazy...
%     fov = FOV(1)*10; % mm
%     mask = true(matrix_size(1),matrix_size(2));
%     sizeMask = size(mask);
%     nufft_args = {sizeMask, [6 6], 2*sizeMask, sizeMask/2, 'table', 2^12, 'minmax:kb'};
% 
% %     if slices > 1
% %         % just being lazy and guessing
% %         csm_est_frames = 1:interleaves;
% %         omega = trajectory_nominal.*(pi/max(max(max(trajectory_nominal))));
% %         
% %         omega = reshape(omega,length(csm_est_frames)*samples2*2,2);
% %         csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
% %         
% %         kx = (omega(:,1)./pi).*matrix_size(1)/2;
% %         ky = (omega(:,2)./pi).*matrix_size(2)/2;
% %         ksp = [kx(:) ky(:)]./fov;
% %         G = Gmri(ksp, mask, 'fov', fov, 'basis', {'dirac'}, 'nufft', nufft_args);
% %         csm_weights = abs(mri_density_comp(ksp, 'pipe','G',G.arg.Gnufft)); % figure, plot(csm_weights(1:samples))
% %         
% %     end
%     
% for iSLC = 1:slices
%     
%     slice_ind = TSE_acq_index(find((raw_data.head.idx.slice(TSE_acq_index)+1) == iSLC));
%     data_tse = zeros(samples, length(TSE_acq_index)/slices, channels);
%     % or %     data_tse = zeros(samples, length(slice_ind), channels);
%     
%     % for i = 1:length(raw_data.data)
%     %     %     solid_int = raw_data.head.idx.kspace_encode_step_1(i)+1;
%     %     d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
%     %     d = reshape(d, samples, 1, channels);
%     %     d = ismrm_apply_noise_decorrelation_mtx(d,dmtx);
%     %     data_tse(:,i,:) = squeeze(d);
%     % end
%     
%     for i = 1:length(slice_ind)
%         %     solid_int = raw_data.head.idx.kspace_encode_step_1(i)+1;
%         d = complex(raw_data.data{slice_ind(i)}(1:2:end), raw_data.data{slice_ind(i)}(2:2:end));
%         d = reshape(d, samples, 1, channels);
%         d = ismrm_apply_noise_decorrelation_mtx(d,dmtx);
%         data_tse(:,i,:) = squeeze(d);
%     end
%     
%     data_tse = reshape(data_tse, [samples, user_opts.TSEf, (length(TSE_acq_index)/slices)/user_opts.TSEf, channels]);
%     
%     % handle ramp samples 
%     if user_opts.ramp_samples == 0
%         data_tse = data_tse(0 + [1:(samples2*2)], :,:,:);
%     else
%         % ramp samples (discard)
%         data_tse = data_tse((samples/2 + [-(samples2-1):(samples2)]), :,:,:);
%     end
%     
% %     samples = samples2*2;
%     
%     % reconstruct each echo
%     img_tse = zeros([matrix_size user_opts.TSEf]);
%     img_tse_phase = img_tse;
%     
%     for iTSE = 1:user_opts.TSEf
%         
%         % chop off ADC ringing
%         spiral_start = floor(user_opts.delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
%         %         data_temp = zeros(samples2, (length(raw_data.data)/slices)/user_opts.TSEf, channels);
%         %         data_temp(1:(1+samples2-spiral_start),:,:,:,:) = squeeze(data_tse(spiral_start:samples2,iTSE,:,:));
%         
%         data_temp = reshape(data_tse(:,iTSE,:,:), [samples2*2*(length(TSE_acq_index)/slices)/user_opts.TSEf, channels]);
%         %     data_temp = reshape(data_tse(:,iTSE,1:20,:), [samples*20, channels]);
%         
% %         if slices == 1 % this is here in case the order gets jiggled.. 
%             csm_est_frames = slice_ind(iTSE:user_opts.TSEf:(length(TSE_acq_index)/slices));
%             %     csm_est_frames = iTSE:user_opts.TSEf:(length(raw_data.data)/slices)/2;
%             
%             omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1((csm_est_frames)))+1,:)*(pi/max(max(max(trajectory_nominal))));
%             %     csm_weights=  DCF_voronoi_RR(double(omega),1,0)';
%             
%             omega = reshape(omega,length(csm_est_frames)*samples2*2,2);
%             csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
%             
%             kx = (omega(:,1)./pi).*matrix_size(1)/2;
%             ky = (omega(:,2)./pi).*matrix_size(2)/2;
%             ksp = [kx(:) ky(:)]./fov;
%             G = Gmri(ksp, mask, 'fov', fov, 'basis', {'dirac'}, 'nufft', nufft_args);
%             csm_weights = abs(mri_density_comp(ksp, 'pipe','G',G.arg.Gnufft)); % figure, plot(csm_weights(1:samples))
% %         end
%     
%         x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st)/numel(repmat(csm_weights,[1,channels]));
%         csm = ismrm_estimate_csm_walsh( x );
%         ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
%         
%         %     img_tse(:,:,iTSE) = sqrt(sum(x.*conj(x),3));
%         img_tse(:,:,iTSE) = abs( sum( squeeze( x ) .* ccm_roemer_optimal, 3) );
%         img_tse_phase(:,:,iTSE) = angle( sum( x , 3) );
%         
%     end
%     
%     mag_RT.img_tse(:,:,:,iSLC) = img_tse;
%     mag_RT.img_tse_phase(:,:,:,iSLC) = img_tse_phase;
% end
% 
% % reorder slices
% if slices == 1
%     montage_RR(img_tse);
% else
%     montage_RR(RR_slice_order(squeeze(mag_RT.img_tse(:,:,1,:))));
% end
% 
% end
% %         montage_RR(img_tse_phase,[-pi pi]);

end