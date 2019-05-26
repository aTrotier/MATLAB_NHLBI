function [img_s] = recon_spiral_SE_av_SNR(dfile,  nfile, SpiDes)
% [mag_RT, seheader, complex_out] = recon_spiral_SE_beta(dfile,  nfile, SpiDes)
% Spin-echo sequence reconstruction (uses user_float field to determine TSE
% factor, FOV, matrix etc).
% R Ramasawmy NHLBI Aug 2018

%% Set up

dfile = RR_run_on_mac(dfile); % incorporate with NHLBI toolbox
nfile = RR_run_on_mac(nfile);

raw_data = h5read(dfile,'/dataset/data');
iRD_s = read_h5_header(dfile);

samples = double(raw_data.head.number_of_samples(1));
gangles = double(max(raw_data.head.idx.kspace_encode_step_1))+1;
pe2 = 1+double(max(raw_data.head.idx.kspace_encode_step_2));
averages = (1 + double(max(raw_data.head.idx.average)));
slices = (1 + double(max(raw_data.head.idx.slice)));
contrasts = double(max(raw_data.head.idx.contrast))+1;
phases = double(max(raw_data.head.idx.phase))+1;
reps = (1 + double(max(raw_data.head.idx.repetition)));
sets = (1 + double(max(raw_data.head.idx.set)));
channels = double(raw_data.head.active_channels(1));
interleaves = gangles/averages;

pseudoRep_slices = [];
if nargin < 3
    VDSf = 100;
    delayFactor = 0;
    pseudoRep = 0;
else
    delayFactor = SpiDes(1);
    VDSf = SpiDes(2);
    pseudoRep = SpiDes(3);
    % pseudoRep = 0, no PR
    %             n, perform PR on slice n  
    %             <0, perform PR on all slices
    
    if pseudoRep < 0 % specify slice(s);isempty(pseudoRep)
        pseudoRep = 1;
        pseudoRep_slices = 1:slices;
    
    elseif pseudoRep > 0 % specify slice(s);
        pseudoRep_slices = pseudoRep;
        pseudoRep = 1;
        if pseudoRep_slices > slices
            pseudoRep_slices = RR_slice_order(round(slices/2));
            warning(['PR slice > #slices, using slice ' num2str(pseudoRep_slices)])
        end
    end
end

% VDSf = 100;
dt = raw_data.head.sample_time_us(1)*1e-6;
matrix = iRD_s.encoding.reconSpace.matrixSize.x;
matrix_size = [matrix matrix];

disp(['Samples: ' num2str(samples)])
disp(['Interleaves: ' num2str(interleaves)])
disp(['Channels: ' num2str(channels)])
disp(['BW: ' num2str(dt)])
disp(['Matrix: ' num2str(matrix)])
disp(['VDS: ' num2str(VDSf)])
disp(['Reps: ' num2str(reps)])
disp(['Averages: ' num2str(averages)])
disp(['Slices: ' num2str(slices)])

seheader.samples = samples;
seheader.dt = dt;
seheader.number_aqs = length(raw_data.data);
seheader.averages = averages;
seheader.channels = channels;

%% Noise
if isempty(nfile)
    dmtx = diag(ones(1,channels));
else
    dmtx = ismrm_dmtx_RR(nfile, dt);
end

%% Build Nominal Fully Sampled traj and gradients
% temp = raw_data.head.user_float(:,1);
% vds factor % temp = raw_data.head.user_int(6,1);

FOV = iRD_s.encoding.reconSpace.fieldOfView_mm.x/10;
% TSEf = double(temp(8)); if TSEf==0; TSEf = 1; end; clear temp;

disp(['FOV: ' num2str(FOV)])
% disp(['TSE factor: ' num2str(TSEf)])
FOV = [FOV -1*FOV*(1 - VDSf/100)];

smax = 14414.4; %
% smax = 3669.72;
krmax = 1/(2*(FOV(1)/matrix_size(1)));
[k,g] = vds(smax, 2.4, dt, interleaves, FOV, krmax); close;

%% Rotate Fibonacci steps
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

%% GIRF corrections
trajectory_nominal_u = trajectory_nominal;
%
% [trajectory_nominal, weights_G, st_G] = spiral_set_up_GIRF_GA(raw_data.head, FOV, matrix);

% R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix
R = [raw_data.head.phase_dir(:,1), raw_data.head.read_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix

gradients_store = gradients_nominal;
if delayFactor > 0
    spiral_start = floor(delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
    gradients_nominal = cat(1,zeros([spiral_start gangles 2]), gradients_store(1:(samples-spiral_start),:,:));
end
tRR = 0;

trajectory_nominal = apply_GIRF(gradients_nominal, dt, R, tRR );

trajectory_nominal = trajectory_nominal(:,:,1:2);

%% Grab data
kspace = zeros([samples interleaves 1 1 slices contrasts phases reps sets channels]);
% disp(['Kspace dims: ' num2str(size(kspace))])

for ii = 1:length(raw_data.data)
    
    d1 = raw_data.data{ii};
    d2 = complex(d1(1:2:end), d1(2:2:end));
    d3 = reshape(d2, samples, channels); %  RE & IM (2)
    
    if nargin > 1
        d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
    end
    
    kspace(:,...
        raw_data.head.idx.kspace_encode_step_1(ii)+1, ...
        raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
        1, .... % fix for this implementation % raw_data.head.idx.average(ii)+1, ....
        raw_data.head.idx.slice(ii)+1 , ...
        raw_data.head.idx.contrast(ii)+1, ...
        raw_data.head.idx.phase(ii)+1, ...
        raw_data.head.idx.repetition(ii)+1, ...
        raw_data.head.idx.set(ii)+1, ...
        :) = d3;
    
end

% memory management
clear raw_data;

%% recon

mag_RT = zeros([matrix_size slices]);

omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
omega = reshape(omega,size(trajectory_nominal,2)*samples2,2);
csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
%         gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);
%         gradients_nominal2 = reshape(gradients_nominal2, length(csm_est_frames)*size(trajectory_nominal,1),2);
meyer_weights = 0;
if meyer_weights
gradients_nominal2 = gradients_nominal;
gradients_nominal2 = reshape(gradients_nominal2, size(trajectory_nominal,2)*size(trajectory_nominal,1),2);
grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
kk = complex(omega(:,1),omega(:,2));
csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.

% GIRF handling with axes switch
a = corr(trajectory_nominal(:,:,1)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
if abs(a(1)) < 0.9
    a = corr(trajectory_nominal(:,:,2)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
    if abs(a(1)) < 0.9
        disp('eh!?');
    end
    kk = complex(omega(:,2),omega(:,1));
    csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
end
csm_weights = csm_weights/sum(csm_weights(:));
end
csm_weights = DCF_voronoi_RR(double(trajectory_nominal),0,0);
csm_weights = csm_weights/sum(csm_weights(:));
% % figure,plot([csm_weights(1:samples) csm_weights2(1:samples)])

for iSlices = 1:slices;%slice_vec
    disp(['Slice ' num2str(iSlices)]);
    %                    [samps ints pe2 ave sli con pha rep set cha]
    data = squeeze(kspace(  : ,  : ,  1 , : ,iSlices, 1 , 1 , 1 , 1 , : ));
    
    %         data = reshape(data,length(csm_est_frames)*samples2,channels);
    data = reshape(data,size(kspace,2)*samples2,channels);
    
    x = nufft_adj(data.*repmat(csm_weights,[1,channels]), csm_st);
    img_coil = squeeze(x);
    
    % B1-combined version
    csm = ismrm_estimate_csm_walsh( img_coil );
    ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
    mag_RT(:,:,iSlices) = abs( sum( squeeze( img_coil ) .* ccm_roemer_optimal, 3) );
    
    montage_RR(mag_RT(:,:,iSlices));
    
    % +++++MS catch for pseudo-rep, do middle slice of ten++++++
    % assumed order: [6 1 7 2 8 3 9 4 10 5]
    if ismember(iSlices, pseudoRep_slices) && pseudoRep
        
        disp('Running pseudo-rep');
        pseudo_reps = 100;
        img_pr=zeros([matrix_size pseudo_reps]);
        tic;
        for i = 1:pseudo_reps
            RR_loop_count(i,pseudo_reps);
            data_pr = data + complex(randn(size(data)),randn(size(data)));
            x = nufft_adj(data_pr.*repmat(csm_weights,[1,channels]), csm_st);
            %                 x = cg_RR(data_pr(:), csm_st, I, D, csm, repmat(csm_weights2, [length(data)/length(csm_weights2) 1]), 5);
            img_coil = squeeze(x);
            csm = ismrm_estimate_csm_walsh( img_coil );
            ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
            img_pr(:,:,i) = abs( sum( squeeze( img_coil ) .* ccm_roemer_optimal, 3) );
            
        end
        toc;
        
        g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
        g(g < eps) = 1;
        snr = mean(img_pr,3)./g;
        %     figure, subplot(1,2,1); imshow([snr],[0 50]); colorbar; colormap('parula'), subplot(1,2,2); imshow([g],[]); colorbar; colormap('parula')
        
        %     SNR_s.img(:,:,iSlices)   = img_pr;
%         SNR_s.SNR(:,:,iSlices)   = snr;
%         SNR_s.G_map(:,:,iSlices) = g;
        img_s.snr(:,:,iSlices) = snr;
    end
    % +++++++++++
    
    % CG version
    %         %         csm  = ismrm_estimate_csm_walsh(squeeze(img_coil)); %montage_RR(abs(csm));
    %         % %         I = ones(matrix_size(1)); D = repmat(csm_weights, [1 size(csm,3)]);
    %         I = ones(matrix_size(1)); D = repmat(csm_weights, [length(data)/length(csm_weights) size(csm,3)]);
    %
    %         % %         x = cg_RR(data(:), csm_st, I, D, csm, csm_weights, 5);
    %         x = cg_RR(data(:), csm_st, I, D, csm, repmat(csm_weights, [length(data)/length(csm_weights) 1]), 5);
    %
    %         % test = cg_RR2(data(:), csm_st, I, D, csm, csm_weights,15); implay_RR(test(:,:,2:end));
    %         mag_RT(:,:,iSlices) = sqrt(sum(x.*conj(x),3));%/max(abs(x(:)));
    %
end

img_s.img = mag_RT;
img_s.header = iRD_s;

end
