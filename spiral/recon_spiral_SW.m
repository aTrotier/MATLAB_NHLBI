function [all_res] = recon_spiral_SW(dfile, nfile, SW, cgiter, SpiDes)
% function [all_res] = recon_spiral_SW(dfile, nfile, SW, cgiter, SpiDes)
% [MAG, PHA, res, csm] = RR_manual_spiral_CGsense(datapath)
%% Set up

dfile = RR_run_on_mac(dfile); % incorporate with NHLBI toolbox
nfile = RR_run_on_mac(nfile);

raw_data = h5read(dfile,'/dataset/data');
iRD_s = read_h5_header(dfile);

gangles = double(max(raw_data.head.idx.kspace_encode_step_1)+1); % figure, plot(1+double(raw_data.head.idx.kspace_encode_step_1))
samples = double(raw_data.head.number_of_samples(1));
channels = double(raw_data.head.active_channels(1));
gareps = double(length(raw_data.head.idx.kspace_encode_step_1));
sets = (1 + double(max(raw_data.head.idx.set)));


if nargin < 5
    % custom configuration; iceprogram parameters
    
    interleaves = 8;% double(max(raw_data.head.idx.kspace_encode_step_1))+1;
    matrix = iRD_s.encoding.reconSpace.matrixSize.x;
  
else
    interleaves = SpiDes(1);
end

dt = 1e-6*double(raw_data.head.sample_time_us(1));
matrix_size = [matrix matrix];

%% Noise
if isempty(nfile)
    dmtx = diag(ones(1,channels));
else
    dmtx = ismrm_dmtx_RR(nfile, dt);
end

%% Build Nominal Fully Sampled traj and gradients

VDSf = 100; % double(raw_data.head.user_int(6,1));
FOV = iRD_s.encoding.reconSpace.fieldOfView_mm.x/10; 
FOV = [FOV -1*FOV*(1 - VDSf/100)]; disp(['FOV: ' num2str(FOV)])

krmax = 1/(2*(FOV(1)/matrix_size(1)));
[k,g] = vds(14414.4, 2.4, dt, interleaves, FOV, krmax); close;


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
% R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     
R = [raw_data.head.phase_dir(:,1), raw_data.head.read_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     

gradients_store = gradients_nominal;
delayFactor = 0; % hard-fix for the moment
if delayFactor > 0
    spiral_start = floor(delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
    gradients_nominal = cat(1,zeros([spiral_start interleaves 2]), gradients_store(1:(samples2-spiral_start),:,:));
end
tRR = 0;
% tRR = 0;

sR.R = R;
sR.T = iRD_s.acquisitionSystemInformation.systemFieldStrength_T;
trajectory_nominal = apply_GIRF(gradients_nominal, dt, sR, tRR );
% trajectory_nominal = apply_GIRF(gradients_nominal, dt, R, tRR );
trajectory_nominal = trajectory_nominal(:,:,1:2);

% figure, plot(trajectory_nominal(:,:,1), trajectory_nominal(:,:,2), 'b-');

%% Estimate coil sensitivity
csm_est_frames = 2*gangles;

clear data;
for i = 1:csm_est_frames
    d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
    d = reshape(d, samples, 1, channels);
    data(:,i,:) = squeeze(d);
end

data = ismrm_apply_noise_decorrelation_mtx(data(1:samples2,:,:),dmtx);
   
omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(1:csm_est_frames))+1,:)*(pi/max(max(max(trajectory_nominal))));
gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(1:csm_est_frames))+1,:);

omega = reshape(omega,csm_est_frames*samples2,2);
csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
gradients_nominal2 = reshape(gradients_nominal2,csm_est_frames*size(trajectory_nominal,1),2);
grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
kk = complex(omega(:,1),omega(:,2));
csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
data = reshape(data,csm_est_frames*samples2,channels);
x = nufft_adj(data.*repmat(csm_weights,[1,channels]), csm_st);
img_coil = squeeze(x);
csm_img=sqrt(sum(img_coil.*conj(img_coil),3));
% montage_RR(csm_img);

%% Mckenzie CSM
csm  = ismrm_estimate_csm_mckenzie(squeeze(img_coil)); %montage_RR(abs(csm));
% csm = double(ismrm_apply_noise_decorrelation_mtx(csm,dmtx)); %montage_RR(abs(csm));

%% Sliding window recon
%  reps = 10;
% SW = 4;

reps= (gareps/2-SW+1);
clear cg_iter n_cg_iter x
cg_iter = zeros(matrix, matrix, reps,sets);
phase_diff = zeros(matrix, matrix, reps);

for i = 1:reps % 1:gareps 
     disp(i)
     temp_window = (i*sets - 1) + (0:(sets*SW-1));
     for iSet = 1:sets
        si = find(1+raw_data.head.idx.set(temp_window) == iSet);
        fi = temp_window(si);
        xi = double(raw_data.head.idx.kspace_encode_step_1(fi)) + 1;
     
          clear data; solid_int = 0;
        for j = fi
            %      i
            d = complex(raw_data.data{j}(1:2:end), raw_data.data{j}(2:2:end));
            d = reshape(d, samples, 1, channels);
            
            solid_int = solid_int + 1;
            data(:,solid_int,:) = squeeze(d);
        end
        
        data = reshape(data(1:samples2,:,:),SW*samples2,channels);
        data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
        
        omega = trajectory_nominal(:,xi,:)*(pi/(max(trajectory_nominal(:))));
        omega = reshape(omega,length(xi)*size(trajectory_nominal,1),2); % figure, plot(omega(:,1),omega(:,2))
        
        %Calculate density compensation for nominal trajectories%
        gradients_nominal2 = reshape(gradients_nominal(:,xi,:),length(xi)*size(trajectory_nominal,1),2);
        % %     gradients_nominal2 = reshape(gradients_nominal(:,xi,:),length(xi)*length(k),2);
        grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
        kk = complex(omega(:,1),omega(:,2));
        weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
        
        st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
        
        %% CGRR (Pruessman 2001 Implementation)
        I = ones(matrix_size(1)); D = repmat(weights, [1 size(csm,3)]);
        
%         [x] = cg_RR(data(:), st, I, D, csm, weights, 15);
%         cg_iter(:,:,i) = sqrt(sum(x.*conj(x),3));
        
%         x = nufft_adj(data.*repmat(weights,[1,channels]), st);
%         n_cg_iter(:,:,i) = sqrt(sum(x.*conj(x),3));
%       [mcr, mcr2]  = cg_RR2(data(:), st, I, D, csm, weights, 20); implay_RR(mcr(:,:,2:end));
      
        [x(:,:,iSet)] = cg_RR(data(:), st, I, D, csm, weights, cgiter);
        cg_iter(:,:,i, iSet) = sqrt(sum(x.*conj(x),3))/max(abs(x(:)));
        
    end
    phase_diff(:,:,i) = atan2(imag(x(:,:,2).*conj(x(:,:,1))), real(x(:,:,2).*conj(x(:,:,1))));
     
end

% implay_RR(cg_iter(:,:,:,1)); % implay_RR(phase_diff);
all_res.mag_RT = cg_iter(:,:,:,1); %clear cg_iter;
all_res.pdiff_RT = phase_diff; %clear phase_diff;



end

