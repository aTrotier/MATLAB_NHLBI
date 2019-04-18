function [mag_RT, seheader, SNR_s] = recon_spiral_SE_beta(dfile,  nfile, SpiDes)
% [mag_RT, seheader, complex_out] = recon_spiral_SE_beta(dfile,  nfile, SpiDes)
% Spin-echo sequence reconstruction (uses user_float field to determine TSE
% factor, FOV, matrix etc).
% R Ramasawmy NHLBI Aug 2018

%% Set up
raw_data = h5read(dfile,'/dataset/data');
iRD_s = read_h5_header(dfile);

disp(['Reconstructing: ' iRD_s.measurementInformation.protocolName]);

figure, 
subplot(3,2,1); plot(1+double(raw_data.head.idx.kspace_encode_step_1)); title('kspace step 1')
subplot(3,2,2); plot(1+double(raw_data.head.idx.average)); title('average')
subplot(3,2,3); plot(1+double(raw_data.head.idx.set)); title('set')
subplot(3,2,4); plot(1+double(raw_data.head.idx.slice)); title('slice')
subplot(3,2,5); plot(1+double(raw_data.head.idx.repetition)); title('repetition')
subplot(3,2,6); plot(1+double(raw_data.head.idx.phase)); title('phase')

% gangles = double(max(raw_data.head.idx.kspace_encode_step_1)+1); 
samples = double(raw_data.head.number_of_samples(1));
channels = double(raw_data.head.active_channels(1));
interleaves = double(max(raw_data.head.idx.kspace_encode_step_1))+1;
sets = (1 + double(max(raw_data.head.idx.set)));
reps = (1 + double(max(raw_data.head.idx.repetition)));
averages = (1 + double(max(raw_data.head.idx.average)));
slices = (1 + double(max(raw_data.head.idx.slice)));

if nargin < 3
%     VDSf = 100;
delayFactor = 0;
else
%     VDSf = SpiDes(2);
delayFactor = SpiDes;
end

temp = raw_data.head.user_float(:,1);
VDSf = temp(6);

matrix = iRD_s.encoding.reconSpace.matrixSize.x;

dt = raw_data.head.sample_time_us(1)*1e-6;
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

%% Noise checks
disp(['Dependency ID: ' iRD_s.measurementInformation.measurementDependency.measurementID]);
asi_names = fieldnames(iRD_s.acquisitionSystemInformation);
n_iRD = read_h5_header(nfile);
nasi_names = fieldnames(n_iRD.acquisitionSystemInformation);

solid_int = 0;
coil_label = cell(channels,4);
coil_flag = zeros(channels,2);
for i = 1:length(asi_names)
    if regexp(asi_names{i}, 'coilLabel')
%         [asi_names{i} ' ' nasi_names{i}]
        solid_int = solid_int + 1;
        coil_label{solid_int,1} = iRD_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilNumber;
        coil_label{solid_int,2} = iRD_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilName;
        
        coil_label{solid_int,3} = n_iRD.acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilNumber;
        coil_label{solid_int,4} = n_iRD.acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilName;
        
%         coil_label{solid_int,5} = coil_label{solid_int,1} == coil_label{solid_int,3};
%         coil_label{solid_int,6} = regexp(coil_label{solid_int,2}, coil_label{solid_int,4});
         coil_flag(solid_int,1) = coil_label{solid_int,1} == coil_label{solid_int,3};
         coil_flag(solid_int,2) = regexp(coil_label{solid_int,2}, coil_label{solid_int,4});
    end
end

Data_CoilNum = coil_label(:,1);
Data_CoilName = coil_label(:,2);
Noise_CoilNum = coil_label(:,3);
Noise_CoilName = coil_label(:,4);

if sum(coil_flag(:,1)) ~= channels
    disp('### WARNING ### ! Coil order mismatch! (not critical?) '); disp(' ');
    disp(table(Data_CoilNum, Data_CoilName, Noise_CoilNum, Noise_CoilName)); disp(' ');
elseif sum(coil_flag(:,2)) ~= channels
    disp('### WARNING ###  ! Coil name mismatch!'); disp(' ');
    disp(table(Data_CoilNum, Data_CoilName, Noise_CoilNum, Noise_CoilName)); disp(' ');
end
clear Data_CoilName Data_CoilNum Noise_CoilName Noise_CoilNum;


%% Noise
if isempty(nfile)
    dmtx = diag(ones(1,channels));
else
    dmtx = ismrm_dmtx_RR(nfile, dt);
end

%% Build Nominal Fully Sampled traj and gradients

% vds factor % temp = raw_data.head.user_int(6,1);
% FOV = 25.6;
% FOV = double(temp(6));  if FOV== 0; FOV = 25.6; end
FOV = iRD_s.encoding.reconSpace.fieldOfView_mm.x/10;

TSEf = double(temp(8)); if TSEf==0; TSEf = 1; end; clear temp;

disp(['FOV: ' num2str(FOV)])
disp(['TSE factor: ' num2str(TSEf)])
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
trajectory_nominal = zeros(samples2,interleaves,2);
gradients_nominal =  zeros(samples2,interleaves,2);

neg = -1;
for solid_int= 1:interleaves
    rot = (solid_int-1)*(2*pi/interleaves);
    trajectory_nominal(:,solid_int,2) = neg*-( real(k(1:samples2)) *cos(rot) + imag(k(1:samples2)) *sin(rot));
    trajectory_nominal(:,solid_int,1) = neg*-(-real(k(1:samples2)) *sin(rot) + imag(k(1:samples2)) *cos(rot));
    gradients_nominal(:,solid_int,2)  = neg*-( real(g(1:samples2)) *cos(rot) + imag(g(1:samples2)) *sin(rot));
    gradients_nominal(:,solid_int,1)  = neg*-(-real(g(1:samples2)) *sin(rot) + imag(g(1:samples2)) *cos(rot));
end

%% 8-shot poet read
% % for i = 1:interleaves
% % %     temp = cumsum(xx2(i,:));
% % temp = xx2(i,:);
% % grad_poet(:,i,2) = interp1(linspace(0,1,length(temp)),temp,linspace(0,1,samples2));
% %     temp = yy2(i,:);
% % grad_poet(:,i,1) = interp1(linspace(0,1,length(temp)),temp,linspace(0,1,samples2));
% %     
% % end
% % figure, plot(grad_poet(:,:,1), grad_poet(:,:,2))
% % 
% % R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     
% % [trajectory_nominal2, gradients_nominal2] = apply_GIRF(grad_poet, dt, R); 
% % figure,plot(trajectory_nominal2(:,:,1),trajectory_nominal2(:,:,2))
% % trajectory_nominal = trajectory_nominal2(:,:,1:2);

%% GIRF corrections
trajectory_nominal_u = trajectory_nominal;
% [trajectory_nominal, weights_G, st_G] = spiral_set_up_GIRF(raw_data.head, FOV, matrix);
R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     
% tRR = -1*rem(delayFactor*1e5,dt)/dt;

% chop off ADC ringing
gradients_store = gradients_nominal;
if delayFactor > 0
    spiral_start = floor(delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
    gradients_nominal = cat(1,zeros([spiral_start interleaves 2]), gradients_store(1:(samples-spiral_start),:,:));
end
% tRR = -1*((delayFactor*(1e-5)/dt)-round(delayFactor*(1e-5)/dt));
tRR = 0;

sR.R = R;
sR.T = iRD_s.acquisitionSystemInformation.systemFieldStrength_T;
trajectory_nominal = apply_GIRF(gradients_nominal, dt, sR, tRR );
    % trajectory_nominal = apply_GIRF(gradients_nominal, dt, R, tRR );

    trajectory_nominal = trajectory_nominal(:,:,1:2);

% spiral_start = floor(delayFactor*(1e-5)/dt);

% ============
% Visualise the delay factor (scaled by factor of 100)
% 
% temp1 = dt*1e7;
%  figure, hold on,
%  plot([0:temp1:500], ones(size([0:temp1:500])), 'bo');
%  plot(spiral_start*temp1, 1, 'mo');
%  plot(repmat([0:100:500],[2 1]), [0.9*ones(size([0:100:500])); 1.1*ones(size([0:100:500]))], 'g-')
% % ============
% 
% % ===========
% % Testing adding zeros to the beginning of the traj
% 
% gradients_store = gradients_nominal;
% gradients_nominal = cat(1,zeros([spiral_start interleaves 2]), gradients_store(1:(samples-spiral_start),:,:));
% trajectory_nominal = apply_GIRF(gradients_nominal, dt, R, tRR );
% trajectory_nominal = trajectory_nominal(:,:,1:2);
% figure, plot([gradients_nominal(:,1,1) gradients_nominal(:,1,2)])
% 

% % data = zeros(samples, interleaves, channels, averages, reps);
% % for i = 1:length(raw_data.data)
% %     solid_int = raw_data.head.idx.kspace_encode_step_1(i)+1;
% %     av = raw_data.head.idx.average(i)+1;
% %     rep = raw_data.head.idx.repetition(i)+1;
% %     d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
% %     d = reshape(d, samples, 1, channels);
% %     d = ismrm_apply_noise_decorrelation_mtx(d,dmtx);
% %     data(:,solid_int,:,av,rep) = squeeze(d);
% % end
% % 
% % data_store = data;
% % data = mean(data, 4);
% % data = squeeze(data);
% % 
% % omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
% % omega = reshape(omega,interleaves*samples2,2);
% % csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
% % 
% % csm_weights = DCF_voronoi_RR(double(trajectory_nominal),0,0); % figure, plot(csm_weights);
% % 
% % iRep = 1;
% %     data_temp = reshape(data(:,:,:,iRep),interleaves*samples2,channels);
% %     x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st);
% %     img_coil = squeeze(x);
% %     
% %     csm = ismrm_estimate_csm_walsh( img_coil );
% %     ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
% %     img_ccm(:,:,iRep) = abs( sum( squeeze( img_coil ) .* ccm_roemer_optimal, 3) );
% %     
% %     montage_RR(img_ccm);

% ============

% % dt
% % 2*1e-5/dt
% % 
% % rem(1e-5,dt)
% % rem(2e-5,dt)
% % 
% % tRR_vec = -2:0.1:2;
% % R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     
% % 
% % for i = 1:length(tRR_vec)
% %     trajectory_nominal = apply_GIRF(gradients_nominal, dt, R, tRR_vec(i) );
% %     trajectory_nominal = trajectory_nominal(:,:,1:2);
% %     
% %     omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
% %     omega = reshape(omega,interleaves*samples2,2);
% %     csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
% %     
% %     csm_weights = DCF_voronoi_RR(double(trajectory_nominal),0,0); % figure, plot(csm_weights);
% %     
% %     data_temp = reshape(data,interleaves*samples2,channels);
% %     x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st);
% %     img_coil(:,:,i) = sqrt(sum(x.*conj(x),3));
% %     
% % end
% % 
% % implay_RR(img_coil);
% img_coil(:,:,[5:14,25:34]) = rot90_stack_RR(img_coil(:,:,[5:14,25:34]),2);
% img_coil(:,:,[2:11,22:31]) = rot90_stack_RR(img_coil(:,:,[2:11,22:31]),2);

% %  trajectory_nominal = trajectory_nominal_u;
%   figure,plot(trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)), trajectory_nominal_u(:,:,2)./max(trajectory_nominal_u(:)), 'b-')
%    hold on, plot(trajectory_nominal(:,:,2)./max(trajectory_nominal(:)), trajectory_nominal(:,:,1)./max(trajectory_nominal(:)), 'r-'); axis image; zoom(4);
% 
% w_u = DCF_voronoi_RR(trajectory_nominal_u);
% w_c = DCF_voronoi_RR(double(trajectory_nominal));
% 
% figure, hold on, 
% plot(w_u./median(w_u(:)), 'b-');
% plot(w_c./median(w_c(:)), 'r-');

%% Plot Trajectories
% figure, plot(trajectory_nominal(:,:,1), trajectory_nominal(:,:,2), 'b-'); axis image;

%% Slice loop?

if slices > 1
%     slice_vec = 1:slices;
%     if nargout ==3
%         slice_vec = 3;
%     end
    mag_RT = zeros([matrix_size slices]);
    mag_RT = zeros([matrix_size 1]);
  
    if reps > 5
        reps2 = 5;
    else
        reps2 = reps;
    end
    
    for iSlices = 1:slices;%slice_vec
        disp(['Slice ' num2str(iSlices)]);
        csm_est_frames = reps2*interleaves*averages*slices;
        % csm_est_frames = 1:csm_est_frames;
        csm_est_frames = iSlices:TSEf*slices:csm_est_frames;
        
        clear data; solid_int = 0;
        
        for i = csm_est_frames
            solid_int = solid_int + 1;
            d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
            d = reshape(d, samples, 1, channels);
            data(:,solid_int,:) = squeeze(d);
        end
        data = data(1:samples2,:,:);
        data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
        
        omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(trajectory_nominal))));
        gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);
        
        omega = reshape(omega,length(csm_est_frames)*samples2,2);
        csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
%         gradients_nominal2 = reshape(gradients_nominal2, length(csm_est_frames)*size(trajectory_nominal,1),2);
%         grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
%         kk = complex(omega(:,1),omega(:,2));
%         csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
% % %    
% % % %         GIRF handling with axes switch
%         a = corr(trajectory_nominal(:,:,1)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
%         if abs(a(1)) < 0.9
%             a = corr(trajectory_nominal(:,:,2)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
%             if abs(a(1)) < 0.9
%                 disp('eh!?');
%             end
%             kk = complex(omega(:,2),omega(:,1));
%             csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
%         end % figure, plot(csm_weights)
% % %     
% % %         %
% % %         data = reshape(data, samples2, length(csm_est_frames), channels); %data_hold = data;
% % %         data(1000:end,:,:) = 0;
% % %         %    
        data = reshape(data,length(csm_est_frames)*samples2,channels);
      
        x = nufft_adj(data.*repmat(csm_weights,[1,channels]), csm_st)/numel(repmat(csm_weights,[1,channels]));
        img_coil = squeeze(x);
%         csm_img=sqrt(sum(img_coil.*conj(img_coil),3));
%         csm  = ismrm_estimate_csm_walsh(squeeze(img_coil));
%         I = ones(matrix_size(1)); D = repmat(csm_weights, [1 size(csm,3)]);
%         x1 = cg_RR(data(:), csm_st, I, D, csm, csm_weights, 5);

         csm_weights2 = DCF_voronoi_RR(double(trajectory_nominal)); %close all; %figure, plot(csm_weights2)
%         x = nufft_adj(data.*repmat(csm_weights2,[length(data)/length(csm_weights2),channels]), csm_st)/numel(repmat(csm_weights2,[length(data)/length(csm_weights2),channels]));
%         img_coil = squeeze(x);
%         csm_img2=sqrt(sum(img_coil.*conj(img_coil),3));
        
        csm  = ismrm_estimate_csm_walsh(squeeze(img_coil));
        I = ones(matrix_size(1)); D = repmat(csm_weights2, [length(data)/length(csm_weights2) size(csm,3)]);
        x2 = cg_RR(data(:), csm_st, I, D, csm, repmat(csm_weights2, [length(data)/length(csm_weights2) 1]), 5);
        
%         figure, imshow([csm_img./mean(csm_img(:)) csm_img2./mean(csm_img2(:))],[]),
% figure, imshow([log(abs(ismrm_transform_kspace_to_image(fliplr(csm_img)./mean(csm_img(:))))+1) log(abs(ismrm_transform_kspace_to_image(rot90(csm_img2)./mean(csm_img2(:))))+1)],[0 1] ); colormap('parula');

%         figure, imshow([abs(x1)./mean(abs(x1(:))) abs(x2)./mean(abs(x2(:)))],[]),
        
%         ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels));
%         im_roemer_optimal = abs(sum(img_coil .* ccm_roemer_optimal, 3));
%         figure, imshow([abs(x1)./mean(abs(x1(:))) abs(x2)./mean(abs(x2(:))); abs(im_roemer_optimal)./mean(abs(im_roemer_optimal(:))) csm_img./mean(csm_img(:))],[0 4]),
%         
%         figure, imshow([abs(x2)./mean(abs(x1(:))) - abs(im_roemer_optimal)./mean(abs(im_roemer_optimal(:)));
%                         abs(x2)./mean(abs(x1(:))) - csm_img./mean(csm_img(:))],[-2 2]),
%       
        % +++++MS catch for pseudo-rep, do middle slice of ten++++++
        % assumed order: [6 1 7 2 8 3 9 4 10 5]
         if iSlices == 3 && nargout == 3
            disp('Running pseudo-rep');
            
            pseudo_reps = 100;
            img_pr=zeros([matrix_size pseudo_reps]);
            tic;
            
%             I = ones(matrix_size(1)); D = repmat(csm_weights, [1 size(csm,3)]);
            
            for i = 1:pseudo_reps
                RR_loop_count(i,pseudo_reps);
                data_pr = data + complex(randn(size(data)),randn(size(data)));
%                 x = nufft_adj(data_pr.*repmat(csm_weights2,[1,channels]), csm_st)/numel(repmat(csm_weights2,[1,channels]));
                x = cg_RR(data_pr(:), csm_st, I, D, csm, repmat(csm_weights2, [length(data)/length(csm_weights2) 1]), 5);
                
                img_coil = squeeze(x);
                img_pr(:,:,i)=sqrt(sum(img_coil.*conj(img_coil),3));
            end
            toc;
            
            g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
            g(g < eps) = 1;
            snr = mean(img_pr,3)./g;
            figure, imshow([snr],[0 50]); colorbar; colormap('parula')
            SNR_s.img   = img_pr;
            SNR_s.SNR   = snr;
            SNR_s.G_map = g;
            
         end
        % +++++++++++
        
        
% % %         dirPath = '\\hl-share.nhlbi.nih.gov\DIRHome\RamasawmyR\Scan Data\2018\180911_NV_brain\h5\';
% % %         [b0_ax_ss]=cartesian_b0_map2([dirPath 'meas_MID00258_FID04273_gre_TE410_TR50_fa25_m320.h5'], [dirPath 'meas_MID00259_FID04274_gre_TE510_TR50_fa25_m320.h5'], [dirPath 'noise219.h5'], 1); close all;
% % %         trajstuff.dt        = dt;
% % %         trajstuff.traj      = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);
% % %         trajstuff.weights   = repmat(csm_weights,[1,channels]);
% % %         trajstuff.matrix    = matrix;
% % %         
% % %         b0_map_temp= fliplr(rot90(b0_ax_ss.B0_filt));
% % %         figure, imshow(b0_map_temp,[-150 150]); colormap('parula');
% % %         temp_b= spiral_deblur_RR(reshape(data,samples2,length(csm_est_frames),channels), trajstuff, b0_map_temp);

%         
csm_img=sqrt(sum(img_coil.*conj(img_coil),3));
montage_RR(csm_img);
        
% %         mag_RT(:,:,iSlices) = csm_img;
        
    % CG version
        csm  = ismrm_estimate_csm_walsh(squeeze(img_coil)); %montage_RR(abs(csm));
        
        I = ones(matrix_size(1)); D = repmat(csm_weights, [1 size(csm,3)]);
        x = cg_RR(data(:), csm_st, I, D, csm, csm_weights, 5);
            % test = cg_RR2(data(:), csm_st, I, D, csm, csm_weights,15); implay_RR(test(:,:,2:end));
        mag_RT(:,:,iSlices) = sqrt(sum(x.*conj(x),3));%/max(abs(x(:)));
        
    end
    
    return;
end

%% TSE recon
if TSEf > 0
    
% %     [g_x_M0,g_y_M0] = vds_M0(smax,2.4,dt,interleaves,FOV,krmax); 
% %     gradients_M0 =  zeros(length(g_x_M0),interleaves,2);
% % 
% %     for i = 1:interleaves
% %         rot = (i-1)*(2*pi/interleaves);
% %         gradients_M0(:,i,1)  = neg*-( g_x_M0 *cos(rot) + g_y_M0 *sin(rot));
% %         gradients_M0(:,i,2)  = neg*-(-g_x_M0 *sin(rot) + g_y_M0 *cos(rot));
% %     end
% %     %     figure, plot(gradients_M0(:,:,1),gradients_M0(:,:,2))
% %     R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ];
% %     [traj_M0_GIRF, grad_M0_GIRF] = apply_GIRF(gradients_M0, dt, R);
% %     traj_M0_GIRF = traj_M0_GIRF(:,:,1:2);
% %     grad_M0_GIRF = grad_M0_GIRF(:,:,1:2);
% %     %     figure, plot(traj_M0_GIRF(:,:,1),traj_M0_GIRF(:,:,2))
% % % %     figure, 
% % % %     for i = 1:16
% % % %         subplot(4,4,i); hold on, plot(grad_M0_GIRF(:,i,1), 'b-'); plot(grad_M0_GIRF(:,i,2), 'r-'); title(num2str(i));
% % % %     end
% % % %     
% %     k_end_M0 = squeeze(traj_M0_GIRF(end,:,:));
% %     k_end_M0_180 = -1*k_end_M0;

% ############## average it all together to match IO ###########################################
 
data = zeros(samples, interleaves, channels, averages, reps);
for i = 1:length(raw_data.data)
    solid_int = raw_data.head.idx.kspace_encode_step_1(i)+1;
    av = raw_data.head.idx.average(i)+1;
    rep = raw_data.head.idx.repetition(i)+1;
    d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
    d = reshape(d, samples, 1, channels);
    d = ismrm_apply_noise_decorrelation_mtx(d,dmtx);
    data(:,solid_int,:,av,rep) = squeeze(d);
end
% 
% figure, 
% subplot(1,2,1); colorful_plots(abs(data(1:50,1:5,1,1))); title('Abs - First 5 readouts');
% subplot(1,2,2); colorful_plots(unwrap(angle(data(1:50,1:5,1,1)))); title('Angle - First 5 readouts');


data_store = data;

% chop off ADC ringing
% spiral_start = floor(delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
% % data = zeros(samples2, interleaves, channels,averages,reps);
% % data(1:(1+samples2-spiral_start),:,:,:,:) = data_store(spiral_start:samples2,:,:,:,:);
% % data_store = data;

data = mean(data, 4);
data = squeeze(data);
% data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);

omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
omega = reshape(omega,interleaves*samples2,2);
csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);

csm_weights = DCF_voronoi_RR(double(trajectory_nominal),0,0); % figure, plot(csm_weights);

% csm_filt= hamming(samples2*2);
% csm_filt= repmat(csm_filt((1+samples2):end),[1 interleaves channels]);
% data = data.*csm_filt;

for iRep = 1:reps
    data_temp = reshape(data(:,:,:,iRep),interleaves*samples2,channels);
    x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st);
    img_coil = squeeze(x);
    
    csm = ismrm_estimate_csm_walsh( img_coil );
    ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
    img_ccm(:,:,iRep) = abs( sum( squeeze( img_coil ) .* ccm_roemer_optimal, 3) );
    
    % montage_RR(img_ccm);
    
    csm  = ismrm_estimate_csm_walsh(squeeze(img_coil));
    I = ones(matrix_size(1)); D = repmat(reshape(csm_weights, [length(csm_weights), 1]), [length(data_temp)/length(csm_weights) size(csm,3)]);
    x = cg_RR(data_temp(:), csm_st, I, D, csm, csm_weights, 5);
    img_coil = squeeze(x);
    img_cg(:,:,iRep) = sqrt(sum(img_coil.*conj(img_coil),3));
end

if  nargout == 3
    disp('Running pseudo-rep');
    
    pseudo_reps = 100;
    img_pr=zeros([matrix_size pseudo_reps]);
    
    tic;
    
    %             I = ones(matrix_size(1)); D = repmat(csm_weights, [1 size(csm,3)]);
    
    for i = 1:pseudo_reps
        RR_loop_count(i,pseudo_reps); 
        data_pr = data_store + complex(randn(size(data_store)),randn(size(data_store)));
        data_pr = mean(data_pr,4);
        data_pr = reshape(data_pr,[interleaves*samples2 channels]);
        %         x = cg_RR(data_pr(:), csm_st, I, D, csm, repmat(csm_weights2, [length(data)/length(csm_weights2) 1]), 5);
        %         img_coil = squeeze(x);
        %         img_pr(:,:,i)=sqrt(sum(img_coil.*conj(img_coil),3));
        
        x = nufft_adj(data_pr.*repmat(csm_weights,[1,channels]), csm_st);
        img_coil = squeeze(x);
        img_temp = abs( sum( squeeze( img_coil ) .* ccm_roemer_optimal, 3) );
        img_pr(:,:,i) = img_temp;
        
    end
    toc;
    
    g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
    g(g < eps) = 1;
    snr = mean(img_pr,3)./g;
%     figure, imshow([snr],[0 50]); colorbar; colormap('parula')
    SNR_s.img   = img_pr;
    SNR_s.SNR   = snr;
    SNR_s.G_map = g;

    %%
%     I = ones(matrix_size(1)); D = repmat(csm_weights, [1 size(csm,3)]);
%     tic;
%     for i = 1:pseudo_reps
%         RR_loop_count(i,pseudo_reps);
%         data_pr = data_store + complex(randn(size(data_store)),randn(size(data_store)));
%         data_pr = mean(data_pr,4);
%         data_pr = reshape(data_pr,[interleaves*samples2 channels]);
%         x = cg_RR(data_pr(:), csm_st, I, D, csm, csm_weights, 3);
%         img_coil = squeeze(x);
%         img_pr2(:,:,i)=sqrt(sum(img_coil.*conj(img_coil),3));
%         
%     end
%     toc;
%     
%     g = std(abs(img_pr2 + max(abs(img_pr2(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
%     g(g < eps) = 1;
%     snr2 = mean(img_pr2,3)./g;
%     figure, imshow([snr snr2],[0 50]); colorbar; colormap('parula')
%     SNR_s2.img   = img_pr2;
%     SNR_s2.SNR   = snr;
%     SNR_s2.G_map = g;
%     
end

% figure, imshow(img_cg,[])
montage_RR([img_ccm./mean(img_ccm(:)) img_cg./mean(img_cg(:))],[0 4]);

%% KWIK filter
% figure, plot(raw_data.head.idx.kspace_encode_step_1)
T2_map_fit = 0;
if T2_map_fit
    
    data_tse = zeros(samples, length(raw_data.data), channels);
    for i = 1:length(raw_data.data)
        %     solid_int = raw_data.head.idx.kspace_encode_step_1(i)+1;
        d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
        d = reshape(d, samples, 1, channels);
        d = ismrm_apply_noise_decorrelation_mtx(d,dmtx);
        data_tse(:,i,:) = squeeze(d);
    end
    data_tse = reshape(data_tse, [samples, TSEf, length(raw_data.data)/TSEf, channels]);
    
    for iTSE = 1:TSEf
        
        % chop off ADC ringing
        spiral_start = floor(delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
        data_temp = zeros(samples2, length(raw_data.data)/TSEf, channels);
        data_temp(1:(1+samples2-spiral_start),:,:,:,:) = squeeze(data_tse(spiral_start:samples2,iTSE,:,:));
        
        data_temp = reshape(data_temp, [samples*length(raw_data.data)/TSEf, channels]);
        
        csm_est_frames = iTSE:TSEf:length(raw_data.data);
        omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(trajectory_nominal))));
        gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);
        
        omega = reshape(omega,length(csm_est_frames)*samples2,2);
        csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
        %     gradients_nominal2 = reshape(gradients_nominal2, length(csm_est_frames)*size(trajectory_nominal,1),2);
        %     grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
        %     kk = complex(omega(:,2),omega(:,1));
        %     csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
        csm_weights = DCF_voronoi_RR(double(trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)),0,0);
        
        x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st)/numel(repmat(csm_weights,[1,channels]));
        %     img_tse(:,:,iTSE) = sqrt(sum(x.*conj(x),3));
        img_tse(:,:,iTSE) = abs( sum( squeeze( x ) .* ccm_roemer_optimal, 3) );
        
    end
    montage_RR(img_tse); implay_RR(img_tse)
    
    % =======================
    % T2 fitting
    % =======================
    lft = fittype('a*exp(-x/b)');
    fo = fitoptions(lft);
    fo.Lower = [0 0];
    fo.Upper = [10*max(img_tse(:)) 1000];
    TE_vec = iRD_s.sequenceParameters.TE:iRD_s.sequenceParameters.TE:iRD_s.sequenceParameters.TE*TSEf;
    
    temp_roi = imclose(img_ccm > mean(img_ccm(:)),strel('disk',10)) ;  figure, imshow(img_ccm, []); hold on, contour(temp_roi, 'b');
    [tempr, tempc] = find(temp_roi);
    T2_map =zeros (size(img_ccm));
    T2_map2 =zeros (size(img_ccm));
    
    % for i = 1:length(tempr)
    %     T2_map(tempr(i), tempc(i)) = img_tse(tempr(i),tempc(i),1);
    % end
    % figure, imshow(T2_map,[])
    
    for i = 1:length(tempr)
        temp = squeeze(img_tse(tempr(i),tempc(i),:));
        fo.StartPoint = [temp(1) 100];
        [PFit] = fit(TE_vec', temp, lft,fo);
        T2_map(tempr(i), tempc(i)) = PFit.b;
        
        temp = squeeze(test_img_tse(tempr(i),tempc(i),:));
        fo.StartPoint = [temp(1) 100];
        [PFit] = fit(TE_vec', temp, lft,fo);
        T2_map2(tempr(i), tempc(i)) = PFit.b;
    end
    figure, imshow([T2_map T2_map2],[]); colorbar;
    save('temp_mat.mat');
end

% ==== intensity-based-filtering ====
% 
% weighting_vec = sin(linspace(1*pi/16, 15*pi/16, TSEf-1));
% wv2 = zeros(1,24);
% wv2(24/4+1:(3*24/4)-1) = weighting_vec; figure, plot(wv2)
% 
% for i = 1:TSEf
%     weighting_vec2 = wv2((13-i) + (0:(TSEf-1)));
%     wv_store(:,i) = weighting_vec2;
%     weighting_matrix = repmat(permute(weighting_vec2,[3 1 2]),matrix_size);
%     test_img_tse(:,:,i) = mean(img_tse.*weighting_matrix,3);
% end
% figure,colorful_plots(wv_store(:,8));
% implay_RR([rot90_stack_RR(dev.nrr(img_tse),3) rot90_stack_RR(dev.nrr(test_img_tse),3); rot90_stack_RR(dev.nrr(repmat(mean(img_tse,3),[1 1 TSEf])),3) dev.nrr(fliplr(rot90_stack_RR(repmat(crt_ax.img(:,:,2),[1 1 TSEf]),2)))]);
% 
% % a01 = img_tse(:,:,1:6);
% % a02 = test_img_tse(:,:,1:6);
% % a03 = img_tse(:,:,7:12);
% % a04 = test_img_tse(:,:,7:12);
% % figure, imshow([dev.nrr(a01(:,:));dev.nrr(a02(:,:));dev.nrr(a03(:,:));dev.nrr(a04(:,:))],[]);
% 
% % weighting_vec2 = sin(linspace(0*pi/8, pi/2, 12));
% 
% weighting_matrix = repmat(permute(weighting_vec,[3 1 2]),matrix_size);
% weighting_matrix2 = repmat(permute(weighting_vec2,[3 1 2]),matrix_size);
% 
% figure, plot([weighting_vec' weighting_vec2'])
% montage_RR([dev.nrr(mean(img_tse,3)) dev.nrr(mean(img_tse.*weighting_matrix,3)) dev.nrr(mean(img_tse.*weighting_matrix2,3))])


%
% temp_weights = csm_weights(1:samples);
%
% for iwvs = 1:12
%     csm_weights_orig = reshape(csm_weights,samples,interleaves);
%     
%     data = zeros(samples, interleaves, channels, averages, reps);
%     for i = 1:length(raw_data.data)
%         solid_int = raw_data.head.idx.kspace_encode_step_1(i)+1;
%         av = raw_data.head.idx.average(i)+1;
%         rep = raw_data.head.idx.repetition(i)+1;
%         d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
%         
%         wvi = mod(i - 1,TSEf)+1;
%         temp_weights2 = repmat(csm_weights_orig(:,raw_data.head.idx.kspace_encode_step_1(i)+1),[1 1 channels])*wv_store(wvi,iwvs);
%         
%         d = reshape(d, samples, 1, channels).*temp_weights2;
%         d = ismrm_apply_noise_decorrelation_mtx(d,dmtx);
%         data(:,solid_int,:,av,rep) = squeeze(d);
%     end
%     
%     data_store = data(1:samples2,:,:,:,:);
%     
%     spiral_start = floor(delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
%     data = zeros(samples2, interleaves, channels,averages,reps);
%     data(1:(1+samples2-spiral_start),:,:,:,:) = data_store(spiral_start:samples2,:,:,:,:);
%     data_store = data;
%     
%     data = mean(data, 4);
%     data = squeeze(data);
%     
%     iRep = 1;
%     data_temp = reshape(data(:,:,:,iRep),interleaves*samples2,channels);
%     x = nufft_adj(data_temp, csm_st);
%     img_coil = squeeze(x);
%     
%     csm = ismrm_estimate_csm_walsh( img_coil );
%     ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
%     img_ccm_kf(:,:,iwvs) = abs( sum( squeeze( img_coil ) .* ccm_roemer_optimal, 3) );
% end
% montage_RR(img_ccm_kf);
% 
% implay_RR([rot90_stack_RR(dev.nrr(img_tse),3) rot90_stack_RR(dev.nrr(test_img_tse),3) rot90_stack_RR(dev.nrr(img_ccm_kf),3); rot90_stack_RR(dev.nrr(repmat(mean(img_tse,3),[1 1 TSEf])),3) zeros(512, 512, 12) dev.nrr(fliplr_stack_RR(rot90_stack_RR(repmat(crt_ax.img(:,:,2),[1 1 TSEf]),2)))]);
% 
% 
% weight_test.img_ccm = img_ccm;

%% ############### recon separate TEs ########################################
% % % %%  
% % % csm_img = zeros([matrix_size TSEf]);
% % %     img_coil_tse = zeros([matrix_size channels TSEf]);
% % %     for iTSE = 1:TSEf
% % %         csm_est_frames = reps*interleaves*averages;
% % %         csm_est_frames = iTSE:TSEf:csm_est_frames;
% % % % %         double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1
% % %         clear data; solid_int = 0;
% % %         
% % %         for i = csm_est_frames
% % %             solid_int = solid_int + 1;
% % %             d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
% % %             d = reshape(d, samples, 1, channels);
% % %             data(:,solid_int,:) = squeeze(d);
% % %         end
% % %         data = data(1:samples2,:,:);
% % %         data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
% % % 
% % %         omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(trajectory_nominal))));
% % %         gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);
% % %         
% % %         % add on k-point to omega
% % % % %         k_corr = zeros(size(omega));
% % % % %         if iTSE > 1
% % % % %             csm_est_frames_p = iTSE-1:TSEf:(reps*interleaves*averages);
% % % % %             previous_arms = double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames_p))+1;
% % % % %             k_corr = reshape(k_end_M0_180(previous_arms,:),[1 length(csm_est_frames) 2]);
% % % % %             k_corr = repmat(k_corr, [samples 1 1]);
% % % % %         end
% % % % %         omega = omega + k_corr;
% % %         
% % %         omega = reshape(omega,length(csm_est_frames)*samples2,2);
% % %         csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
% % %         gradients_nominal2 = reshape(gradients_nominal2, length(csm_est_frames)*size(trajectory_nominal,1),2);
% % %         grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
% % %         kk = complex(omega(:,1),omega(:,2));
% % %         csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
% % %         
% % %         % - - -
% % % %         GIRF handling with axes switch
% % %         a = corr(trajectory_nominal(:,:,1)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
% % %         if abs(a(1)) < 0.9
% % %             a = corr(trajectory_nominal(:,:,2)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
% % %             if abs(a(1)) < 0.9
% % %                 disp('eh!?');
% % %             end
% % %             kk = complex(omega(:,2),omega(:,1));
% % %             csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
% % %         end
% % %         
% % %         % - - -  
% % %         data = reshape(data,length(csm_est_frames)*samples2,channels);
% % % % %         x = nufft_adj(data.*repmat(csm_weights,[1,channels]), csm_st);
% % % % %         img_coil = squeeze(x);
% % % % %          csm  = ismrm_estimate_csm_walsh(squeeze(img_coil));
% % % % %         I = ones(matrix_size(1)); D = repmat(csm_weights, [length(data)/length(csm_weights) size(csm,3)]);
% % % % %         x = cg_RR(data(:), csm_st, I, D, csm, csm_weights, 5);
% % % % %          img_coil = squeeze(x);
% % % % %         csm_img(:,:,iTSE)=sqrt(sum(img_coil.*conj(img_coil),3));
% % % %         img_coil_tse(:,:,:,iTSE)= img_coil;
% % %             
% % % %         double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames(1:length(csm_est_frames)/unique_traj_order) ))+1
% % %         
% % % %         data = reshape(data,length(csm_est_frames)*samples2,channels);
% % %         % voronoi traj 
% % %         unique_traj_order = averages/TSEf;
% % %         omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames(1:length(csm_est_frames)/unique_traj_order)))+1,:)*(pi/max(max(max(trajectory_nominal))));
% % %         
% % %         csm_weights2 = DCF_voronoi_RR(double(omega)); close all; %figure, plot(csm_weights2)
% % %         x = nufft_adj(data.*repmat(csm_weights2,[length(data)/length(csm_weights2),channels]), csm_st)/numel(repmat(csm_weights2,[length(data)/length(csm_weights2),channels]));
% % %         img_coil = squeeze(x);
% % %         
% % %         csm  = ismrm_estimate_csm_walsh(squeeze(img_coil));
% % %         I = ones(matrix_size(1)); D = repmat(csm_weights2, [length(data)/length(csm_weights2) size(csm,3)]);
% % %         x2 = cg_RR(data(:), csm_st, I, D, csm, repmat(csm_weights2, [length(data)/length(csm_weights2) 1]), 5);
% % %        
% % %         csm_img(:,:,iTSE)=sqrt(sum(x2.*conj(x2),3));
% % % %       
% % %    
% % %     end

%% original with no correction:     
%     csm_img = zeros([matrix_size TSEf]);
%     img_coil_tse = zeros([matrix_size channels TSEf]);
%     for iTSE = 1:TSEf
%         csm_est_frames = reps*interleaves*averages;
%         csm_est_frames = iTSE:TSEf:csm_est_frames;
%         clear data; solid_int = 0;
%         
%         for i = csm_est_frames
%             solid_int = solid_int + 1;
%             d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
%             d = reshape(d, samples, 1, channels);
%             data(:,solid_int,:) = squeeze(d);
%         end
%         data = data(1:samples2,:,:);
%         data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
%         
%         omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(trajectory_nominal))));
%         gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);
%         
%         omega = reshape(omega,length(csm_est_frames)*samples2,2);
%         csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
%         gradients_nominal2 = reshape(gradients_nominal2, length(csm_est_frames)*size(trajectory_nominal,1),2);
%         grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
%         kk = complex(omega(:,1),omega(:,2));
%         csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
%         
%         data = reshape(data,length(csm_est_frames)*samples2,channels);
%         x = nufft_adj(data.*repmat(csm_weights,[1,channels]), csm_st);
%         img_coil = squeeze(x);
%         img_coil_tse(:,:,:,iTSE)= img_coil;
%         csm_img(:,:,iTSE)=sqrt(sum(img_coil.*conj(img_coil),3));
%    
%     end
%% GIRF sweep
% 
% tRR_vec = -2:0.1:2;
% 
% img_ccm_gs = zeros([matrix_size length(tRR_vec)]);
% for j = 1:length(tRR_vec)
%     trajectory_nominal = apply_GIRF(gradients_nominal, dt, R, tRR_vec(j) );
%     trajectory_nominal = trajectory_nominal(:,:,1:2);
%     
%     omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
%     omega = reshape(omega,interleaves*samples2,2);
%     csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
%     
%     csm_weights = DCF_voronoi_RR(double(trajectory_nominal),0,0); % figure, plot(csm_weights);
%     
%     x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st);
%     img_coil = squeeze(x);
%     
%     csm = ismrm_estimate_csm_walsh( img_coil );
%     ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
%     img_ccm_gs(:,:,j) = abs( sum( squeeze( img_coil ) .* ccm_roemer_optimal, 3) );
%     
% end



%%  view and return  
%     montage_RR(csm_img);
    % figure, imshow(csm_img(:,:),[])
    
    mag_RT.img_ccm = img_ccm;
    mag_RT.img_cg = img_cg;
%     mag_RT.img_ccm_gs = img_ccm_gs;
    return
    
  %  img_coil = mean(img_coil_tse,4);
    
end

% % csm_est_frames = [61 46 31 16];
% % % double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1
% % % ### grab data ###
% % figure, plot(real(data(:,:,1)))
% % legend({'1','2','3','4'})

%% Long term averages
% % 
% % for iAv = 1:averages
% % csm_est_frames = find(ismember(raw_data.head.idx.average+1,iAv))';
% % clear data; solid_int = 0;
% % for i = csm_est_frames
% %     solid_int = solid_int + 1;
% %     d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
% %     d = reshape(d, samples, 1, channels);
% %     data(:,solid_int,:) = squeeze(d);
% % end
% % data = data(1:samples2,:,:);
% % data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
% % 
% % data_bucket(:,:,:,iAv) = data;
% % 
% % if iAv < 2
% % omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(trajectory_nominal))));
% % gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);
% % 
% % omega = reshape(omega,length(csm_est_frames)*samples2,2);
% % csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
% % gradients_nominal2 = reshape(gradients_nominal2, length(csm_est_frames)*size(trajectory_nominal,1),2);
% % grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
% % kk = complex(omega(:,1),omega(:,2));
% % csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
% % 
% % % GIRF handling with axes switch
% % a = corr(trajectory_nominal(:,:,1)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
% % if abs(a(1)) < 0.9
% %     a = corr(trajectory_nominal(:,:,2)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
% %     if abs(a(1)) < 0.9
% %         disp('eh!?');
% %     end
% %     kk = complex(omega(:,2),omega(:,1));
% %     csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
% % end
% % % plot(csm_weights)
% % end
% % data = reshape(data,length(csm_est_frames)*samples2,channels);
% % x = nufft_adj(data.*repmat(csm_weights,[1,channels]), csm_st);
% % img_coil = squeeze(x);
% % csm_img=sqrt(sum(img_coil.*conj(img_coil),3));
% % LTA(:,:,iAv) = squeeze(csm_img);
% % 
% % end
% % 
% % size(data_bucket) 
% % data2 = mean(data_bucket,4);
% % size(data2) 
% % x = nufft_adj(reshape(data2,length(csm_est_frames)*samples2,channels).*repmat(csm_weights,[1,channels]), csm_st);
% % img_coil = squeeze(x);
% % csm_img=sqrt(sum(img_coil.*conj(img_coil),3));
% % montage_RR(csm_img);
% % 
% % implay_RR(LTA);
% % figure, 
% % subplot(1,3,1); imshow(mean(LTA,3),[0 4*mean(LTA(:))])
% % subplot(1,3,2); imshow(std(LTA,[],3),[])
% % subplot(1,3,3); imshow(mean(LTA,3)./std(LTA,[],3),[])

%% Estimate coil sensitivity
if reps > 5
    reps2 = 5;
else
    reps2 = reps;
end
% csm_est_frames = length(raw_data.data);%
csm_est_frames = reps2*interleaves*averages*slices;
% csm_est_frames = 1:csm_est_frames;e
csm_est_frames = 1:TSEf*slices:csm_est_frames;


csm_est_frames = 1:32;
clear data; solid_int = 0;
% csm_est_frames = (find(ismember(raw_data.head.idx.average, 0)))';
for i = csm_est_frames
    solid_int = solid_int + 1;
    d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
    d = reshape(d, samples, 1, channels);
    data(:,solid_int,:) = squeeze(d);
end
data = data(1:samples2,:,:);
data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);

omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(trajectory_nominal))));
gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);

omega = reshape(omega,length(csm_est_frames)*samples2,2);
csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
gradients_nominal2 = reshape(gradients_nominal2, length(csm_est_frames)*size(trajectory_nominal,1),2);
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
% plot(csm_weights)

% assuming fully-sampled, repeated set
csm_weights = DCF_voronoi_RR(double(trajectory_nominal)); close all; %figure, plot(csm_weights2)
csm_weights = repmat(csm_weights,[size(data,2)/interleaves 1]);

data = reshape(data,length(csm_est_frames)*samples2,channels);

x = nufft_adj(data.*repmat(csm_weights,[1,channels]), csm_st);
img_coil = squeeze(x);
csm_img=sqrt(sum(img_coil.*conj(img_coil),3));

montage_RR(csm_img);

%% separate average recon

data = zeros(samples, interleaves, channels, averages, reps);
for i = 1:length(raw_data.data)
    solid_int = raw_data.head.idx.kspace_encode_step_1(i)+1;
    av = raw_data.head.idx.average(i)+1;
    rep = raw_data.head.idx.repetition(i)+1;
    d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
    d = reshape(d, samples, 1, channels);
    d = ismrm_apply_noise_decorrelation_mtx(d,dmtx);
    data(:,solid_int,:,av,rep) = squeeze(d);
end

data_store = data;

spiral_start = floor(delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
data = zeros(samples2, interleaves, channels,averages,reps);
data(1:(1+samples2-spiral_start),:,:,:,:) = data_store(spiral_start:samples2,:,:,:,:);

data = squeeze(data);

omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
omega = reshape(omega,interleaves*samples2,2);
csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);

csm_weights = DCF_voronoi_RR(double(trajectory_nominal),0,0); % figure, plot(csm_weights);

for i = 1:averages

    data_temp = reshape(data(:,:,:,i),interleaves*samples2,channels);
    x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st);
    img_coil = squeeze(x);
    img(:,:,i) = sqrt(sum(img_coil.*conj(img_coil),3));
    
end

% implay_RR(img);
% figure, plot(abs(squeeze(data(:,1,:,1))))



%%

%%

spiral_out.traj_girf = trajectory_nominal;
spiral_out.traj_nom = trajectory_nominal_u;
spiral_out.csm_img = csm_img;
spiral_out.grad_nom = gradients_nominal;

% ----------------------


[y_m0, x_m0, m0_info]=vds_M0(smax,2.4,dt,interleaves,FOV,krmax,1); close all;
gradients_M0 =  zeros(length(x_m0),interleaves,2);

neg= -1;
for i = 1:interleaves
    rot = (i-1)*(2*pi/interleaves);
    gradients_M0(:,i,1)  = neg*-( x_m0 *cos(rot) + y_m0 *sin(rot));
    gradients_M0(:,i,2)  = neg*-(-x_m0 *sin(rot) + y_m0 *cos(rot));
end

gradients_M0 = gradients_M0(:,mod(1-(1:interleaves),interleaves)+1,[1 2]);

figure, 
subplot(2,1,1); plot(gradients_M0(:,1,1), 'b-'), hold on, plot(gradients_M0(:,1,2), 'r-'); legend({'G_x','G_y'})
subplot(2,1,2); plot(gradients_M0(:,2,1), 'b-'), hold on, plot(gradients_M0(:,2,2), 'r-'); legend({'G_x','G_y'})

R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     
[traj_IO_girf, grad_IO_girf] = apply_GIRF(gradients_M0, 10e-6, R); 
% figure,plot(traj_IO_girf(:,:,1),traj_IO_girf(:,:,2))
% figure,hold on,plot(grad_IO_girf(:,1,1),'b-'),plot(grad_IO_girf(:,1,2),'r-')


extract_vec = 1:m0_info.n_spiral;
figure, plot(traj_IO_girf(:,1,1), 'b-'); hold on; plot(extract_vec,traj_IO_girf(extract_vec,1,1), 'r-');
figure, plot(grad_IO_girf(:,1,1), 'b-'); hold on; plot(extract_vec,grad_IO_girf(extract_vec,1,1), 'r-');

traj_IO_girf_e = double(traj_IO_girf(extract_vec,:,1:2));
traj_IO_girf_i = zeros(samples2,interleaves,2);
% interpolate
for ax = 1:2
    for iInt = 1:interleaves
        traj_IO_girf_i(:,iInt,ax) = interp1(linspace(0,1,length(traj_IO_girf_e)), traj_IO_girf_e(:,iInt,ax),linspace(0,1,samples2) );
    end
end


traj_IO_girf_2 = zeros(samples2,interleaves,2);

for ax = 1:2
    for iInt = 1:interleaves
        temp = traj_IO_girf_e(:,iInt,ax);
        temp2= fftshift(fft(fftshift(squeeze(temp))));
        temp3= zeros(samples2,1);
        
        temp3(samples2/2 + ( -floor(length(temp)/2):floor(length(temp)/2) )) = temp2;
        
        traj_IO_girf_2(:,iInt,ax) = abs(fftshift(ifft(ifftshift(temp3))));
    end
end
traj_store = trajectory_nominal; 
trajectory_nominal = traj_store;
trajectory_nominal = traj_IO_girf_i(:,:,1:2);
trajectory_nominal = -1*traj_IO_girf_i(:,:,[2 1]);

% ------------------------------

figure, subplot(2,3,1), plot(spiral_out.traj_girf(:,1,1),spiral_out.traj_girf(:,1,2)); title('GIRF OUT'); axis image;
        subplot(2,3,2), plot(spiral_out.traj_nom(:,1,1), spiral_out.traj_nom(:,1,2)); title('NOM OUT'); axis image;
        subplot(2,3,3), plot(trajectory_nominal(:,1,1),  trajectory_nominal(:,1,2)); title('GIRF M0'); axis image;
        subplot(2,3,4), plot(spiral_out.traj_girf(:,8,1),spiral_out.traj_girf(:,8,2)); title('GIRF OUT'); axis image;
        subplot(2,3,5), plot(spiral_out.traj_nom(:,8,1), spiral_out.traj_nom(:,8,2)); title('NOM OUT'); axis image;
        subplot(2,3,6), plot(trajectory_nominal(:,8,1),  trajectory_nominal(:,8,2)); title('GIRF M0'); axis image;

figure, subplot(2,3,1), plot(spiral_out.traj_girf(:,1,1),spiral_out.traj_girf(:,1,2)); title('GIRF OUT'); axis image;
        subplot(2,3,2), plot(spiral_out.traj_nom(:,1,1), spiral_out.traj_nom(:,1,2)); title('NOM OUT'); axis image;
        subplot(2,3,3), plot(trajectory_nominal(:,1,1),  trajectory_nominal(:,1,2)); title('GIRF M0'); axis image;
        subplot(2,3,4), plot(spiral_out.traj_girf(:,8,1),spiral_out.traj_girf(:,8,2)); title('GIRF OUT'); axis image;
        subplot(2,3,5), plot(spiral_out.traj_nom(:,8,1), spiral_out.traj_nom(:,8,2)); title('NOM OUT'); axis image;
        subplot(2,3,6), plot(trajectory_nominal(:,8,1),  trajectory_nominal(:,8,2)); title('GIRF M0'); axis image;
        
        
        
figure, plot(spiral_out.traj_girf(:,1,1),spiral_out.traj_nom(:,1,2))
figure, hold on
for i = 1:32
    plot(spiral_out.traj_girf(:,i,1),trajectory_nominal(:,i,1))
end

figure, 
subplot(1,2,1), hold on
plot(spiral_out.traj_girf(:,1,1)./max(spiral_out.traj_girf(:)),'b-'),
plot(trajectory_nominal(:,1,1)./max(trajectory_nominal(:)),'r-'),
subplot(1,2,2), hold on
plot(spiral_out.traj_girf(:,1,2)./max(spiral_out.traj_girf(:)),'b-'),
plot(trajectory_nominal(:,1,2)./max(trajectory_nominal(:)),'r-'),


figure, 
subplot(1,2,1), hold on
plot(spiral_out.traj_girf(2:end,1,1)./max(spiral_out.traj_girf(:)),'b-'),
plot(trajectory_nominal(:,1,1)./max(trajectory_nominal(:)),'r-'),
subplot(1,2,2), hold on
plot(spiral_out.traj_girf(2:end,1,2)./max(spiral_out.traj_girf(:)),'b-'),
plot(trajectory_nominal(:,1,2)./max(trajectory_nominal(:)),'r-'),


%% filtering play

% temp = 0.5+ nuttallwin(size(x,1))*nuttallwin(size(x,1))';
% % temp = 0.5+ hann(size(x,1))*hann(size(x,1))';
% temp(find(temp>1)) = 1; % figure, imshow(temp,[0 1])
% for i = 1:channels
%     yyy(:,:,i) = fftshift(fft2(fftshift(x(:,:,i))));
%     yy2(:,:,i) = temp.*yyy(:,:,i);
%     xx2(:,:,i) = fftshift(ifft2(fftshift(yy2(:,:,i))));
% end
% 
% img_coil = squeeze(xx2);
% csm_img2=sqrt(sum(img_coil.*conj(img_coil),3));
% 
% figure,imshow([csm_img./mean(csm_img(:)) csm_img2./mean(csm_img2(:))],[0 4]);
% 

%% Calculate CSM

% McKenzie
% csm  = ismrm_estimate_csm_mckenzie(squeeze(img_coil)); 
% Walsh
csm  = ismrm_estimate_csm_walsh(img_coil); 
% % % csm = double(ismrm_apply_noise_decorrelation_mtx(csm,dmtx)); 

%montage_RR(abs(csm));

%% pseudo-rep // ismrm_non_cartesian_sense
if nargout == 3
    % % reg_img=   imgaussfilt(double(csm_img > (mean(csm_img(:)) )),5)  ;
    % % reg_img(find(reg_img==0))=min(reg_img(find(reg_img>0)));% figure, imshow( reg_img )
    % % [img,snr,g,noise_psf]=ismrm_non_cartesian_sense(double(data),double(omega)/(2*pi), double(csm_weights)*(pi*(0.5^2))/sum(double(csm_weights)), csm, reg_img, 0.3,10);
    % % whos img snr g noise_psf
    % % figure,
    % % subplot(2,2,1); imshow(abs(img),[]); title('img');
    % % subplot(2,2,2); imshow(abs(snr),[]); colorbar; title('SNR');
    % % subplot(2,2,3); imshow(abs(g),[]); colorbar; title('g');
    % % subplot(2,2,4); imshow(abs(noise_psf),[]); title('noise_psf');
    
    csm_est_frames = find(double(raw_data.head.idx.repetition) + 1 == reps);
    
    clear data; solid_int = 0;
    
    for i = reshape(csm_est_frames,1,length(csm_est_frames))
        solid_int = solid_int + 1;
        d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
        d = reshape(d, samples, 1, channels);
        data(:,solid_int,:) = squeeze(d);
    end
    data = data(1:samples2,:,:);
    data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
    
    omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(trajectory_nominal))));
    gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);
    
    omega = reshape(omega,length(csm_est_frames)*samples2,2);
    csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
    gradients_nominal2 = reshape(gradients_nominal2, length(csm_est_frames)*size(trajectory_nominal,1),2);
    grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
    kk = complex(omega(:,1),omega(:,2));
    csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
    
    data = reshape(data,length(csm_est_frames)*samples2,channels);
    
    disp('Running pseudo-rep');
    pseudo_reps = 100;
    img_pr=zeros([matrix_size pseudo_reps]);
    tic;
    for i = 1:pseudo_reps
        data_pr = data + complex(randn(size(data)),randn(size(data)));
        x = nufft_adj(data_pr.*repmat(csm_weights,[1,channels]), csm_st);
        img_coil = squeeze(x);
        img_pr(:,:,i)=sqrt(sum(img_coil.*conj(img_coil),3));
    end
    toc;
    g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
    g(g < eps) = 1;
    snr = mean(img_pr,3)./g;
    
%     clear data_t
%     for iInt = 1:interleaves
%         fi = find(double(raw_data.head.idx.kspace_encode_step_1) + 1==iInt);
%         clear data_int; solid_int = 0;
%         for i = fi'
%             solid_int = solid_int + 1;
%             d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
%             d = reshape(d, samples, 1, channels);
%             data_int(:,solid_int,:) = squeeze(d);
%         end
%         data_t(:,iInt,:) = mean(data_int,2);
%     end
%     data_t = ismrm_apply_noise_decorrelation_mtx(data_t,dmtx);
%     
%     
%     % data_t = reshape(data,samples,averages,interleaves,channels); size(data_t)
%     % data_t = squeeze(mean(data_t,2)); figure, plot(squeeze(real(data_t(:,:,1)))); xlim([0 15])
%     data_t = reshape(data_t,interleaves*samples2,channels);
%     
%     omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
%     gradients_nominal2 = gradients_nominal;
%     
%     omega = reshape(omega,interleaves*samples2,2);
%     csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
%     gradients_nominal2 = reshape(gradients_nominal2, interleaves*size(trajectory_nominal,1),2);
%     grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
%     kk = complex(omega(:,1),omega(:,2));
%     csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
%     x = nufft_adj(data_t.*repmat(csm_weights,[1,channels]), csm_st);
%     
%     x = nufft_adj(data_t.*weights_G, csm_st);
%     
%     
%     tic;
%     for i = 1:pseudo_reps
%         data_pr = data_t + complex(randn(size(data_t)),randn(size(data_t)));
%         x = nufft_adj(data_pr.*weights_G, csm_st);
%         img_coil = squeeze(x);
%         img_pr(:,:,i)=sqrt(sum(img_coil.*conj(img_coil),3));
%     end
%     toc;
%     
%     g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
%     g(g < eps) = 1;
%     snr = mean(img_pr,3)./g;
    
    figure, imshow([snr],[0 50]); colorbar; colormap('parula')
    SNR_s.img   = img_pr;
    SNR_s.SNR   = snr;
    SNR_s.G_map = g;
    
    mag_RT = [];
    return;
end

%%

clear cg_iter n_cg_iter x
disp('Repetition loop'); % fprintf('\n');

% % % pseudo replica on the last loop
% % if nargout ==3
% %     start_ind = reps;
% % else

start_ind = 1;
% % end


for i = start_ind:reps
    %     disp(i)
    
    % fprintf('%d / %d',i,reps);
    RR_loop_count(i,reps);
    
    temp_window = i;
    
    
    fi_temp = find(ismember(1 + raw_data.head.idx.repetition, temp_window));
    xi_temp = double(raw_data.head.idx.kspace_encode_step_1(fi_temp)) + 1;
    
    for iSet = 1:sets
        si = find(1+raw_data.head.idx.set(fi_temp) == iSet);
        fi = fi_temp(si);
        xi = xi_temp(si);
        
        clear data; solid_int = 0;
        for j = fi'
            d = complex(raw_data.data{j}(1:2:end), raw_data.data{j}(2:2:end));
            d = reshape(d, samples, 1, channels);
            
            solid_int = solid_int + 1;
            data(:,solid_int,:) = squeeze(d);
        end
        
        
        omega = trajectory_nominal(:,xi,:)*(pi/(max(trajectory_nominal(:))));
        omega = reshape(omega,length(xi)*size(trajectory_nominal,1),2); % figure, plot(omega(:,1),omega(:,2))
        
        %Calculate density compensation for nominal trajectories%
        gradients_nominal2 = reshape(gradients_nominal(:,xi,:),length(xi)*size(trajectory_nominal,1),2);
        % %     gradients_nominal2 = reshape(gradients_nominal(:,xi,:),length(xi)*length(k),2);
        grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
        kk = complex(omega(:,1),omega(:,2));
        weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
        
% %         % GIRF handling with axes switch
% %         a = corr(trajectory_nominal(:,:,1)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
% %         if abs(a(1)) < 0.9
% %             a = corr(trajectory_nominal(:,:,2)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
% %             if abs(a(1)) < 0.9
% %                 disp('eh!?');
% %             end
% %             kk = complex(omega(:,2),omega(:,1));
% %             weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
% %         end
% %         % plot(csm_weights)
        
        weights = DCF_voronoi_RR(double(trajectory_nominal(:,xi,:))); close all; %figure, plot(csm_weights2)
%         weights = repmat(weights,[interleaves/size(data,2) 1]);
 
        data = data(1:samples2,:,:);
        data = reshape(data,length(fi)*samples2,channels);
        data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
       
        st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
        
        %% CGRR (Pruessman 2001 Implementation)
        I = ones(matrix_size(1)); D = repmat(weights, [1 size(csm,3)]);
        
        [x(:,:,iSet)] = cg_RR(data(:), st, I, D, csm, weights, 5);
 cg_iter(:,:,i) = x; %         cg_iter(:,:,i, iSet) = sqrt(sum(x.*conj(x),3));%/max(abs(x(:)));

    end
    %     fprintf(repmat('\b', [1 length([num2str(i) ' / ' num2str(reps)])]));
end
fprintf('\n');
mag_RT = squeeze(cg_iter(:,:,:,1));
% mag_RTn = squeeze(n_cg_iter(:,:,:,1));

%% attempt B0 correction post 

% 
% close all
% figure, imshow(abs(  sqrt(sum(x.*conj(x),3))    ),[]);
% %  weights_G, st_G
% 
% %x = phantom(320); whos x
% scale = sqrt(prod(prod(st_G.Kd))/numel(weights_G(:,1 )));
% datab = (nufft(x,st_G)./(sqrt(prod(st_G.Kd))))*scale;size(datab)
% 
% datab = reshape(datab, samples2,interleaves,1);
% 
% trajstuff.dt        = dt;
% trajstuff.traj      = trajectory_nominal;
% trajstuff.weights   = repmat(weights_G(1:length(datab),1),[interleaves 1]);
% trajstuff.matrix    = matrix;
% 
% spiral_deblur_RR(datab,trajstuff, b0_tmp.B0_filt);
% 
% % datab = reshape(datab,length(fi)*samples2,1);
% % 
% xxx = (nufft_adj(datab.*(weights_G(1:length(datab),1)), st_G)./(sqrt(prod(st.Kd))))*scale;
% % % % 
% % % 
% % % img_coil = squeeze(xxx);
% % % csm_img=sqrt(sum(img_coil.*conj(img_coil),3));
% % % montage_RR(csm_img./mean(csm_img(:)),[0 4]);
% % % 
% % % 
% % % %++++++++++++++++++
% % % 
% % % data = reshape(data, samples2,interleaves,length(csm_est_frames));
% % % 
% % % trajstuff.dt        = dt;
% % % trajstuff.traj      = trajectory_nominal(:,(raw_data.head.idx.kspace_encode_step_1(csm_est_frames)+1),:);
% % % trajstuff.weights   = repmat(csm_weights,[1,channels]);
% % % trajstuff.matrix    = matrix;
% % % 
% % % spiral_deblur_RR(data,trajstuff, ax_TR100_FA25.b0_map.B0_filt);

%% POET READ

load('d_POET_SE_out_24s_res256_bw250.mat'); 
temp_vec = 1:interleaves*floor(length(poet_data)/interleaves);
x = poet_data(temp_vec,3); x = reshape(x, length(temp_vec)/interleaves,interleaves);
y = poet_data(temp_vec,4); y = reshape(y, length(temp_vec)/interleaves,interleaves);

adc = poet_data(temp_vec,1); adc = reshape(adc, length(temp_vec)/interleaves,interleaves);
figure, plot(x,y)
figure, hold on,
plot(x(:,1)),plot(adc(:,1)),

a = adc(:,1);
a_ind1 = find(a, 1, 'first');
a_ind2 = find(a, 1, 'last');
xnew = x(a_ind1:a_ind2);
ynew = y(a_ind1:a_ind2);

clear gradients_M0
ind = linspace( round(length(xnew)/(samples2)/2 ),length(xnew), samples2);
x_m0 = ynew( round(ind) );
y_m0 = xnew( round(ind) );

figure, subplot(3,1,1); plot(x_m0), subplot(3,1,2); plot(y_m0), subplot(3,1,3); plot(cumsum(x_m0), cumsum(y_m0)); axis image;

for i = 1:interleaves
    rot = -1*(i-1)*(2*pi/interleaves);
    gradients_M0(:,i,1)  = ( x_m0 *cos(rot) + y_m0 *sin(rot));
    gradients_M0(:,i,2)  = (-x_m0 *sin(rot) + y_m0 *cos(rot));
end

figure, plot(gradients_M0(:,:,1), gradients_M0(:,:,2))

R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     
[traj_POET_girf, grad_POET_girf] = apply_GIRF(gradients_M0(:,:,[1 2]) , dt, R); 

% for i = 1:interleaves
%     traj_POET_girf(:,i,1) = traj_POET_girf(:,i,1) - traj_POET_girf(end,i,1);
%     traj_POET_girf(:,i,2) = traj_POET_girf(:,i,2) - traj_POET_girf(end,i,2);
% end

trajectory_POET_GIRF = traj_POET_girf(:,:,1:2); 

trajectory_POET_GIRF(:,:,1) = cumsum(gradients_M0(:,:,1));
trajectory_POET_GIRF(:,:,2) = cumsum(gradients_M0(:,:,2));

figure, plot(trajectory_POET_GIRF(:,:,1),trajectory_POET_GIRF(:,:,2))

omega = trajectory_POET_GIRF*(pi/max(max(max(trajectory_POET_GIRF))));
omega = reshape(omega,interleaves*samples2,2);
csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);

csm_weights = DCF_voronoi_RR(double(trajectory_POET_GIRF),0,0);

x = nufft_adj(data.*repmat(csm_weights,[1,channels]), csm_st);
img_coil = squeeze(x);
img_nominal_GIRF = sqrt(sum(img_coil.*conj(img_coil),3));
figure('Name', 'POET+GIRF Image'), imshow(img_nominal_GIRF,[])



%% other display
%  figure, imshow(mag_RT(:,:,end)./mean(mag_RT(:)),[0 4])
%  figure, imshow(mean(mag_RT,3),[])
%   figure, imshow([mag_RT(:,:,end)./max(mag_RT(:)) mag_RTn(:,:,end)./max(mag_RTn(:))],[])

end
