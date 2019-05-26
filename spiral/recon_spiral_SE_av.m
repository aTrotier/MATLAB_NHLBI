function [mag_RT, seheader, SNR_s] = recon_spiral_SE_av(dfile,  nfile, SpiDes)
% [mag_RT, seheader, complex_out] = recon_spiral_SE_beta(dfile,  nfile, SpiDes)
% Spin-echo sequence reconstruction (uses user_float field to determine TSE
% factor, FOV, matrix etc).
% R Ramasawmy NHLBI Aug 2018

%% Set up
raw_data = h5read(dfile,'/dataset/data');
ismrmrd_s = read_h5_header(dfile);

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
gangles = double(max(raw_data.head.idx.kspace_encode_step_1))+1;
sets = (1 + double(max(raw_data.head.idx.set)));
reps = (1 + double(max(raw_data.head.idx.repetition)));
averages = (1 + double(max(raw_data.head.idx.average)));
slices = (1 + double(max(raw_data.head.idx.slice)));
interleaves = gangles/averages;

if nargin < 3
    %     VDSf = 100;
    %     matrix = 320;
else
    %     matrix = SpiDes(1);
    delayFactor = SpiDes;
    %     VDSf = SpiDes(2);
end

VDSf = 100;
dt = raw_data.head.sample_time_us(1)*1e-6;
matrix = ismrmrd_s.encoding.reconSpace.matrixSize.x;
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
temp = raw_data.head.user_float(:,1);
% vds factor % temp = raw_data.head.user_int(6,1);
% FOV = 25.6;
FOV = double(temp(6));  if FOV== 0; FOV = 25.6; end
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

% if slices > 1
if slices > 0
    %     slice_vec = 1:slices;
    %     if nargout ==3
    %         slice_vec = 3;
    %     end
    mag_RT = zeros([matrix_size slices]);
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
        
        data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
        data_store = data;
        
        % chop off ADC ringing
% %         spiral_start = floor(delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
% %         data = zeros(samples2, gangles, channels);
% %         data(1:(1+samples2-spiral_start),:,:) = data_store(spiral_start:samples2,:,:);
        
        omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(trajectory_nominal))));
        gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);
        
        omega = reshape(omega,length(csm_est_frames)*samples2,2);
        csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
        gradients_nominal2 = reshape(gradients_nominal2, length(csm_est_frames)*size(trajectory_nominal,1),2);
        grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
        kk = complex(omega(:,1),omega(:,2));
        csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
        
% % %         csm_weights2= DCF_voronoi_RR(double(trajectory_nominal),0,0); 
        
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
%         csm_weights2 = csm_weights2/sum(csm_weights2(:));
        
        % figure,plot([csm_weights(1:samples) csm_weights2(1:samples)])
        
        data = reshape(data,length(csm_est_frames)*samples2,channels);
        
        x = nufft_adj(data.*repmat(csm_weights,[1,channels]), csm_st);
        img_coil = squeeze(x);
        montage_RR(sqrt(sum(img_coil.*conj(img_coil),3)));
        
% % %         x = nufft_adj(data.*repmat(csm_weights2,[length(data)/length(csm_weights2),channels]), csm_st)/numel(repmat(csm_weights2,[length(data)/length(csm_weights2),channels]));
% % %         img_coil = squeeze(x);
% % %         montage_RR(sqrt(sum(img_coil.*conj(img_coil),3)));
% %         
        csm  = ismrm_estimate_csm_walsh(squeeze(img_coil));
        
        %
        %         if 0
        %             load('OCELOT_studies\d20181120_NV.mat','b0_sag')
        %             b0_map_temp = RR_slice_order(b0_sag.B0_map,5);
        %
        %             temp = ismrm_transform_image_to_kspace(mag_RT);
        %             temp = ismrm_transform_kspace_to_image(temp(160-80+1:160+80,160-80+1:160+80));
        %             mask = abs(temp) > 0.75*mean(abs(temp(:))); figure, imshow(mask);
        %
        %             montage_RR(b0_map_temp.*mask,[-150 150]); colormap('parula');
        %
        %             trajstuff.dt        = dt;
        %             trajstuff.traj      = omega;
        %             trajstuff.weights   = repmat(csm_weights2,[length(data)/length(csm_weights2),channels]);
        %             trajstuff.matrix    = matrix;
        %
        %             spiral_deblur_RR(reshape(data,samples2,length(csm_est_frames),channels),trajstuff, b0_map_temp.*mask);
        %
        %         end
        
        % +++++MS catch for pseudo-rep, do middle slice of ten++++++
        % assumed order: [6 1 7 2 8 3 9 4 10 5]
        if iSlices == 3 && nargout == 3
            % if nargout == 3
            I = ones(matrix_size(1)); D = repmat(csm_weights2, [length(data)/length(csm_weights2) size(csm,3)]);
            
            disp('Running pseudo-rep');
            pseudo_reps = 100;
            img_pr=zeros([matrix_size pseudo_reps]);
            tic;
            for i = 1:pseudo_reps
                RR_loop_count(i,pseudo_reps);
                data_pr = data + complex(randn(size(data)),randn(size(data)));
                %                 x = nufft_adj(data_pr.*repmat(csm_weights,[1,channels]), csm_st);
                x = cg_RR(data_pr(:), csm_st, I, D, csm, repmat(csm_weights2, [length(data)/length(csm_weights2) 1]), 5);
                img_coil = squeeze(x);
                img_pr(:,:,i)=sqrt(sum(img_coil.*conj(img_coil),3));
            end
            toc;
            g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
            g(g < eps) = 1;
            snr = mean(img_pr,3)./g;
            %     figure, imshow([snr],[0 50]); colorbar; colormap('parula')
            %     SNR_s.img(:,:,iSlices)   = img_pr;
            SNR_s.SNR(:,:,iSlices)   = snr;
            SNR_s.G_map(:,:,iSlices) = g;
            
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
        
        
        % %         csm_img=sqrt(sum(img_coil.*conj(img_coil),3));
        
        % %         % montage_RR(csm_img);
        
        % %         mag_RT(:,:,iSlices) = csm_img;
        
        % CG version
        %         csm  = ismrm_estimate_csm_walsh(squeeze(img_coil)); %montage_RR(abs(csm));
        
        % %         I = ones(matrix_size(1)); D = repmat(csm_weights, [1 size(csm,3)]);
        I = ones(matrix_size(1)); D = repmat(csm_weights, [length(data)/length(csm_weights) size(csm,3)]);
        
        % %         x = cg_RR(data(:), csm_st, I, D, csm, csm_weights, 5);
        x = cg_RR(data(:), csm_st, I, D, csm, repmat(csm_weights, [length(data)/length(csm_weights) 1]), 5);
        
        % test = cg_RR2(data(:), csm_st, I, D, csm, csm_weights,15); implay_RR(test(:,:,2:end));
        mag_RT(:,:,iSlices) = sqrt(sum(x.*conj(x),3));%/max(abs(x(:)));
        
    end
    
    return;
end


%% Estimate coil sensitivity
if reps > 5
    reps2 = 5;
else
    reps2 = reps;
end
csm_est_frames = reps2*interleaves*averages*slices;
% csm_est_frames = 1:csm_est_frames;
csm_est_frames = 1:TSEf*slices:csm_est_frames;
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

data = reshape(data,length(csm_est_frames)*samples2,channels);
x = nufft_adj(data.*repmat(csm_weights,[1,channels]), csm_st);
img_coil = squeeze(x);
csm_img=sqrt(sum(img_coil.*conj(img_coil),3));

montage_RR(csm_img./mean(csm_img(:)),[0 4]);

%% Calculate CSM

% McKenzie
% csm  = ismrm_estimate_csm_mckenzie(squeeze(img_coil));
% Walsh
csm  = ismrm_estimate_csm_walsh(squeeze(img_coil));
% % % csm = double(ismrm_apply_noise_decorrelation_mtx(csm,dmtx));

%montage_RR(abs(csm));

%% pseudo-rep // ismrm_non_cartesian_sense
if nargout == 3
    
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

% pseudo replica on the last loop
if nargout ==3
    start_ind = reps;
else
    start_ind = 1;
end


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
        
        data = data(1:samples2,:,:);
        data = reshape(data,length(fi)*samples2,channels);
        data = ismrm_apply_noise_decorrelation_mtx(data,dmtx);
        
        omega = trajectory_nominal(:,xi,:)*(pi/(max(trajectory_nominal(:))));
        omega = reshape(omega,length(xi)*size(trajectory_nominal,1),2); % figure, plot(omega(:,1),omega(:,2))
        
        %Calculate density compensation for nominal trajectories%
        gradients_nominal2 = reshape(gradients_nominal(:,xi,:),length(xi)*size(trajectory_nominal,1),2);
        % %     gradients_nominal2 = reshape(gradients_nominal(:,xi,:),length(xi)*length(k),2);
        grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
        kk = complex(omega(:,1),omega(:,2));
        weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
        
        % GIRF handling with axes switch
        a = corr(trajectory_nominal(:,:,1)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
        if abs(a(1)) < 0.9
            a = corr(trajectory_nominal(:,:,2)./max(trajectory_nominal(:)),trajectory_nominal_u(:,:,1)./max(trajectory_nominal_u(:)));
            if abs(a(1)) < 0.9
                disp('eh!?');
            end
            kk = complex(omega(:,2),omega(:,1));
            weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
        end
        % plot(csm_weights)
        
        st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
        
        %% CGRR (Pruessman 2001 Implementation)
        I = ones(matrix_size(1)); D = repmat(weights, [1 size(csm,3)]);
        
        %         %         [x] = cg_RR(data(:), st, I, D, csm, weights, 15);
        %         %         cg_iter(:,:,i) = sqrt(sum(x.*conj(x),3));
        
        %         xt = nufft_adj(data.*repmat(weights,[1,channels]), st);
        %         n_cg_iter(:,:,i) = sqrt(sum(xt.*conj(xt),3));
        %         %     test = cg_RR2(data(:), st, I, D, csm, weights, 15);
        %         implay_RR(test(:,:,2:end))
        
        [x(:,:,iSet)] = cg_RR(data(:), st, I, D, csm, weights, 5);
        cg_iter(:,:,i, iSet) = sqrt(sum(x.*conj(x),3));%/max(abs(x(:)));
        
        
        
    end
    %     fprintf(repmat('\b', [1 length([num2str(i) ' / ' num2str(reps)])]));
end
fprintf('\n');
mag_RT = squeeze(cg_iter(:,:,:,1));
% mag_RTn = squeeze(n_cg_iter(:,:,:,1));

%% debug testing
% data = reshape(data,length(csm_est_frames)*samples2,channels);
% 
% 
% tRR_vec = [-2:0.1:2];
% R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     
% 
% for i = 1:length(tRR_vec)
%     [traj2] = apply_GIRF(gradients_nominal, dt, R, tRR_vec(i));
%     traj2= traj2(:,:,1:2);
%     
%     omega = traj2(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(traj2))));
%     omega = reshape(omega,length(csm_est_frames)*samples2,2);
%     csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
%     csm_weights2= DCF_voronoi_RR(double(traj2)); close; close;
%     
%     x = nufft_adj(data.*repmat(csm_weights2,[1,channels]), csm_st);
%     img_coil = squeeze(x);
%     img_coil_trr(:,:,i) = sqrt(sum(img_coil.*conj(img_coil),3));
%     
% end
% 
% img_coil_trr(:,:,1:5) = rot90_stack_RR(img_coil_trr(:,:,1:5),1);
% img_coil_trr(:,:,6:15) = rot90_stack_RR(img_coil_trr(:,:,6:15),3);
% img_coil_trr(:,:,16:25) = rot90_stack_RR(img_coil_trr(:,:,16:25),1);
% img_coil_trr(:,:,26:35) = rot90_stack_RR(img_coil_trr(:,:,26:35),3);
% img_coil_trr(:,:,36:41) = rot90_stack_RR(img_coil_trr(:,:,36:41),1);
% implay_RR(img_coil_trr)
% 
% montage_RR(img_coil_trr)
% 
% 
% 
% tRR_vec = [-3:0.1:-2];
% for i = 1:length(tRR_vec)
%     [traj2] = apply_GIRF(gradients_nominal, dt, R, tRR_vec(i));
%     traj2= traj2(:,:,1:2);
%     
%     omega = traj2(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(traj2))));
%     omega = reshape(omega,length(csm_est_frames)*samples2,2);
%     csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
%     csm_weights2= DCF_voronoi_RR(double(traj2)); close; close;
%     
%     x = nufft_adj(data.*repmat(csm_weights2,[1,channels]), csm_st);
%     img_coil = squeeze(x);
%     img_coil_trr2(:,:,i) = sqrt(sum(img_coil.*conj(img_coil),3));
%     
% end
% 
% montage_RR(img_coil_trr2)
% img_coil_trr2(:,:,1:5) = rot90_stack_RR(img_coil_trr2(:,:,1:5),3);
% img_coil_trr2(:,:,6:11) = rot90_stack_RR(img_coil_trr2(:,:,6:11),1);
% implsy_RR(cat(3,img_coil_trr2(:,:,1:10), img_coil_trr))

%% other display
%  figure, imshow(mag_RT(:,:,end)./mean(mag_RT(:)),[0 4])

end
