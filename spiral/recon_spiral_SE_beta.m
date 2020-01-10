function [mag_RT, seheader, SNR_s] = recon_spiral_SE_beta(dfile,  nfile, SpiDes)
% [mag_RT, seheader, complex_out] = recon_spiral_SE_beta(dfile,  nfile, SpiDes)
% Spin-echo sequence reconstruction (uses user_float field to determine TSE
% factor, FOV, matrix etc).
% R Ramasawmy NHLBI Aug 2018

%% Set up
raw_data = h5read(dfile,'/dataset/data');
iRD_s = read_h5_header(dfile);

disp(['Reconstructing: ' iRD_s.measurementInformation.protocolName]);

% figure,
% subplot(3,2,1); plot(1+double(raw_data.head.idx.kspace_encode_step_1)); title('kspace step 1')
% subplot(3,2,2); plot(1+double(raw_data.head.idx.average)); title('average')
% subplot(3,2,3); plot(1+double(raw_data.head.idx.set)); title('set')
% subplot(3,2,4); plot(1+double(raw_data.head.idx.slice)); title('slice')
% subplot(3,2,5); plot(1+double(raw_data.head.idx.repetition)); title('repetition')
% subplot(3,2,6); plot(1+double(raw_data.head.idx.phase)); title('phase')

% gangles = double(max(raw_data.head.idx.kspace_encode_step_1)+1);
samples = double(raw_data.head.number_of_samples(1));
channels = double(raw_data.head.active_channels(1));
interleaves = double(max(raw_data.head.idx.kspace_encode_step_1))+1;
sets = (1 + double(max(raw_data.head.idx.set)));
reps = (1 + double(max(raw_data.head.idx.repetition)));
averages = (1 + double(max(raw_data.head.idx.average)));
slices = (1 + double(max(raw_data.head.idx.slice)));

if nargin < 3
    VDSf = 100;
    delayFactor = 0;
    gmax = 2.4;
else
    VDSf = SpiDes(1);
    delayFactor = 0;
    gmax=SpiDes(2);
end

temp = raw_data.head.user_float(:,1);
% VDSf = 100;%temp(6);

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

%% Noise
if isempty(nfile)
    dmtx = diag(ones(1,channels));
    
else
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
    
    %%
    dmtx = ismrm_dmtx_RR(nfile, dt);
    
end

%% Build Nominal Fully Sampled traj and gradients

% vds factor % temp = raw_data.head.user_int(6,1);
% FOV = 25.6;
% FOV = double(temp(6));  if FOV== 0; FOV = 25.6; end
FOV = iRD_s.encoding.reconSpace.fieldOfView_mm.x/10;

TSEf = double(temp(8)); if TSEf==0; TSEf = 1; end; clear temp;
TSEf = 4;

disp(['FOV: ' num2str(FOV)])
disp(['TSE factor: ' num2str(TSEf)])
FOV = [FOV -1*FOV*(1 - VDSf/100)];

smax = 14414.4; %
% smax = 3669.72;
krmax = 1/(2*(FOV(1)/matrix_size(1)));
[k,g] = vds(smax,gmax, dt, interleaves, FOV, krmax); close;

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
gradients_nominal = gradients_io;

%%
samples = double(raw_data.head.number_of_samples(1));

data_tse = zeros(samples, length(raw_data.data), channels);
for i = 1:length(raw_data.data)
    %     solid_int = raw_data.head.idx.kspace_encode_step_1(i)+1;
    d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
    d = reshape(d, samples, 1, channels);
    d = ismrm_apply_noise_decorrelation_mtx(d,dmtx);
    data_tse(:,i,:) = squeeze(d);
end


data_tse = reshape(data_tse, [samples, TSEf, length(raw_data.data)/TSEf, channels]);

% ramp samples (discard)
data_tse = data_tse((samples/2 + [-samples2:(samples2-1)]), :,:,:);

% data_tse = data_tse(0 + [1:(samples2*2)], :,:,:);

samples = samples2*2;

for iTSE = 1:TSEf
    
    % chop off ADC ringing
    spiral_start = floor(delayFactor*(1e-5)/dt); if spiral_start == 0; spiral_start = 1; end;
    %         data_temp = zeros(samples2, length(raw_data.data)/TSEf, channels);
    %         data_temp(1:(1+samples2-spiral_start),:,:,:,:) = squeeze(data_tse(spiral_start:samples2,iTSE,:,:));
    
    data_temp = reshape(data_tse(:,iTSE,:,:), [samples*length(raw_data.data)/TSEf, channels]);
    
    csm_est_frames = iTSE:TSEf:length(raw_data.data);
    omega = trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)*(pi/max(max(max(trajectory_nominal))));
    gradients_nominal2 = gradients_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:);
    
    omega = reshape(omega,length(csm_est_frames)*samples2*2,2);
    csm_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);
    gradients_nominal2 = reshape(gradients_nominal2, length(csm_est_frames)*size(trajectory_nominal,1),2);
    grad = complex(gradients_nominal2(:,1),gradients_nominal2(:,2));
    kk = complex(omega(:,1),omega(:,2));
    csm_weights = abs(grad(:)) .* abs(sin(angle(grad(:))-angle(kk(:))))*10; %Estimating weights from Meyer et al. Magn Reson Med. 1992 Dec;28(2):202-13.
    %         csm_weights = DCF_voronoi_RR(double(trajectory_nominal(:,double(raw_data.head.idx.kspace_encode_step_1(csm_est_frames))+1,:)),0,0);
    
    x = nufft_adj(data_temp.*repmat(csm_weights,[1,channels]), csm_st)/numel(repmat(csm_weights,[1,channels]));
    csm = ismrm_estimate_csm_walsh( x );
    ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
    
    %     img_tse(:,:,iTSE) = sqrt(sum(x.*conj(x),3));
    img_tse(:,:,iTSE) = abs( sum( squeeze( x ) .* ccm_roemer_optimal, 3) );
    img_tse_phase(:,:,iTSE) = angle( sum( x , 3) );
    
end

mag_RT.img_tse = img_tse;
mag_RT.img_tse_phase = img_tse_phase;


montage_RR(img_tse);

%         montage_RR(img_tse_phase,[-pi pi]);

end