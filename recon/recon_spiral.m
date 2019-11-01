function [img_s] = recon_spiral(dfile,  nfile, SpiDes)
% function [mag_RT] = recon_spiral(dfile,  nfile, SpiDes)
% R Ramasawmy NHLBI April 2019
% SpiDes = [VDSf <100> delayFactor <0> pseudoRep <0>]
% i.e. archimedean spiral, no post-adc wait-time, and no pseudo-rep 

%% Set up

dfile = RR_run_on_mac(dfile); % incorporate with NHLBI toolbox
nfile = RR_run_on_mac(nfile);

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

interleaves = double(max(raw_data.head.idx.kspace_encode_step_1))+1;
pe2 = 1+double(max(raw_data.head.idx.kspace_encode_step_2));
contrasts = double(max(raw_data.head.idx.contrast))+1;
phases = double(max(raw_data.head.idx.phase))+1;
samples = double(raw_data.head.number_of_samples(1));
channels = double(raw_data.head.active_channels(1));
sets = (1 + double(max(raw_data.head.idx.set)));
reps = (1 + double(max(raw_data.head.idx.repetition)));
averages = (1 + double(max(raw_data.head.idx.average)));
slices = (1 + double(max(raw_data.head.idx.slice)));

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
    %             [], perform PR on all slices
    if isempty(pseudoRep)
        pseudoRep = 1;
        pseudoRep_slices = 1:slices;
    end
    if pseudoRep > 0 % specify slice(s);
        pseudoRep_slices = pseudoRep;
        pseudoRep = 1;
        if pseudoRep_slices > slices
            pseudoRep_slices = RR_slice_order(round(slices/2));
            warning(['PR slice > #slices, using slice ' num2str(pseudoRep_slices)])
        end
    end
        
end

matrix = iRD_s.encoding.reconSpace.matrixSize.x;

dt = raw_data.head.sample_time_us(1)*1e-6;
matrix_size = [matrix matrix];

disp(' ');disp('### Experiment Dimensions ###');disp(' ');
Experiment_parameters = {'Samples', 'Interleaves', 'PE2', 'Averages', 'Slices', 'Contrasts', 'Phases', 'Repetitions', 'Sets', 'Channels'}';
Value = [samples interleaves pe2 averages slices contrasts phases reps sets channels]';
disp(table( Experiment_parameters,Value )); clear Experiment_parameters Value; disp(' ');

 s_header.samples = samples;
 s_header.dt = dt;
 s_header.number_aqs = length(raw_data.data);
 s_header.averages = averages;
 s_header.channels = channels;

 iRD_s.s_header = s_header;
 
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

FOV = iRD_s.encoding.reconSpace.fieldOfView_mm.x/10;

FOV = [FOV -1*FOV*(1 - VDSf/100)]; disp(['FOV: ' num2str(FOV)])

% smax = 14414.4; % smax = 3669.72;
% gmax = 2.4;

traj_setup.gMax = iRD_s.encoding.userParameterDouble.value;
traj_setup.sMax = iRD_s.encoding.userParameterDouble_1.value;

krmax = 1/(2*(FOV(1)/matrix_size(1)));

[k,g] = vds(traj_setup.sMax, traj_setup.gMax, dt, interleaves, FOV, krmax); close;
% [k,g] = vds(smax, gmax, dt, interleaves, FOV, krmax); close;

%% Rotate 
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

%% GIRF corrections
trajectory_nominal_u = trajectory_nominal;
% R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     
R = [raw_data.head.phase_dir(:,1), raw_data.head.read_dir(:,1), raw_data.head.slice_dir(:,1)  ]; %Rotation matrix     

gradients_store = gradients_nominal;
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

%% Collect all data
kspace = complex(zeros([samples interleaves pe2 averages slices contrasts phases reps sets channels],'single'));
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
        raw_data.head.idx.average(ii)+1, ....
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

% === assuming same encoding ===
omega = trajectory_nominal*(pi/max(max(max(trajectory_nominal))));
omega = reshape(omega,interleaves*samples2,2);
recon_st = nufft_init(omega, matrix_size, [6 6],matrix_size.*2, matrix_size./2);

% spiral-out weight estimation
recon_weights = DCF_voronoi_RR(double(trajectory_nominal),0,0); % figure, plot(recon_weights);

% data_store = kspace;

if pseudoRep
    for slc = pseudoRep_slices
        %                    [samps ints pe2 ave sli con pha rep set cha]
        pr_k = squeeze(kspace(  : ,  : ,  1 , : , slc , 1 , 1 , 1 , 1 , : ));
        
        if length(size(pr_k)) == 3
            pr_k = reshape(pr_k,[size(pr_k,1) size(pr_k,2) 1 size(pr_k,3)]);
        end
        av_dim = 3;
        
        pseudo_reps = 100; disp(['Running ' num2str(pseudo_reps) ' pseudo-reps']);
        
        img_pr = zeros([matrix_size pseudo_reps]);
        
        for i = 1:pseudo_reps
            RR_loop_count(i,pseudo_reps);
            data_pr = pr_k + complex(randn(size(pr_k)),randn(size(pr_k)));
            data_pr = squeeze(mean(data_pr,av_dim));
            data_pr = reshape(data_pr,interleaves*samples2,channels);
            
            x = squeeze(nufft_adj(data_pr.*repmat(recon_weights,[1,channels]), recon_st))/sqrt(prod(recon_st.Kd));
            % calculate coil combine
            csm = ismrm_estimate_csm_walsh(x);
            ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
            
            img_pr(:,:,i) = abs( sum( x .* ccm_roemer_optimal , 3) );
        end
        
        g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
        g(g < eps) = 1;
        snr(:,:,slc) = mean(img_pr,3)./g;
        
    end
    
    img_s.snr = snr;
end

kspace = mean(kspace,4); % pseudo-rep will need to preceed this step

figure,
make_dev;
cCoil_imgs = zeros([matrix_size pe2 1 slices contrasts phases reps sets ]);
% % timer_vec = zeros(1,slices*contrasts*phases*reps*sets); tcounter = 0;

if pe2 > 1 
    % stack of spirals recon
%     kspace = ismrm_transform_kspace_to_image(kspace, 3);
    kspace = fftshift( ifft( ifftshift(kspace), pe2, 3) );

end
    
for par = 1:pe2
    for slc = 1:slices
        for coc = 1:contrasts
            for phc = 1:phases
                for repc = 1:reps
                    for setc = 1:sets
                        % %                     tic; tcounter = tcounter + 1;
                        data_temp = squeeze(kspace(1:samples2,:,par,1,slc,coc,phc,repc,setc,:));
                        
                        data_temp = reshape(data_temp,interleaves*samples2,channels);
                        x = nufft_adj(data_temp.*repmat(recon_weights,[1,channels]), recon_st)/sqrt(prod(recon_st.Kd));
                        img_coil = squeeze(x);
                        
                        csm = ismrm_estimate_csm_walsh( img_coil );
                        ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
                        cCoil_imgs(:,:,par,1,slc,coc,phc,repc,setc)= abs( sum( squeeze( img_coil ) .* ccm_roemer_optimal, 3) );
                        
                        % %                     timer_vec(tcounter) = toc;
                        
                        imshow(dev.nrr(cCoil_imgs(:,:,par,1,slc,coc,phc,repc,setc)),[0 4]); 
                        drawstring = [];
                        if pe2 > 1
                            drawstring = [drawstring 'Slice ' num2str(par) ' '];
                        end
                        if slices > 1
                            drawstring = [drawstring 'Slice ' num2str(slc) ' '];
                        end
                        if contrasts > 1
                            drawstring = [drawstring 'Contrast ' num2str(coc) ' '];
                        end
                        if phases > 1
                            drawstring = [drawstring 'Phase ' num2str(phc) ' '];
                        end
                        if reps > 1
                            drawstring = [drawstring 'Repetition ' num2str(repc) ' '];
                        end
                        if sets > 1
                            drawstring = [drawstring 'Set ' num2str(setc) ' '];
                        end
                        title(drawstring); drawnow;
                    end
                end
            end
        end
    end
end

% % disp(['Average spiral recon time: ' num2str(mean(timer_vec)) ' seconds.']);
imgs = squeeze(cCoil_imgs); close; 
% montage_RR(imgs);

%% return components
img_s.imgs = imgs;
img_s.header = iRD_s;


end
