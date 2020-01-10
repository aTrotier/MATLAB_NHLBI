%%-------------------------------%%
%%----MRI reconstruction code----%%
%%-------------------------------%%
% function [img_s,kspace,header] = recon_cartesian(file, <nfile>, <SNR_flag>)
% Fully-sampled Cartesian recon (handles partial-fourier asymmetric-echo
% and interpolation).
% input:    file = ISMRMRD .h5 data file
%           nfile = ISMRMRD .h5 noise dependency file (optional)
% output:   img_s = [structure]
%           img  = reconstructed image
%           header = ISMRMRD .xml header + modified output
%           snr = pseudo-reconstructed SNR
%           kspace = raw k-space (prewhitened if noise-file is present)
%
% R Ramasawmy Dec 2018 NHLBI

function [img_s,kspace] = recon_cartesian_flow(file, nfile, SNR_flag)
%% Read data file
if nargin == 1
    nfile = [];
end

% file = RR_run_on_mac(file); % incorporate with NHLBI toolbox
% nfile = RR_run_on_mac(nfile);

make_dev;
make_nhlbi_toolbox;

file = nhlbi_toolbox.run_path_on_sys(file);
nfile = nhlbi_toolbox.run_path_on_sys(nfile);

%%
raw_data= h5read(file, '/dataset/data');
ismrmrd_s = read_h5_header(file); disp(' ');disp('### Protocol Name ###');disp(ismrmrd_s.measurementInformation.protocolName);disp(' ');
img_s.timestamp = datetime(now, 'ConvertFrom', 'datenum'); 
header = ismrmrd_s;

samples = double(raw_data.head.number_of_samples(1));
% asymmetric echo?
asym_e = 0; if (samples~=ismrmrd_s.encoding.encodedSpace.matrixSize.x); asym_e=1; end;
echo_vec = (ismrmrd_s.encoding.encodedSpace.matrixSize.x-samples+1):ismrmrd_s.encoding.encodedSpace.matrixSize.x;
channels = double(raw_data.head.active_channels(1));

pe1 = ismrmrd_s.encoding.reconSpace.matrixSize.y; % pe1 = 1+double(max(raw_data.head.idx.kspace_encode_step_1));
pe2 = ismrmrd_s.encoding.reconSpace.matrixSize.z; % pe2 = 1+double(max(raw_data.head.idx.kspace_encode_step_2));

averages = double(max(raw_data.head.idx.average))+1;
slices = double(max(raw_data.head.idx.slice))+1;
contrasts = double(max(raw_data.head.idx.contrast))+1;
phases = double(max(raw_data.head.idx.phase))+1;
reps = double(max(raw_data.head.idx.repetition)+1);
sets = double(max(raw_data.head.idx.set))+1;
dt = raw_data.head.sample_time_us(1)*1e-6;

im_header.samples = samples;
im_header.dt = dt;
im_header.number_aqs = length(raw_data.data);
im_header.averages = averages;
im_header.channels = channels;

header.im_header = im_header;

%% Print out to window

if (samples < ismrmrd_s.encoding.encodedSpace.matrixSize.x); disp('Asymmetric Echo'); end;
% disp(['BW: ' num2str(dt)])

IM_R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ];
 
disp(' ');disp('### Acquisition Dimensions ###');disp(' ');
header_info = {'Encoded_Res','Encoded_FOV','Recon_Res','Recon_FOV'}';
X_dim = [ismrmrd_s.encoding.encodedSpace.matrixSize.x ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.x ismrmrd_s.encoding.reconSpace.matrixSize.x ismrmrd_s.encoding.reconSpace.fieldOfView_mm.x]';
Y_dim = [ismrmrd_s.encoding.encodedSpace.matrixSize.y ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.y ismrmrd_s.encoding.reconSpace.matrixSize.y ismrmrd_s.encoding.reconSpace.fieldOfView_mm.y]';
Z_dim = [ismrmrd_s.encoding.encodedSpace.matrixSize.z ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.z ismrmrd_s.encoding.reconSpace.matrixSize.z ismrmrd_s.encoding.reconSpace.fieldOfView_mm.z]';
disp(table(header_info, X_dim, Y_dim, Z_dim)); clear header_info X_dim Y_dim Z_dim;

disp(' ');disp('### Experiment Dimensions ###');disp(' ');
Experiment_parameters = {'RO', 'PE1', 'PE2', 'Averages', 'Slices', 'Contrasts', 'Phases', 'Reps', 'Sets', 'Channels'}';
Value = [ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 pe2 averages slices contrasts phases reps sets channels]';
disp(table( Experiment_parameters,Value )); clear Experiment_parameters Value; disp(' ');

if sets < 2
   warning('Only a single set acquisition, reconstructing mag + phase'); 
end

%% Calculate Noise Decorrelation Matrix

if isempty(nfile)
    disp('No noise adjustment: ');
    dmtx = diag(ones(1,channels));
else
    dmtx = nhlbi_toolbox.noise_adjust(nfile, ismrmrd_s, raw_data.head.sample_time_us(1)*1e-6, nhlbi_toolbox);
end

disp(' ');

%% Recon

ismrmrd_s.userParameters.userParameterLong_1.name, ismrmrd_s.userParameters.userParameterLong_1.value
bin_phases = ismrmrd_s.userParameters.userParameterLong_1.value;

[bin_data, NominalInterval] = physio_Binning(raw_data.head.physiology_time_stamp, bin_phases); % figure, plot(raw_data.head.physiology_time_stamp(1,:), 'r-');
   
kspace = complex(zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 pe2 1 slices contrasts bin_phases reps sets channels],'single'));

for iPha = 1:bin_phases
    temp = bin_data{iPha};

    % == average each line, each set together ==
    temp_lines = single(raw_data.head.idx.kspace_encode_step_1(temp)) + 1;
    temp_sets = single(raw_data.head.idx.set(temp)) + 1;
    
    for iSet = 1:sets
        s_ind = find(temp_sets == iSet);
        temp_lines2 = temp_lines(s_ind);
        
        for iPE1 = 1:pe1
            
            % == probably a cleaner way than using for-loops and finds ==
            current_line = find(temp_lines2 == iPE1);
            
            if isempty(current_line)
                d4 = zeros(samples, channels);
            else
                d4 = zeros(samples,channels,length(current_line));
                
                for iCL = 1:length(current_line)
                    
                    ii =  temp(s_ind(current_line(iCL)));
                    d1 = raw_data.data{ii};
                    d2 = complex(d1(1:2:end), d1(2:2:end));
                    d3 = reshape(d2, samples, channels); %  RE & IM (2)
                    d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
                    d4(:,:,iCL) = d3;
                end
                d4 = squeeze(mean(d4,3));
                
            end
            
        kspace(echo_vec,...
            iPE1, ... % raw_data.head.idx.kspace_encode_step_1(ii)+1
            raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
            1, .... %raw_data.head.idx.average(ii)+1
            raw_data.head.idx.slice(ii)+1 , ...
            raw_data.head.idx.contrast(ii)+1, ...
            iPha, ...
            raw_data.head.idx.repetition(ii)+1, ...
            iSet, ... % raw_data.head.idx.set(ii)+1
            :) = d4;
        
        end
    end
    
end

% size(kspace)


%% SNR
if nargin < 3
    SNR_flag = 0;
end

if SNR_flag
    
    if SNR_flag > 0 % specify slice(s);
        pseudoRep_phases = SNR_flag;
    else
        pseudoRep_phases = 1:bin_phases;
    end
    
    %     snr = zeros(ismrmrd_s.encoding.reconSpace.matrixSize.x, ismrmrd_s.encoding.reconSpace.matrixSize.y, slices);
    snr = zeros( ismrmrd_s.encoding.reconSpace.matrixSize.x, size(kspace,2), length(pseudoRep_phases) );
    
    
    for i = 1:length(pseudoRep_phases)
        
        iPhase = pseudoRep_phases(i);
        %         snr(:,:,i) = cartesian_pseudoreps(squeeze(kspace(:,:,:,:,1,1,iPhase,1,1,:)));
        
        temp = bin_data{iPhase};
        
        temp_lines = single(raw_data.head.idx.kspace_encode_step_1(temp)) + 1;
        temp_sets = single(raw_data.head.idx.set(temp)) + 1;
        
        pseudo_reps = 100;
        disp(['Running ' num2str(pseudo_reps) ' pseudo-reps']);
        
        img_pr = zeros(samples, size(kspace,2), pseudo_reps);
        
        for iPR = 1:pseudo_reps
%             RR_loop_count(iPR,pseudo_reps);
            
            % x y z 1-average (with binning method) [assume slice,
            % contrasts and reps = 1], and choose set = 1;
            kspace_pr = complex(zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 pe2 1 slices contrasts 1 reps 1 channels],'single'));
            
            s_ind = find(temp_sets == 1);
            temp_lines2 = temp_lines(s_ind);
            
            for iPE1 = 1:pe1
                current_line = find(temp_lines2 == iPE1);
                
                if isempty(current_line)
                    d4 = complex(randn(samples, channels), randn(samples, channels));
                    
                else
                    d4 = zeros(samples,channels,length(current_line));
                    
                    for iCL = 1:length(current_line)
                        
                        ii =  temp(s_ind(current_line(iCL)));
                        d1 = raw_data.data{ii};
                        d2 = complex(d1(1:2:end), d1(2:2:end));
                        d3 = reshape(d2, samples, channels); %  RE & IM (2)
                        d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
                        d4(:,:,iCL) = d3;
                    end
                    
                    d4 = d4 + complex(randn(size(d4)),randn(size(d4)));
                    
                    d4 = squeeze(mean(d4,3));
                    
                end
                
                kspace_pr(echo_vec,...
                    iPE1, ... % raw_data.head.idx.kspace_encode_step_1(ii)+1
                    raw_data.head.idx.kspace_encode_step_2(ii)+1, ...
                    1, ....
                    raw_data.head.idx.slice(ii)+1 , ...
                    raw_data.head.idx.contrast(ii)+1, ...
                    1, ...
                    raw_data.head.idx.repetition(ii)+1, ...
                    1, ... % raw_data.head.idx.set(ii)+1
                    :) = d4;
                
            end
            
            
            kspace_pr = squeeze(kspace_pr);
            
            if ismrmrd_s.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 == 1
                
                x = ismrm_transform_kspace_to_image(kspace_pr,[1 2]);
                
            else
                % detect which lines are reference embedded:
                sampling_pattern = kspace(:,:,1,1,1,1,1,1,1,1); % use first coil
                sampling_pattern = abs(sampling_pattern) > 0;  %
                sampling_pattern_1 = sampling_pattern(1,:);
                acc = ismrmrd_s.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1;
                
                temp0 = find(sampling_pattern_1);
                temp1 = diff(temp0);
                temp2 = temp0(find(temp1==1));
                sampling_pattern = single(sampling_pattern);
                
                sampling_pattern(:,temp2) = sampling_pattern(:,temp2)*3;
                
                %                 img_pr(:,:,iPR) = abs(ismrm_cartesian_GRAPPA( kspace_pr ,sampling_pattern, acc));
                
                [~,ksp] = GRAPPA(kspace_pr,sampling_pattern, acc);
                % if we want to pad the kspace (asym echo etc)
                x = ismrm_transform_kspace_to_image(ksp,[1 2]);
            end
            
            csm = ismrm_estimate_csm_walsh(x);
            ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels));
            
            img_pr(:,:,iPR) = abs(sum(x .* ccm_roemer_optimal, 3));
                
            
            
        end
        
        img_pr = img_pr( size(img_pr,1)/4:size(img_pr,1)/4 +size(img_pr,1)/2 -1, :,:);
        
        g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
        g(g < eps) = 1;
        snr(:,:,i) = mean(img_pr,3)./g;
        
    end
%     figure, imagesc(snr(:,:),[0 50]); colorbar;
end


%% transform to image-space

kspace = mean(kspace,4);% kspace = squeeze(kspace);cCoil_imgs = ismrm_transform_kspace_to_image(kspace);

cCoil_imgs = zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 pe2 1 slices contrasts bin_phases reps sets channels],'single');

asym_e = 0; % RR override for laziness..
if asym_e == 1
    cCoil_imgs = zeros(size(kspace));
    for slc = 1:slices
        for coc = 1:contrasts
            for phc = 1:phases
                for repc = 1:reps
                    for setc = 1:sets
                        for i = 1:channels
                            
                            
                            % POCS recon (scaling?) - being lazy and doing
                            % it per coil, though the recon can handle
                            % Nc x Ny x Nx x Nz data
                            
                            temp = squeeze(kspace(:,:,:,1,slc,coc,phc,repc,setc,i));
                            if length(size(temp))==2
                                temp = permute(temp,[2 1]);
                                pocs_recon = pocs(temp,10); % pocs_recon = pocs(temp,10, true);
                                cCoil_imgs(:,:,1,1,slc,coc,phc,repc,setc,i) = permute(pocs_recon,[2 1]);
                            else
                                temp = permute(temp,[2 1 3]);
                                pocs_recon = pocs(temp,10);
                                cCoil_imgs(:,:,:,1,slc,coc,phc,repc,setc,i) = permute(pocs_recon,[2 1 3]);
                            end
                            
                            
                            
                        end
                    end
                end
            end
        end
    end
else
    
    if ismrmrd_s.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 == 1
        cCoil_imgs = ismrm_transform_kspace_to_image(kspace,[1 2]);
        
    else
        
        disp(['PAT factor: ' num2str(ismrmrd_s.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1) ' -- GRAPPA type: ' ismrmrd_s.encoding.parallelImaging.calibrationMode, ' -- num reference lines: '  num2str(ismrmrd_s.userParameters.userParameterLong_3.value) ])
        %         ismrmrd_s.userParameters.userParameterLong_3.name
        %          ismrmrd_s.userParameters.userParameterLong_3.value
        
        % detect which lines are reference embedded:
        sampling_pattern = squeeze(kspace(:,:,:,1,1,1,1,1,1,1));
        sampling_pattern = abs(sampling_pattern) > 0;  %
        %         figure, imshow(sampling_pattern)
        
        sampling_pattern_1 = sampling_pattern(1,:);
        acc = ismrmrd_s.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1;
        
        temp0 = find(sampling_pattern_1);
        temp1 = diff(temp0);
        temp2 = temp0(find(temp1==1));
        sampling_pattern = single(sampling_pattern);
        
        sampling_pattern(:,temp2) = sampling_pattern(:,temp2)*3;
                   
        % === GRAPPA recon  === >> 
        for iPha = 1:bin_phases
            for iSet = 1:sets
                
%                 cCoil_imgs(:,:,:,:,:,:,iPha,:,iSet) = ismrm_cartesian_GRAPPA(squeeze( kspace(:,:,:,:,:,:,iPha,:,iSet,:)),sampling_pattern, acc) ;

                [~,ksp] = GRAPPA(squeeze( kspace(:,:,:,:,:,:,iPha,:,iSet,:)),sampling_pattern, acc);
                % if we want to pad the kspace (asym echo etc)
                cCoil_imgs(:,:,:,:,:,:,iPha,:,iSet,:) = ismrm_transform_kspace_to_image(ksp,[1 2]);

            end
        end
        % === GRAPPA recon  === <<
        
    end
    
end



%% remove OS
OverSampling = ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.x./ismrmrd_s.encoding.reconSpace.fieldOfView_mm.x;
if OverSampling == 2
    
    xdim = ismrmrd_s.encoding.encodedSpace.matrixSize.x;
    temp = reshape(cCoil_imgs, [xdim prod(size(cCoil_imgs))/xdim]);
    temp = temp( (xdim/4): (xdim/4 + xdim/2 -1), : );
    new_dims = size(cCoil_imgs); new_dims(1) = ismrmrd_s.encoding.reconSpace.matrixSize.x;
    cCoil_imgs = reshape(temp, new_dims);
    
end


%% 2D Phase Difference
% size(cCoil_imgs)

% ccd = cCoil_imgs(:,:,:,:,:,:,:,:,1,:).*conj(cCoil_imgs(:,:,:,:,:,:,:,:,2,:)); size(ccd); CD = squeeze(sum(ccd,10)); size(CD)
% PD = angle(CD);
% implay_RR(PD)

%     % Optimised coil combination (Roemer): Using data as CSM
PD_ccm = zeros(size(cCoil_imgs,1), size(cCoil_imgs,2), bin_phases);
img_ccm = zeros(size(cCoil_imgs,1), size(cCoil_imgs,2), bin_phases);

% if ismrmrd_s.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 == 1
if sets == 2
    
    for iPha = 1:bin_phases
        csm = ismrm_estimate_csm_walsh( squeeze(mean(cat(9, cCoil_imgs(:,:,:,:,:,:,iPha,:,1,:), cCoil_imgs(:,:,:,:,:,:,iPha,:,2,:)  ),9)) );
        ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened data
        CCM1 = (sum(squeeze(cCoil_imgs(:,:,:,:,:,:,iPha,:,1,:)) .* ccm_roemer_optimal, 3));
        CCM2 = (sum(squeeze(cCoil_imgs(:,:,:,:,:,:,iPha,:,2,:)) .* ccm_roemer_optimal, 3));
        PD_ccm(:,:,iPha) = angle(CCM1.*conj(CCM2));
        
        temp = squeeze(cCoil_imgs(:,:,:,:,:,:,iPha,:,1,:));
        img_ccm(:,:,iPha) = abs( sum( squeeze( temp ) .* ccm_roemer_optimal, 3) );
    end
    
    % implay_RR([PD PD_ccm])
    
elseif sets == 1
    % cheeky attempt at phase recon.. 
    for iPha = 1:bin_phases
        temp = squeeze(cCoil_imgs(:,:,:,:,:,:,iPha,:,1,:));
        
        csm = ismrm_estimate_csm_walsh(temp);
        ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened data
        
        CCM1 = ( sum( temp .* ccm_roemer_optimal, 3 ) );
        
        PD_ccm(:,:,iPha) = angle(CCM1);
        img_ccm(:,:,iPha) = abs( CCM1 );
        
    end
    
end

% else
%     for iPha = 1:bin_phases
%         
%         CCM1 = squeeze(cCoil_imgs(:,:,:,:,:,:,iPha,:,1,:));
%         CCM2 = squeeze(cCoil_imgs(:,:,:,:,:,:,iPha,:,2,:));
%         
%         PD_ccm(:,:,iPha) = angle(CCM1.*conj(CCM2));
%         img_ccm(:,:,iPha) = abs(CCM1);
%         
%     end
% end

% dev.implay_flow(img_ccm, PD_ccm);

%% Return variables

img_s.img = img_ccm;
img_s.pd  = PD_ccm;
img_s.header = header;

img_s.header.im_header.NominalInterval = NominalInterval;

if SNR_flag
    img_s.snr = snr;
end

end

function [ccc] = binning_stats(raw_data, bin_data)

%% % Output stats: Count number of arms per frame
debug = 0;
gangles = 1 + double(max(raw_data.head.idx.kspace_encode_step_1));

if debug
    plot_bin_fig = figure;
    plot_bin_tabgp = uitabgroup(plot_bin_fig); clear a tab;
end


for i = 1:30
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
    
%     [counts] = hist(t_arms0, 0.5+(1:gangles));
    [counts] = hist(t_arms0, (1:gangles));
    ccc(:,i) = counts;
end

%  figure, imagesc(ccc, [0 max(ccc(:))])
%  temp_str = ['Min #unique arms: ' num2str(min(ccc(:))) ' Max: ' num2str(max(ccc(:))) ]; % interleaves/min(bbb(:))
%  xlabel('Frames'); ylabel('Gangles'); title(temp_str); shg;

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


function snr = cartesian_pseudoreps(crt_k)
% SNR_s = cartesian_psuedoreps(crt_k, slice)
% Expects pre-whitened 2D k-space data [RO,PE,AV,COILS]
% This can be expanded to generic dimensions, but will require header info
% on dimension - better used as a subfunction within recon_cartesian
% R Ramasawmy NHLBI 2019

if length(size(crt_k)) == 3
    crt_k = reshape(crt_k,[size(crt_k,1) size(crt_k,2) 1 size(crt_k,3)]);
end
% size(crt_k)
av_dim = 3;

pseudo_reps = 100;
disp(['Running ' num2str(pseudo_reps) ' pseudo-reps']);

ncoils = size(crt_k,ndims(crt_k));
img_pr = zeros(size(crt_k,1), size(crt_k,2), pseudo_reps);

for i = 1:pseudo_reps
    RR_loop_count(i,pseudo_reps);
    data_pr = crt_k + complex(randn(size(crt_k)),randn(size(crt_k)));
    data_pr = squeeze(mean(data_pr,av_dim));
    x = ismrm_transform_kspace_to_image(data_pr,[1 2]);
    
    csm = ismrm_estimate_csm_walsh(x);
    ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(ncoils));
    
    img_pr(:,:,i) = abs(sum(x .* ccm_roemer_optimal, 3));
end

% img_pr_orig = img_pr;
img_pr = img_pr( size(crt_k,1)/4:size(crt_k,1)/4 +size(crt_k,1)/2 -1, :,:);

g = std(abs(img_pr + max(abs(img_pr(:)))),[],3); %Measure variation, but add offset to create "high snr condition"
g(g < eps) = 1;
snr = mean(img_pr,3)./g;
%
% figure,
% subplot(1,2,1); imshow(snr,[0 50]); colorbar; colormap('parula')
% subplot(1,2,2); imshow(std(img_pr,[],3),[]); colorbar; colormap('parula')
%
% SNR_s.img_pr = img_pr;
% SNR_s.snr = snr;
% SNR_s.gr = g;

end