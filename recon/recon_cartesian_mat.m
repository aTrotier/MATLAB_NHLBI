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

function [img_s,kspace] = recon_cartesian_mat(data, header, SNR_flag)

%%
raw_data= data(2);
ismrmrd_s = header(2);
disp(' ');disp('### Protocol Name ###');disp(ismrmrd_s.measurementInformation.protocolName);disp(' ');

samples = double(raw_data.head.number_of_samples(1));
% asymmetric echo?
% asym_e = 0; if (samples~=ismrmrd_s.encoding.encodedSpace.matrixSize.x); asym_e=1; end;
% echo_vec = (ismrmrd_s.encoding.encodedSpace.matrixSize.x-samples+1):ismrmrd_s.encoding.encodedSpace.matrixSize.x;
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

ismrmrd_s.im_header = im_header;

%% Print out to window

if (samples < ismrmrd_s.encoding.encodedSpace.matrixSize.x); disp('Asymmetric Echo'); end;
% disp(['BW: ' num2str(dt)])

% IM_R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ];
 
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

%% Calculate Noise Decorrelation Matrix

if length(data)==1
    disp('No noise adjustment: ');
    dmtx = diag(ones(1,channels));
else
    disp('Required noise ID: ');
    disp(ismrmrd_s.measurementInformation.measurementDependency.measurementID);
    
    %% Noise checks
    
        asi_names = fieldnames(header(2).acquisitionSystemInformation);
        nasi_names = fieldnames(header(1).acquisitionSystemInformation);
    
        solid_int = 0;
        coil_label = cell(channels,4);
        coil_flag = zeros(channels,2);
        for i = 1:length(asi_names)
            if regexp(asi_names{i}, 'coilLabel')
                %         [asi_names{i} ' ' nasi_names{i}]
                solid_int = solid_int + 1;
                coil_label{solid_int,1} = ismrmrd_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilNumber;
                coil_label{solid_int,2} = ismrmrd_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilName;
    
                coil_label{solid_int,3} = header(1).acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilNumber;
                coil_label{solid_int,4} = header(1).acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilName;
    
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
    %
    noise_test = data(1);
    %     iRD_s = header(1);
    % disp('Noise ID: ');
    % disp(iRD_s.measurementInformation.measurementID);
    
    Siemens_rBW = ismrmrd_s.acquisitionSystemInformation.relativeReceiverNoiseBandwidth; % Siemens_rBW = 0.79;
    
    n_samples = double(noise_test.head.number_of_samples(1));
    n_channels = double(noise_test.head.active_channels(1));
    % assuming Siemens using 2 averages:
    noise_ind = find(noise_test.head.idx.average==1, 1, 'last');
    
    %% Plot noise STD
    nt2 = zeros(n_samples, noise_ind, n_channels);
    for i = 1:noise_ind
        nt2(:,i,:)=  double(reshape(complex(noise_test.data{i}(1:2:end), noise_test.data{i}(2:2:end)), [n_samples, 1, n_channels ]));
    end
    nt3 = reshape(nt2,size(nt2,1)*size(nt2,2), size(nt2,3));
    n_stats = std(nt3,[],1);
    
    figure, hold on, title('Coil Noise SD (m +/- std)'); xlabel('Coils'); ylabel('Noise SD');
    plot(1:n_channels, n_stats, 'bo');
    plot([ones(n_channels,1).*mean(n_stats)], 'r-');
    plot([ones(n_channels,1).*(mean(n_stats)-std(n_stats)) ones(n_channels,1).*(mean(n_stats)+std(n_stats))], 'r--');
    
    str = {['Mean ' num2str(mean(n_stats))] ['Std  ' num2str(std(n_stats))]};dim = [.2 .5 .3 .3];
    annotation('textbox',dim,'String',str,'FitBoxToText','on')
    
    n_scaling = Siemens_rBW * dt / (noise_test.head.sample_time_us(1)*1e-6);
    dmtx = ismrm_calculate_noise_decorrelation_mtx(nt2, n_scaling ); % figure,imagesc(abs(dmtx));
    
end

disp(' ');

%% Sort data
% navigator extract
% acq_frames = find(single(raw_data.head.number_of_samples) == ismrmrd_s.encoding.encodedSpace.matrixSize.x);
acq_frames = (1:length(raw_data.data))';
    
    kspace = complex(zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 pe2 averages slices contrasts phases reps sets channels],'single'));
    % disp(['Kspace dims: ' num2str(size(kspace))])
    
    for ii = acq_frames'%2:length(raw_data.data) % %
        
        d1 = raw_data.data{ii};
        d2 = complex(d1(1:2:end), d1(2:2:end));
        d3 = reshape(d2, samples, channels); %  RE & IM (2)
        d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
        
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


%% SNR
if nargin < 3
    SNR_flag = 0;
end

if SNR_flag
    if pe2 == 1
        if isempty(SNR_flag)
            pseudoRep_slices = 1:slices;
        end
        if SNR_flag > 0 % specify slice(s);
            pseudoRep_slices = SNR_flag;
            
            if pseudoRep_slices > slices
                pseudoRep_slices = RR_slice_order(round(slices/2));
                warning(['PR slice > #slices, using slice ' num2str(pseudoRep_slices)])
            end
        end
        
        %     snr = zeros(ismrmrd_s.encoding.reconSpace.matrixSize.x, ismrmrd_s.encoding.reconSpace.matrixSize.y, slices);
        snr = zeros(ismrmrd_s.encoding.reconSpace.matrixSize.x, size(kspace,2), slices);
        
        for iSlice = pseudoRep_slices
            %         snr(:,:,iSlice) = cartesian_pseudoreps(squeeze(kspace(:,:,:,:,iSlice,:,:,:,:,:)));
            snr(:,:,iSlice) = cartesian_pseudoreps(squeeze(kspace(:,:,:,:,iSlice,1,1,1,1,:)));
        end
    else
        
        snr = cartesian_3D_pseudoreps(squeeze(kspace));
        
    end
end


%% transform to image-space

kspace = mean(kspace,4);
% no GRAPPA for the mo
if pe2 == 1
    % 2D
    cCoil_imgs = ismrm_transform_kspace_to_image(kspace, [1 2]);    
else
    % 3D
    cCoil_imgs = ismrm_transform_kspace_to_image(kspace, [1 2 3]);
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

%% Roemer coil combination "slice-by-slice"
% if ismrmrd_s.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 == 1
dims = size(cCoil_imgs);
CCM_img =  zeros(dims(1:end-1));

if pe2 ==1
    % for par = 1:pe2
    par = 1;
    for slc = 1:slices
        for coc = 1:contrasts
            for phc = 1:phases
                for repc = 1:reps
                    for setc = 1:sets
                        
                        temp = cCoil_imgs(:,:,par,1,slc,coc,phc,repc,setc,:);
                        csm = ismrm_estimate_csm_walsh( squeeze( temp ) );
                        ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
                        CCM_img(:,:,par,1,slc,coc,phc,repc,setc) = abs( sum( squeeze( temp ) .* ccm_roemer_optimal, 3) );
                      end
                end
            end
        end
    end
    % end
else
    
    for coc = 1:contrasts
        for phc = 1:phases
            for repc = 1:reps
                for setc = 1:sets
                    
                    temp = cCoil_imgs(:,:,:,1,1,coc,phc,repc,setc,:);
                    csm = ismrm_estimate_csm_walsh_3D( squeeze( temp ) );
                    ccm_roemer_optimal = ismrm_compute_ccm_3D(csm, eye(channels)); % with pre-whitened
                    CCM_img(:,:,:,1,1,coc,phc,repc,setc) = abs( sum( squeeze( temp ) .* ccm_roemer_optimal, 4) );
                    
                end
            end
        end
    end
    
end
% else
%     CCM_img = cCoil_imgs;
% end

% Sum-of-Squares coil combine
% SOS_img = squeeze(sqrt(sum(cCoil_imgs.*conj(cCoil_imgs),length(size(cCoil_imgs)) )));
% CCM_img = squeeze(SOS_img);


%% Return variables
kspace = squeeze(kspace);
img = squeeze(CCM_img);

new_dims = size(img);
montage_RR(reshape(img, [new_dims(1) new_dims(2) prod(new_dims)/(new_dims(1)*new_dims(2))]));

% imgc = squeeze(CCM_imgc);
% montage_RR(reshape(imgc, [new_dims(1) new_dims(2) prod(new_dims)/(new_dims(1)*new_dims(2))]),[-pi pi]);

img_s.img = img;
img_s.header = ismrmrd_s;

if SNR_flag
    img_s.snr = snr;
end

end

function snr = cartesian_3D_pseudoreps(crt_k)
% SNR_s = cartesian_psuedoreps(crt_k, slice)
% Expects pre-whitened 3D k-space data [RO,PE,PE2,AV,COILS]
% This can be expanded to generic dimensions, but will require header info
% on dimension - better used as a subfunction within recon_cartesian
% R Ramasawmy NHLBI 2019

if length(size(crt_k)) == 4
    crt_k = reshape(crt_k,[size(crt_k,1) size(crt_k,2) size(crt_k,3) 1 size(crt_k,4)]);
end
% size(crt_k)
av_dim = 4;

pseudo_reps = 100;
disp(['Running ' num2str(pseudo_reps) ' pseudo-reps']);

ncoils = size(crt_k,ndims(crt_k));
img_pr = zeros(size(crt_k,1), size(crt_k,2),size(crt_k,3), pseudo_reps, 'single');

parfor i = 1:pseudo_reps
    RR_loop_count(i,pseudo_reps);
    data_pr = crt_k + complex(randn(size(crt_k)),randn(size(crt_k)));
    data_pr = squeeze(mean(data_pr,av_dim));
    x = ismrm_transform_kspace_to_image(data_pr,[1 2 3]);
    
    csm = ismrm_estimate_csm_walsh_3D(x);
    ccm_roemer_optimal = ismrm_compute_ccm_3D(csm, eye(ncoils));
    img_pr(:,:,:,i) = abs(sum(x .* ccm_roemer_optimal, 4));
end

% img_pr_orig = img_pr;
img_pr = img_pr( size(crt_k,1)/4:size(crt_k,1)/4 +size(crt_k,1)/2 -1, :,:,:);

g = std(abs(img_pr + max(abs(img_pr(:)))),[],4); %Measure variation, but add offset to create "high snr condition"
g(g < eps) = 1;
snr = mean(img_pr,4)./g;
%
% figure,
% subplot(1,2,1); imshow(snr,[0 50]); colorbar; colormap('parula')
% subplot(1,2,2); imshow(std(img_pr,[1 2],3),[]); colorbar; colormap('parula')
%
% SNR_s.img_pr = img_pr;
% SNR_s.snr = snr;
% SNR_s.gr = g;

end

function ccm = ismrm_compute_ccm_3D(csm, noise_matrix)
%
%  ccm = ismrm_compute_ccm(csm, noise_matrix)
%
%  Computes noise-optimal channel combination maps from  coil sensitivity
%  maps and a noise covariance matrix.
%
%  The ccm can be applied to channel-by-channel images as
%
%  im_composite = sum(ccm .* im_channel_by_channel, 3);
%
%  INPUT:
%    - csm [x,y, coil]          : coil sensitivity maps
%    - noise_matrix [coil,coil] : noise covariance matrix
%
%  OUTPUT:
%    - ccm [x,y, coil]          : channel combination maps
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
 % >>>>>>
    
    nx = size(csm, 1);
    ny = size(csm, 2);
    nz = size(csm, 3);
    nc = size(csm, 4);
    
    if( nargin <2 || isempty(noise_matrix) )
        noise_matrix = eye(nc);
    end
    
    csm_matrix = reshape(csm, [nx*ny*nz nc]);
    
    relative_ccm = conj(csm_matrix) * pinv(noise_matrix);
    
    scale_correction = abs(sum(relative_ccm .* csm_matrix, 2));
    
    nonzero_ind = scale_correction>0;
    
    ccm = zeros(size(csm_matrix));
    ccm(nonzero_ind, :) = relative_ccm(nonzero_ind,:) ./ repmat(scale_correction(nonzero_ind), [1 nc]);
    
    ccm = reshape(ccm, [nx ny nz nc]);
    
    % <<<<<<<<<<<<<<<<<<<
end

function [csm] = ismrm_estimate_csm_walsh_3D(img, smoothing)
%
%   [csm] = ismrm_estimate_csm_walsh(img)
%
%   Estimates relative coil sensitivity maps from a set of coil images
%   using the eigenvector method described by Walsh et al. (Magn Reson Med
%   2000;43:682-90.)
%
%   INPUT:
%     - img     [x,y,coil]   : Coil images
%     - smooth  scalar       : Smoothing block size (defaults to 5)
%
%   OUTPUT:
%
%     - csm     [x,y,coil    : Relative coil sensitivity maps
%
%
%   Code is based on an original implementation by Peter Kellman, NHLBI,
%   NIH
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%

if nargin < 2,
    smoothing = 5;
end

ncoils = size(img,ndims(img));

% normalize by root sum of squares magnitude
mag = sqrt(sum(img .* conj(img),ndims(img)));
s_raw=img ./ repmat(mag + eps,[1 1 1 ncoils]); clear mag; 


% compute sample correlation estimates at each pixel location
Rs=(ismrm_correlation_matrix_3D(s_raw));


% apply spatial smoothing to sample correlation estimates (NxN convolution)
if smoothing>1,
	h_smooth=ones(smoothing)/(smoothing^2); % uniform smoothing kernel
	for m=1:ncoils
		for n=1:ncoils
		    Rs(:,:,:,m,n)=convn(Rs(:,:,:,m,n),h_smooth,'same');
		end
	end
end


% compute dominant eigenvectors of sample correlation matrices
[csm,lambda]=ismrm_eig_power_3D(Rs); % using power method


end 


%Utility functions provided by Peter Kellman, NIH.
function [Rs]=ismrm_correlation_matrix_3D(s)
% function [Rs]=correlation_matrix(s);
%
% function correlation_matrix calculates the sample correlation matrix (Rs) for
% each pixel of a multi-coil image s(y,x,coil)
%
% input:
%    s   complex multi-coil image s(y,x,coil)
% output:
%    Rs  complex sample correlation matrices, Rs(y,x,coil,coil)

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

[rows,cols,slices,ncoils]=size(s);
Rs=zeros(rows,cols,slices,ncoils,ncoils, 'single'); % initialize sample correlation matrix to zero
for i=1:ncoils
    for j=1:i-1
		Rs(:,:,:,i,j)=s(:,:,:,i).*conj(s(:,:,:,j));
		Rs(:,:,:,j,i)=conj(Rs(:,:,:,i,j)); % using conjugate symmetry of Rs
    end
	Rs(:,:,:,i,i)=s(:,:,:,i).*conj(s(:,:,:,i));
end

end

function [v,d]=ismrm_eig_power_3D(R)
% function [v,d]=eig_power(R);
%
% vectorized method for calculating the dominant eigenvector based on
% power method. Input, R, is an image of sample correlation matrices
% where: R(y,x,:,:) are sample correlation matrices (ncoil x ncoil) for each pixel
%
% v is the dominant eigenvector
% d is the dominant (maximum) eigenvalue

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

rows=size(R,1);cols=size(R,2);slices=size(R,3);ncoils=size(R,4);
N_iterations=2;
v=ones(rows,cols,slices,ncoils, 'single'); % initialize e.v.

d=zeros(rows,cols,slices, 'single');
for i=1:N_iterations
    v=squeeze(sum(R.*repmat(v,[1 1 1 1 ncoils]),4));
	d=ismrm_rss(v);
    d( d <= eps) = eps;
	v=v./repmat(d,[1 1 1 ncoils]);
end

p1=angle(conj(v(:,:,:,1)));
% (optionally) normalize output to coil 1 phase
v=v.*repmat(exp(sqrt(-1)*p1),[1 1 1 ncoils]);
v=conj(v);

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
img_pr = zeros(size(crt_k,1), size(crt_k,2), pseudo_reps, 'single');

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
% subplot(1,2,2); imshow(std(img_pr,[1 2],3),[]); colorbar; colormap('parula')
%
% SNR_s.img_pr = img_pr;
% SNR_s.snr = snr;
% SNR_s.gr = g;

end