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

function [img_s,kspace] = recon_cartesian(file, nfile, SNR_flag)
%% Read data file
raw_data= h5read(file, '/dataset/data');
ismrmrd_s = read_h5_header(file); disp(' ');disp('### Protocol Name ###');disp(ismrmrd_s.measurementInformation.protocolName);disp(' ');
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

%% Noise checks
disp(['Dependency ID: ' ismrmrd_s.measurementInformation.measurementDependency.measurementID]);
asi_names = fieldnames(ismrmrd_s.acquisitionSystemInformation);
n_iRD = read_h5_header(nfile);
nasi_names = fieldnames(n_iRD.acquisitionSystemInformation);

solid_int = 0;
coil_label = cell(channels,4);
coil_flag = zeros(channels,2);
for i = 1:length(asi_names)
    if regexp(asi_names{i}, 'coilLabel')
%         [asi_names{i} ' ' nasi_names{i}]
        solid_int = solid_int + 1;
        coil_label{solid_int,1} = ismrmrd_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilNumber;
        coil_label{solid_int,2} = ismrmrd_s.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilName;
        
        coil_label{solid_int,3} = n_iRD.acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilNumber;
        coil_label{solid_int,4} = n_iRD.acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilName;
        
%         coil_label{solid_int,5} = coil_label{solid_int,1} == coil_label{solid_int,3};
%         coil_label{solid_int,6} = regexp(coil_label{solid_int,2}, coil_label{solid_int,4});
         coil_flag(solid_int,1) = coil_label{solid_int,1} == coil_label{solid_int,3};
         coil_flag(solid_int,2) = regexp(coil_label{solid_int,2}, coil_label{solid_int,4});
    end
end

if sum(coil_flag(:,1)) ~= channels
    disp('Coil order mismatch!');
elseif sum(coil_flag(:,2)) ~= channels
    disp('Coil name mismatch!');
end


%% Calculate Noise Decorrelation Matrix
if nargin > 1
    if isempty(nfile)
        disp('No noise adjustment: ');
        dmtx = diag(ones(1,channels));
    else
        disp('Required noise ID: ');
        disp(ismrmrd_s.measurementInformation.measurementDependency.measurementID);
        dmtx = ismrm_dmtx_RR(nfile, dt);
    end
else
    disp('No noise adjustment: ');
    dmtx = diag(ones(1,channels));
end
disp(' ');

%% Recon

if(length(unique(double(raw_data.head.number_of_samples))) > 1)
    % non-uniform sample number (grappa reference lines?)
    error('GRAPPA recon needed')
    
else

    kspace = zeros([ismrmrd_s.encoding.encodedSpace.matrixSize.x pe1 pe2 averages slices contrasts phases reps sets channels]);
    % disp(['Kspace dims: ' num2str(size(kspace))])

    for ii = 1:length(raw_data.data)
                
        d1 = raw_data.data{ii};
        d2 = complex(d1(1:2:end), d1(2:2:end));
        d3 = reshape(d2, samples, channels); %  RE & IM (2)
        
        if nargin > 1
            d3 = ismrm_apply_noise_decorrelation_mtx(d3, dmtx);
        end
        
        kspace(echo_vec,...
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
end

%% SNR
if nargin < 3
    SNR_flag = 0;
end

if SNR_flag

    snr = zeros(ismrmrd_s.encoding.reconSpace.matrixSize.x, ismrmrd_s.encoding.reconSpace.matrixSize.y, slices);
    
    for iSlice = 1:slices
%         snr(:,:,iSlice) = cartesian_pseudoreps(squeeze(kspace(:,:,:,:,iSlice,:,:,:,:,:)));
          snr(:,:,iSlice) = cartesian_pseudoreps(squeeze(kspace(:,:,:,:,iSlice,1,1,1,1,:)));
    end
end


%% transform to image-space

kspace = mean(kspace,4);% kspace = squeeze(kspace);cCoil_imgs = ismrm_transform_kspace_to_image(kspace);

cCoil_imgs = zeros(size(kspace));
for slc = 1:slices
    for coc = 1:contrasts
        for phc = 1:phases            
            for repc = 1:reps
                for setc = 1:sets
                    for i = 1:channels
                        
                        if asym_e == 0
                            % FFT 
                            
                            cCoil_imgs(:,:,:,1,slc,coc,phc,repc,setc,i) = ismrm_transform_kspace_to_image(kspace(:,:,:,1,slc,coc,phc,repc,setc,i));
                        
                        else
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

dims = size(cCoil_imgs);
CCM_img =  zeros(dims(1:end-1));
for p2c = 1:pe2
    for slc = 1:slices
        for coc = 1:contrasts
            for phc = 1:phases
                for repc = 1:reps
                    for setc = 1:sets
                        
                           temp = cCoil_imgs(:,:,p2c,1,slc,coc,phc,repc,setc,:);
                           csm = ismrm_estimate_csm_walsh( squeeze( temp ) );
                           ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(channels)); % with pre-whitened
                           CCM_img(:,:,p2c,1,slc,coc,phc,repc,setc) = abs( sum( squeeze( temp ) .* ccm_roemer_optimal, 3) );
                            
                    end
                end
            end
        end
    end
end

% Sum-of-Squares coil combine
% % SOS_img = squeeze(sqrt(sum(cCoil_imgs.*conj(cCoil_imgs),length(size(cCoil_imgs)) )));
% % cCoil_imgs = squeeze(cCoil_imgs);

%% Return variables
kspace = squeeze(kspace);
img = squeeze(CCM_img);

new_dims = size(img); 
montage_RR(reshape(img, [new_dims(1) new_dims(2) prod(new_dims)/(new_dims(1)*new_dims(2))]));

img_s.img = img;
img_s.header = header;

if SNR_flag
    img_s.snr = snr;
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

% estimate csm from data 
ncoils = size(crt_k,ndims(crt_k));
img = ismrm_transform_kspace_to_image(squeeze(mean(crt_k,av_dim)),[1 2]); % montage_RR(abs(img));figure, imshow(sqrt(sum(img.*conj(img),3)),[])
csm = ismrm_estimate_csm_walsh(img); % montage_RR((csm));
% calculate coil combine
ccm_roemer_optimal = ismrm_compute_ccm(csm, eye(ncoils)); % with pre-whitened

img_pr = zeros(size(crt_k,1), size(crt_k,2), pseudo_reps);
for i = 1:pseudo_reps
    data_pr = crt_k + complex(randn(size(crt_k)),randn(size(crt_k)));
    data_pr = squeeze(mean(data_pr,av_dim));
    x = ismrm_transform_kspace_to_image(data_pr,[1 2]);
%     img_pr(:,:,i) = sqrt(sum(x.*conj(x),3)); % NO!
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