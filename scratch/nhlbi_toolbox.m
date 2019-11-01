%  ????? NHLBI ??????
% ????? TOOLBOX ??????

% called by nhlbi_toolbox script
% (This is real nasty "class"-ing.)

% -- Search for local functions --
function fh = nhlbi_toolbox
fh = localfunctions;
end

% -- edit call --
function edit_toolbox
edit(mfilename)
end

% -------- TOOLS (w/ dependent functions) ---------------

function ismrmrd_xml = h5read_xml(h5file)
xml_string = h5read(h5file, '/dataset/xml');
ismrmrd_xml = xml2hdr(xml_string{1});
end

% -------- TOOLS (w/o dependent functions) ---------------
function check_noise_dependency(iRD_s_data, iRD_s_noise)
asi_names = fieldnames(iRD_s_data.acquisitionSystemInformation);
nasi_names = fieldnames(iRD_s_noise.acquisitionSystemInformation);
channels = iRD_s_noise.acquisitionSystemInformation.receiverChannels;

solid_int = 0;
coil_label = cell(channels,4);
coil_flag = zeros(channels,2);
for i = 1:length(asi_names)
    if regexp(asi_names{i}, 'coilLabel')
        %         [asi_names{i} ' ' nasi_names{i}]
        solid_int = solid_int + 1;
        coil_label{solid_int,1} = iRD_s_data.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilNumber;
        coil_label{solid_int,2} = iRD_s_data.acquisitionSystemInformation.(matlab.lang.makeValidName(asi_names{i})).coilName;
        
        coil_label{solid_int,3} = iRD_s_noise.acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilNumber;
        coil_label{solid_int,4} = iRD_s_noise.acquisitionSystemInformation.(matlab.lang.makeValidName(nasi_names{i})).coilName;
        
        %         coil_label{solid_int,5} = coil_label{solid_int,1} == coil_label{solid_int,3};
        %         coil_label{solid_int,6} = regexp(coil_label{solid_int,2}, coil_label{solid_int,4});
        coil_flag(solid_int,1) = coil_label{solid_int,1} == coil_label{solid_int,3};
        coil_flag(solid_int,2) = strcmp(coil_label{solid_int,2}, coil_label{solid_int,4});
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

end

function data_pca = coil_pca(data, coils_out)
% function data_pca = coil_pca(data, coils_out)
% 2D data
reshape_data = 0;
coils_in = size(data,ndims(data));

if ndims(data) > 2
    reshape_data = 1;
    warning('data dims > 2, reshaping - assuming coils in last dim');
    data_in_dims = size(data);
    
    data = reshape(data, [prod(data_in_dims)/coils_in coils_in]);
end


[U,S,V] = svd(data,'econ');

%% Rank Contribution [ academic : if you want to see which rank of motion data]
debug = 0;
if debug
    
    figure, subplot(1,2,1); plot(diag(S))
    
    rank_matrix = zeros(size(U));
    for i = 1:length(S)
        rank_matrix(:,i) = U(:,i)*S(i,i)*V(i,i)';
    end
    subplot(1,2,2); colorful_plots(abs(rank_matrix))
    
end

%% PCA

data_pca = U(:,1:coils_out)*S(1:coils_out,1:coils_out)*V(1:coils_out,1:coils_out)';

if reshape_data
    data_in_dims(end) = coils_out;
    data_pca = reshape(data_pca, data_in_dims);
end

end

function iRD_loop_counter(i_vec, iRD_vec)
% print out current | PE2, slice, contrast, phase, repetition, set |
% assumed average will be handled upstream, and 2D operations.

print_line = ['\n ' repmat('#', [1 80]) '\n PE2: %d/%d | slice: %d/%d | contrast: %d/%d | phase: %d/%d | repetition: %d/%d | set: %d/%d \n ' repmat('#', [1 80]) '\n'];

if (sum(i_vec) > length(i_vec))
    % bodge: print out the line, and remove it twice.. (will get messy when
    % the indices change order, i.e. 9 -> 10 .. but seem ok..)
    a = fprintf(print_line, i_vec(1),iRD_vec(1),i_vec(2),iRD_vec(2),i_vec(3),iRD_vec(3),i_vec(4),iRD_vec(4),i_vec(5),iRD_vec(5),i_vec(6),iRD_vec(6));
    fprintf([repmat('\b', [1 2*a])]);
end

a = fprintf(print_line, i_vec(1),iRD_vec(1),i_vec(2),iRD_vec(2),i_vec(3),iRD_vec(3),i_vec(4),iRD_vec(4),i_vec(5),iRD_vec(5),i_vec(6),iRD_vec(6));

end

function plot_experiment(raw_data)

hl = [];
temp = fieldnames(raw_data.head.idx);
figure('Name', 'Experiment header');
for i = 1:9 % dont plot user vector
    h.(matlab.lang.makeValidName(['x' num2str(i)])) = subplot(3, 3, i);
    plot( 1+double(raw_data.head.idx.(matlab.lang.makeValidName(temp{i}))) ); title(temp{i}, 'Interpreter', 'none')
    hl = [hl, h.(matlab.lang.makeValidName(['x' num2str(i)]))];
end
linkaxes(hl, 'x');

end

function [v_angle] = vector_angle(u,v)
% % https://www.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab
% v_angle = atan2d(norm(cross(u,v)),dot(u,v));

% % "directional"
c = cross(u,v);
v_angle = sign(c(3))*180/pi*atan2(norm(c),dot(u,v));

end

function [psn_s, n_MRD_s] = prescan_normalize(nfile, i_rd,  nhlbi_toolbox)
% beta_function, inputs and outputs will change
% [psn_lr, n_MRD_s] = prescan_normalize(nfile, i_rd,  nhlbi_toolbox)
%
% NHLBI code for Siemens implementation VE11C+
% R Ramasawmy, July 2019
% prescan normalise is 3D sagittal acquisition, foot-to-head and anterior-to-posterior (head-first + supine)
% however the coordinates in the data header think its axial.. so a sagital
% (F-H, A-P) will be assumed for the following reconstruction.

noise_test = h5read(nfile,'/dataset/data');
noise_test_head = noise_test.head;
n_MRD_s = nhlbi_toolbox.h5read_xml(nfile);

disp(['Required Sens Map: ' i_rd.measurementInformation.measurementDependency.measurementID ', Noise ID: ' n_MRD_s.measurementInformation.measurementID]);
nhlbi_toolbox.check_noise_dependency(i_rd,n_MRD_s);

% find noise-scan: assuming Siemens using 2 averages:
noise_ind = find(noise_test.head.idx.average==1, 1, 'last');

if length(noise_test.data) == noise_ind(end)
    % no prescan normalise present
    psn_s.present = 0;
    % dummy return?
    %     B = [repmat(0, [1 9]) 1];
    %     surfit_s = @(B,XYZ) B(1)*XYZ(:,:,:,1).^2 + B(2)*XYZ(:,:,:,2).^2 + B(3)*XYZ(:,:,:,3).^2 + ... [second order] terms
    %         B(4)*XYZ(:,:,:,1).*XYZ(:,:,:,2) + B(5)*XYZ(:,:,:,1).*XYZ(:,:,:,3) +  B(6)*XYZ(:,:,:,2).*XYZ(:,:,:,3) + ... [cross] terms
    %         B(7)*XYZ(:,:,:,1) + B(8)*XYZ(:,:,:,2) + B(9)*XYZ(:,:,:,3) + B(10); % [1inear and 0th terms] terms
    %
    %     psn_s.model_3D          = surfit_s;
    %     psn_s.model_3D_coefs    = B;
    
else
    psn_s.present = 1;
    
    noise_ind = noise_ind+1:length(noise_test.data);
    
    n_samples = double(noise_test.head.number_of_samples(noise_ind(1)));
    n_channels = double(noise_test.head.active_channels(noise_ind));
    % pe_tab_1 = 1+(noise_test.head.idx.kspace_encode_step_1(noise_ind));
    % pe_tab_2 = 1+(noise_test.head.idx.kspace_encode_step_2(noise_ind));
    % nt3_AC = zeros(n_samples, max(pe_tab_1), max(pe_tab_2), n_channels(1)); size(nt3_AC)
    % nt3_BC = zeros(n_samples, max(pe_tab_1), max(pe_tab_2), n_channels(2));  size(nt3_BC)
    
    nt3_AC = zeros(n_samples, n_MRD_s.encoding.encodedSpace.matrixSize.y, n_MRD_s.encoding.encodedSpace.matrixSize.z, n_channels(1)); %size(nt3_AC)
    nt3_BC = zeros(n_samples, n_MRD_s.encoding.encodedSpace.matrixSize.y, n_MRD_s.encoding.encodedSpace.matrixSize.z, n_channels(2)); %size(nt3_BC)
    
    for i = noise_ind
        if noise_test.head.active_channels(i) == 2
            nt3_BC(:,noise_test.head.idx.kspace_encode_step_1(i)+1,noise_test.head.idx.kspace_encode_step_2(i)+1,:)= double(reshape(complex(noise_test.data{i}(1:2:end), noise_test.data{i}(2:2:end)), [n_samples 1 1 2]));
        else
            nt3_AC(:,noise_test.head.idx.kspace_encode_step_1(i)+1,noise_test.head.idx.kspace_encode_step_2(i)+1,:)= double(reshape(complex(noise_test.data{i}(1:2:end), noise_test.data{i}(2:2:end)), [n_samples 1 1 noise_test.head.active_channels(i)]));
        end
    end
    
    img_AC = ismrm_transform_kspace_to_image(nt3_AC);
    img_BC = ismrm_transform_kspace_to_image(nt3_BC);
    img_AC = sqrt(sum(img_AC.*conj(img_AC),4));
    img_BC = sqrt(sum(img_BC.*conj(img_BC),4));
    
    % montage_RR(img_AC);montage_RR(img_BC);
    % montage_RR(img_AC,[0 1]);montage_RR(img_BC,[0 1]);
    scaled_coil_image = img_AC./img_BC;
    mask = img_BC > std(img_BC(:))+mean(img_BC(:)); %montage_RR(mask) % montage_RR(bwareaopen(mask,100))  % montage_RR(imerode(bwareaopen(mask,100), strel('sphere', 1)))
    mask = imerode(bwareaopen(mask,100), strel('sphere', 1)); % shrink the mask to avoid partial volume effects?
    montage_RR(scaled_coil_image.*mask,[0 mean(scaled_coil_image(mask))]);
    
    psn_lr = scaled_coil_image.*mask;
    
    psn_fov = [n_MRD_s.encoding.encodedSpace.fieldOfView_mm.x n_MRD_s.encoding.encodedSpace.fieldOfView_mm.y n_MRD_s.encoding.encodedSpace.fieldOfView_mm.z];
    psn_res = psn_fov./size(psn_lr);
    
    % == export ==
    psn_s.psn_masked = psn_lr;
    psn_s.psn_mask = mask;
    psn_s.psn_image = scaled_coil_image;
    psn_s.fov = psn_fov;
    psn_s.res = psn_res;
    psn_s.R = [0 0 -1; 0 1 0; 1 0 0]; % ?
    
    %% 3D surface fitting (doesn't work as well..?)
    % close all
    matrixSize = size(scaled_coil_image);
    
    % x = 1:matrixSize(1); y = 1:matrixSize(2);  z = 1:matrixSize(3);
    % physical space coordinates for mapping to images
    z = linspace([0.5]*(psn_fov(1) -0.5*psn_res(1)),-[0.5]*(psn_fov(1) -0.5*psn_res(1)), matrixSize(1));
    y = linspace([0.5]*(psn_fov(2) -0.5*psn_res(2)),-[0.5]*(psn_fov(2) -0.5*psn_res(2)), matrixSize(2));
    x = linspace([-0.5]*(psn_fov(3) -0.5*psn_res(3)),[0.5]*(psn_fov(3) -0.5*psn_res(3)), matrixSize(3));
    
    [X, Y, Z] = meshgrid(y, z, x);
    XYZ = zeros(length(z), length(y), length(x), 3);
    XYZ(:,:,:,1) = X; XYZ(:,:,:,2) = Y; XYZ(:,:,:,3) = Z;
    a = zeros(length(find(mask)),3);
    temp=XYZ(:,:,:,1); a(:,1) = temp(find(mask));
    temp=XYZ(:,:,:,2); a(:,2) = temp(find(mask));
    temp=XYZ(:,:,:,3); a(:,3) = temp(find(mask));
    
    opts= optimset('Display', 'none');
    
    % Three dimensional squared polynomial:
    %   ax2 + by2 + cz2 + dxy + exz + fyz + gx + hy + iz + j  =  0
    
    surfit_v = @(B,a) B(1)*a(:,1).^2 + B(2)*a(:,2).^2 + B(3)*a(:,3).^2 + ... [2nd order] terms
        B(4)*a(:,1).*a(:,2) + B(5)*a(:,1).*a(:,3) +  B(6)*a(:,2).*a(:,3) + ... [cross] terms
        B(7)*a(:,1) + B(8)*a(:,2) + B(9)*a(:,3) + B(10); % [1st and 0th terms] terms

    pn_factor = max(psn_lr(:));
    
    initial_guess = [repmat(0, [1 9]) 0.5*pn_factor];
    fit_minimum = [repmat(-10, [1 9]) 0]*pn_factor;
    fit_maximum = [repmat(10, [1 9]) 10]*pn_factor;
    
    %    lsqcurvefit(func handle, initial guess, ROI, original data,                 fit minimum , fit maximum, opts)
    B =  lsqcurvefit(surfit_v,    initial_guess,   a, scaled_coil_image(find(mask)), fit_minimum , fit_maximum, opts);
    
    surfit_s = @(B,XYZ) B(1)*XYZ(:,:,:,1).^2 + B(2)*XYZ(:,:,:,2).^2 + B(3)*XYZ(:,:,:,3).^2 + ... [2nd order] terms
        B(4)*XYZ(:,:,:,1).*XYZ(:,:,:,2) + B(5)*XYZ(:,:,:,1).*XYZ(:,:,:,3) +  B(6)*XYZ(:,:,:,2).*XYZ(:,:,:,3) + ... [cross] terms
        B(7)*XYZ(:,:,:,1) + B(8)*XYZ(:,:,:,2) + B(9)*XYZ(:,:,:,3) + B(10); % [1inear and 0th terms] terms
    
    coil_sens = surfit_s(B,XYZ);
    
    mtemp = montage_RR(mask); close;
    
    montage_RR(scaled_coil_image,[0 pn_factor]); hold on; contour(mtemp,'w');
    montage_RR(coil_sens,[0 pn_factor]); hold on; contour(mtemp,'w');
    
    % == export ==
    psn_s.model_3D          = surfit_s;
    psn_s.model_3D_coefs    = B;
    psn_s.model_3D_cords    = XYZ;
    
    %% Three dimensional cubic polynomial: <<?>>
    %  ax3 + by3 + cz3, ... cubic
    %  dx2y + ex2z + fy2x + gy2z + hz2x + iz2y + jxyz ... cubic cross
    %  kx2 + ly2 + mz2 ... squared
    %  nxy + oyz + pxz ... squared cross
    %  qx + ry + sz + t  =  0 linear and 0th
    
    surfit_v = @(B,a) B(1)*a(:,1).^3 + B(2)*a(:,2).^3 + B(3)*a(:,3).^3 + ... [3rd order] terms
        B(4)*a(:,1).^2.*a(:,2) + B(5)*a(:,1).^2.*a(:,3) + ... [3rd order] cross terms //x2
        B(6)*a(:,2).^2.*a(:,1) + B(7)*a(:,2).^2.*a(:,3) + ... [3rd order] cross terms //y2
        B(8)*a(:,3).^2.*a(:,1) + B(9)*a(:,3).^2.*a(:,2) + B(10)*a(:,1).*a(:,2).*a(:,3) + ... [3rd order] cross terms //z2
        B(11)*a(:,1).^2 + B(12)*a(:,2).^2 + B(13)*a(:,3).^2 + ... [2nd order] terms
        B(14)*a(:,1).*a(:,2)  + B(15)*a(:,2).*a(:,3)  + B(16)*a(:,1).*a(:,3) + ... [2nd order] cross terms
        B(17)*a(:,1) + B(18)*a(:,2) + B(19)*a(:,3) + B(20); % [1st and 0th terms] terms
    
    surfit_s = @(B,XYZ) B(1)*XYZ(:,:,:,1).^3 + B(2)*XYZ(:,:,:,2).^3 + B(3)*XYZ(:,:,:,3).^3 + ... [3rd order] terms
        B(4)*XYZ(:,:,:,1).^2.*XYZ(:,:,:,2) + B(5)*XYZ(:,:,:,1).^2.*XYZ(:,:,:,3) + ... [3rd order] cross terms //x2
        B(6)*XYZ(:,:,:,2).^2.*XYZ(:,:,:,1) + B(7)*XYZ(:,:,:,2).^2.*XYZ(:,:,:,3) + ... [3rd order] cross terms //y2
        B(8)*XYZ(:,:,:,3).^2.*XYZ(:,:,:,1) + B(9)*XYZ(:,:,:,3).^2.*XYZ(:,:,:,2) + B(10)*XYZ(:,:,:,1).*XYZ(:,:,:,2).*XYZ(:,:,:,3) + ... [3rd order] cross terms //z2
        B(11)*XYZ(:,:,:,1).^2 + B(12)*XYZ(:,:,:,2).^2 + B(13)*XYZ(:,:,:,3).^2 + ... [2nd order] terms
        B(14)*XYZ(:,:,:,1).*XYZ(:,:,:,2)  + B(15)*XYZ(:,:,:,2).*XYZ(:,:,:,3)  + B(16)*XYZ(:,:,:,1).*XYZ(:,:,:,3) + ... [2nd order] cross terms
        B(17)*XYZ(:,:,:,1) + B(18)*XYZ(:,:,:,2) + B(19)*XYZ(:,:,:,3) + B(20); % [1st and 0th terms] terms
    %
    %  %   test3 = surfit_s(1:20, XYZ);
    %  %   implay_RR(test3)
    
    initial_guess = [repmat(0, [1 19]) 0.5*pn_factor];
    fit_minimum = [repmat(-10, [1 19]) 0]*pn_factor;
    fit_maximum = [repmat(10, [1 19]) 10]*pn_factor;
    
    %    lsqcurvefit(func handle, initial guess, ROI, original data,                 fit minimum , fit maximum, opts)
    B =  lsqcurvefit(surfit_v,    initial_guess,   a, scaled_coil_image(find(mask)), fit_minimum , fit_maximum, opts);
       
    coil_sens = surfit_s(B,XYZ);

    figure, imagesc([scaled_coil_image(:,:,size(scaled_coil_image,3)/2) coil_sens(:,:,size(scaled_coil_image,3)/2)]) 
    hold on; contour([mask(:,:,size(scaled_coil_image,3)/2) mask(:,:,size(scaled_coil_image,3)/2)],'w');
    
    % == export ==
    psn_s.model3_3D          = surfit_s;
    psn_s.model3_3D_coefs    = B;clc
    psn_s.model3_3D_cords    = XYZ;
    
    %% Three dimensional quartic polynomial <<?>>
% %     %  c1_x4 + c2_y4 + c3_z4, ... quartic
% %     %  c4_x3y + c5_x3z + c6_x2yz + c7_x2y2 ... quartic cross /x
% %     %  c8_y3x + c9_y3z + c10_y2xz + c11_y2z2 ... quartic cross /y
% %     %  c12_z3x + c13_z3y + c14_z2xy + c15_z2x2... quartic cross /z
% %     %  c16_x3 + c17_y3 + c18_z3, ... cubic
% %     %  c19_x2y + c20_x2z +  ... cubic cross x
% %     %  c21_y2x + c22_y2z +  ... cubic cross y
% %     %  c23_z2x + c24_z2y + c25_xyz ... cubic cross z
% %     %  c26_x2 + c27_y2 + c28_z2 ... squared
% %     %  c29_xy + c30_yz + c31_xz ... squared cross
% %     %  c32_x + c33_y + c34_z + c35_  =  0 1st and 0th
% %     
% %      surfit_v = @(B,a) B(1)*a(:,1).^4 + B(2)*a(:,2).^4 + B(3)*a(:,3).^4 + ... [4th order] terms
% %         B(4)* a(:,1).^3.*a(:,2) + B(5)* a(:,1).^3.*a(:,3) + B(6)* a(:,1).^2.*a(:,2).*a(:,3) + B(7)* a(:,1).^2.*a(:,2).^2 + ... [4th order] cross terms //x
% %         B(8)* a(:,2).^3.*a(:,1) + B(9)* a(:,2).^3.*a(:,3) + B(10)*a(:,2).^2.*a(:,1).*a(:,3) + B(11)*a(:,2).^2.*a(:,3).^2 + ... [4th order] cross terms //y
% %         B(12)*a(:,3).^3.*a(:,1) + B(13)*a(:,3).^3.*a(:,2) + B(14)*a(:,3).^2.*a(:,1).*a(:,2) + B(15)*a(:,3).^2.*a(:,1).^2 + ... [4th order] cross terms //z
% %         B(16)*a(:,1).^3 + B(17)*a(:,2).^3 + B(18)*a(:,3).^3 + ... [3rd order] terms    
% %         B(19)*a(:,1).^2.*a(:,2) + B(20)*a(:,1).^2.*a(:,3) + ... [3rd order] cross terms //x
% %         B(21)*a(:,2).^2.*a(:,1) + B(22)*a(:,2).^2.*a(:,3) + ... [3rd order] cross terms //y
% %         B(23)*a(:,3).^2.*a(:,1) + B(24)*a(:,3).^2.*a(:,2) + B(25)*a(:,1).*a(:,2).*a(:,3) + ... [3rd order] cross terms //z
% %         B(26)*a(:,1).^2 + B(27)*a(:,2).^2 + B(28)*a(:,3).^2 + ... [2nd order] terms
% %         B(29)*a(:,1).*a(:,2)  + B(30)*a(:,2).*a(:,3)  + B(31)*a(:,1).*a(:,3) + ... [2nd order] cross terms
% %         B(32)*a(:,1) + B(33)*a(:,2) + B(34)*a(:,3) + B(35);  % [1st and 0th terms] terms
% %     
% %     surfit_s = @(B,XYZ) B(1)*XYZ(:,:,:,1).^4 + B(2)*XYZ(:,:,:,2).^4 + B(3)*XYZ(:,:,:,3).^4 + ... [4th order] terms
% %         B(4)* XYZ(:,:,:,1).^3.*XYZ(:,:,:,2) + B(5)* XYZ(:,:,:,1).^3.*XYZ(:,:,:,3) + B(6)* XYZ(:,:,:,1).^2.*XYZ(:,:,:,2).*XYZ(:,:,:,3) + B(7)* XYZ(:,:,:,1).^2.*XYZ(:,:,:,2).^2 + ... [4th order] cross terms //x
% %         B(8)* XYZ(:,:,:,2).^3.*XYZ(:,:,:,1) + B(9)* XYZ(:,:,:,2).^3.*XYZ(:,:,:,3) + B(10)*XYZ(:,:,:,2).^2.*XYZ(:,:,:,1).*XYZ(:,:,:,3) + B(11)*XYZ(:,:,:,2).^2.*XYZ(:,:,:,3).^2 + ... [4th order] cross terms //y
% %         B(12)*XYZ(:,:,:,3).^3.*XYZ(:,:,:,1) + B(13)*XYZ(:,:,:,3).^3.*XYZ(:,:,:,2) + B(14)*XYZ(:,:,:,3).^2.*XYZ(:,:,:,1).*XYZ(:,:,:,2) + B(15)*XYZ(:,:,:,3).^2.*XYZ(:,:,:,1).^2 + ... [4th order] cross terms //z
% %         B(16)*XYZ(:,:,:,1).^3 + B(17)*XYZ(:,:,:,2).^3 + B(18)*XYZ(:,:,:,3).^3 + ... [3rd order] terms    
% %         B(19)*XYZ(:,:,:,1).^2.*XYZ(:,:,:,2) + B(20)*XYZ(:,:,:,1).^2.*XYZ(:,:,:,3) + ... [3rd order] cross terms //x
% %         B(21)*XYZ(:,:,:,2).^2.*XYZ(:,:,:,1) + B(22)*XYZ(:,:,:,2).^2.*XYZ(:,:,:,3) + ... [3rd order] cross terms //y
% %         B(23)*XYZ(:,:,:,3).^2.*XYZ(:,:,:,1) + B(24)*XYZ(:,:,:,3).^2.*XYZ(:,:,:,2) + B(25)*XYZ(:,:,:,1).*XYZ(:,:,:,2).*XYZ(:,:,:,3) + ... [3rd order] cross terms //z
% %         B(26)*XYZ(:,:,:,1).^2 + B(27)*XYZ(:,:,:,2).^2 + B(28)*XYZ(:,:,:,3).^2 + ... [2nd order] terms
% %         B(29)*XYZ(:,:,:,1).*XYZ(:,:,:,2)  + B(30)*XYZ(:,:,:,2).*XYZ(:,:,:,3)  + B(31)*XYZ(:,:,:,1).*XYZ(:,:,:,3) + ... [2nd order] cross terms
% %         B(32)*XYZ(:,:,:,1) + B(33)*XYZ(:,:,:,2) + B(34)*XYZ(:,:,:,3) + B(35);  % [1st and 0th terms] terms
% %     
% %     initial_guess = [repmat(0, [1 34]) 0.5*pn_factor];
% % %     fit_minimum = [repmat(-2.5, [1 34]) -2.5]*pn_factor;
% % %     fit_maximum = [repmat(2.5, [1 34]) 2.5]*pn_factor;
% %     fit_minimum = [repmat(-0.1, [1 15]) repmat(-1.0, [1 10]) repmat(-5, [1 9]) 0]*pn_factor;
% %     fit_maximum = [repmat( 0.1, [1 15]) repmat( 1.0, [1 10]) repmat( 5, [1 9]) 5]*pn_factor;
% %     
% %     %    lsqcurvefit(func handle, initial guess, ROI, original data,                 fit minimum , fit maximum, opts)
% %     B =  lsqcurvefit(surfit_v,    initial_guess,   a, scaled_coil_image(find(mask)), fit_minimum , fit_maximum, opts);
% %        
% %     coil_sens = surfit_s(B,XYZ);
% % 
% %     figure, imagesc([scaled_coil_image(:,:,size(scaled_coil_image,3)/2) coil_sens(:,:,size(scaled_coil_image,3)/2)]) 
% %     hold on; contour([mask(:,:,size(scaled_coil_image,3)/2) mask(:,:,size(scaled_coil_image,3)/2)],'w');
% %     
% %     % == export ==
% %     psn_s.model4_3D          = surfit_s;
% %     psn_s.model4_3D_coefs    = B;clc
% %     psn_s.model4_3D_cords    = XYZ;
       
    
end

end

function [sysString] = run_path_on_sys(sysString)
% % PC test case
%  testa = '\\hl-share.nhlbi.nih.gov\dirhome\RamasawmyR\Scan Data\2016\160607\h5\R2_pSpiral_Flow_venc100_Constant125_s8x1_fa10_shimv.h5'
% % MAC test case
%  testb = '/Volumes/DIRHome/RamasawmyR/Scan Data/2016/160607/h5/R2_pSpiral_Flow_venc100_Constant125_s8x1_fa10_shimv.h5'
% % LINUX test case
% testc = '/home/ramasawmyr/dirhome/Scan Data/2016/160607/h5/R2_pSpiral_Flow_venc100_Constant125_s8x1_fa10_shimv.h5'

if isunix
    if ismac
    % MAC OPERATION
    
    if(isempty(regexp(sysString, 'Volumes', 'once'))) % Are we already in the right format?
        catString_i = regexpi(sysString, 'dirhome');
        temp = sysString(catString_i:end);
        temp(regexp(temp, '\')) = '/';
        sysString = lower(['/Volumes/' temp]); % mac is case-independent
    
    end
    
    else
        % LINUX/UBUNTU OPERATION
        
        % Configured for Ramasawmy's "home" mount:
        % sudo -t cifs //hl-share.nhlbi.nih.gov/DIRHome/RamasawmyR ~/dirhome/ -o username=ramasawmyr,domain=NIH
        % THIS WILL NEED TO BE GENERICALLY CONFIGURED!
        
        if(isempty(regexp(sysString, '/home/ramasawmyr/dirhome/', 'once'))) % Are we already in the right format?
        
        catString_i = regexpi(sysString, 'scan data');
        temp = sysString(catString_i:end);
        temp(regexp(temp, '\')) = '/';
        sysString = lower(['/home/ramasawmyr/dirhome/' temp]); % case-independent?

        end
    end
else
    % PC OPERATION
   
%      if(isempty(regexp(sysString, '\\hl-share.nhlbi.nih.gov', 'once'))) % Are we already in the right format?
%         catString_i = regexpi(sysString, 'dirhome');
%         temp = sysString(catString_i:end);
%         temp(regexp(temp, '/')) = '\';
%         sysString = lower(['\\hl-share.nhlbi.nih.gov\' temp]); % case-independent
%     
%     end
end


end

function dmtx = noise_adjust(nfile, i_rd,  data_samp_time, nhlbi_toolbox)

if isempty(nfile)
    dmtx = diag(ones(1,channels));
else
    noise_test = h5read(nfile,'/dataset/data');
    n_MRD_s = nhlbi_toolbox.h5read_xml(nfile);
    
    disp(['Required Sens Map: ' i_rd.measurementInformation.measurementDependency.measurementID ', Noise ID: ' n_MRD_s.measurementInformation.measurementID]);
    
    nhlbi_toolbox.check_noise_dependency(i_rd,n_MRD_s);
    
    Siemens_rBW = n_MRD_s.acquisitionSystemInformation.relativeReceiverNoiseBandwidth;
    
    n_samples = double(noise_test.head.number_of_samples(1));
    n_channels = double(noise_test.head.active_channels(1));
    
    % assuming Siemens using 2 averages:
%     noise_ind = find(noise_test.head.idx.average==1, 1, 'last');
    noise_ind = 256;
    nt2 = zeros(n_samples, noise_ind, n_channels);
    for i = 1:noise_ind
        nt2(:,i,:)=  double(reshape(complex(noise_test.data{i}(1:2:end), noise_test.data{i}(2:2:end)), [n_samples, 1, n_channels ]));
    end
    
    n_scaling = Siemens_rBW * data_samp_time / (noise_test.head.sample_time_us(1)*1e-6);
    dmtx = ismrm_calculate_noise_decorrelation_mtx(nt2, n_scaling ); % figure,imagesc(abs(dmtx));
    
end

end

function parpool_setup(requested_pp)
if nargin < 1
    % default pool num
    requested_pp = 16;
end

% check current instances
gcp_info = gcp('NoCreate');

% limit num workers to allowed number
pc_info = parcluster('local');
max_pp = pc_info.NumWorkers;
if requested_pp > max_pp; requested_pp = max_pp; end

% set-up parallel pool
if isempty(gcp_info.isvalid)
    
    parpool('local', requested_pp);
    
% else % <optional>
    %     % boost the number of workers if necessary
    %     numw = gcp_info.NumWorkers;
    %
    %     if numw < requested_pp
    %         delete(gcp('nocreate'));
    %         parpool('local', requested_pp);
    %     end
end

end