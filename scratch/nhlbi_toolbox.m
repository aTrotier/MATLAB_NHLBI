% nhlbi_toolbox.m
% called by nhlbi_toolbox script
% This is real nasty.

% -- Search for local functions --  
function fh = nhlbi_toolbox
    fh = localfunctions;
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
% fig_handle = figure;
figure('Name', 'Experiment header');
subplot(3,3,1); plot(1+double(raw_data.head.idx.kspace_encode_step_1)); title('kspace step 1')
subplot(3,3,2); plot(1+double(raw_data.head.idx.kspace_encode_step_2)); title('kspace step 2')
subplot(3,3,3); plot(1+double(raw_data.head.idx.average)); title('average')
subplot(3,3,4); plot(1+double(raw_data.head.idx.slice)); title('slice')
subplot(3,3,5); plot(1+double(raw_data.head.idx.set)); title('contrast')
subplot(3,3,6); plot(1+double(raw_data.head.idx.phase)); title('phase')
subplot(3,3,7); plot(1+double(raw_data.head.idx.repetition)); title('repetition')
subplot(3,3,8); plot(1+double(raw_data.head.idx.set)); title('set')
subplot(3,3,9); plot(1+double(raw_data.head.idx.set)); title('segment')
% subplot(3,2,3); plot(1+double(raw_data.head.idx.set)); title('user')

end

function [macString] = RR_run_on_mac(pcString)
% PC test case
% testa = '\\hl-share.nhlbi.nih.gov\dirhome\RamasawmyR\Scan Data\2016\160607\h5\R2_pSpiral_Flow_venc100_Constant125_s8x1_fa10_shimv.h5'
% MAC test case 
% testb = '/Volumes/DIRHome/RamasawmyR/Scan Data/2016/160607/h5/R2_pSpiral_Flow_venc100_Constant125_s8x1_fa10_shimv.h5'

if isunix
    if(isempty(regexp(pcString, 'Volumes', 'once'))) % I cant remember what this step is for..
        catString_i = regexpi(pcString, 'dirhome');
        temp = pcString(catString_i:end);
        temp(regexp(temp, '\')) = '/';
        macString = lower(['/Volumes/' temp]); % mac is case-independent
    else
        % assume string is already in mac format
    end
else 
    macString = pcString;    
end

end

function dmtx = noise_adjust(nfile, i_rd,  data_samp_time, nhlbi_toolbox)

if isempty(nfile)
    dmtx = diag(ones(1,channels));
else
    noise_test = h5read(nfile,'/dataset/data');
    iRD_s = nhlbi_toolbox.h5read_xml(nfile);
    
    disp(['Required Sens Map: ' i_rd.measurementInformation.measurementDependency.measurementID ', Noise ID: ' iRD_s.measurementInformation.measurementID]);
    
    nhlbi_toolbox.check_noise_dependency(i_rd,iRD_s);
    
    Siemens_rBW = iRD_s.acquisitionSystemInformation.relativeReceiverNoiseBandwidth;
    
    n_samples = double(noise_test.head.number_of_samples(1));
    n_channels = double(noise_test.head.active_channels(1));
    
    % assuming Siemens using 2 averages:
    noise_ind = find(noise_test.head.idx.average==1, 1, 'last');
    
    nt2 = zeros(n_samples, noise_ind, n_channels);
    for i = 1:noise_ind
        nt2(:,i,:)=  double(reshape(complex(noise_test.data{i}(1:2:end), noise_test.data{i}(2:2:end)), [n_samples, 1, n_channels ]));
    end
    
    n_scaling = Siemens_rBW * data_samp_time / (noise_test.head.sample_time_us(1)*1e-6);
    dmtx = ismrm_calculate_noise_decorrelation_mtx(nt2, n_scaling ); % figure,imagesc(abs(dmtx));
    
end

end