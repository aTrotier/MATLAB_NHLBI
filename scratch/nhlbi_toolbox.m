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


function dmtx = noise_adjust(nfile, i_rd,  data_samp_time, nhlbi_toolbox)

if isempty(nfile)
    dmtx = diag(ones(1,channels));
else
    noise_test = h5read(nfile,'/dataset/data');
    iRD_s = nhlbi_toolbox.h5read_xml(nfile);
    
    disp(['Required Sens Map: ' i_rd.measurementInformation.measurementDependency.measurementID ', Noise ID: ' iRD_s.measurementInformation.measurementID]);
    
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