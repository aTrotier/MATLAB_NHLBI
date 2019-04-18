function [dmtx, n_stats] = ismrm_dmtx_RR(nfile, data_samp_time)
% function dmtx = ismrm_dmtx_RR(nfile, data_samp_time)
% Process noise data and return decorrelation matrix
% R Ramasawmy NHLBI - I don't know if this is even right. 

%% Estimate 
noise_test = h5read(nfile,'/dataset/data');
iRD_s = read_h5_header(nfile);
disp('Noise ID: ');
disp(iRD_s.measurementInformation.measurementID);

Siemens_rBW = iRD_s.acquisitionSystemInformation.relativeReceiverNoiseBandwidth;
% Siemens_rBW = 0.79;

n_samples = double(noise_test.head.number_of_samples(1));
n_channels = double(noise_test.head.active_channels(1)); 
% assuming Siemens using 2 averages:
noise_ind = find(noise_test.head.idx.average==1, 1, 'last');

% % pe_tab = 1+unique(noise_test.head.idx.kspace_encode_step_1(1:noise_ind));
% % nt2 = zeros(n_samples, max(pe_tab), n_channels, 2); 
% % for i = 1:noise_ind
% %     av = noise_test.head.idx.average(i) + 1;
% %     nt2(:,noise_test.head.idx.kspace_encode_step_1(i)+1,:, av)= double(reshape(complex(noise_test.data{i}(1:2:end), noise_test.data{i}(2:2:end)), n_samples,1, n_channels ));
% % end
% % nt2 = mean(nt2, 4); 

%% Plot noise STD
nt2 = zeros(n_samples, noise_ind, n_channels); 
for i = 1:noise_ind
     nt2(:,i,:)=  double(reshape(complex(noise_test.data{i}(1:2:end), noise_test.data{i}(2:2:end)), [n_samples, 1, n_channels ]));
end

if nargout == 2
    nt3 = reshape(nt2,size(nt2,1)*size(nt2,2), size(nt2,3));
    n_stats = std(nt3,[],1);
    
    figure, hold on, title('Coil Noise SD (m +/- std)'); xlabel('Coils'); ylabel('Noise SD');
    plot(1:n_channels, n_stats, 'bo');
    plot([ones(n_channels,1).*mean(n_stats)], 'r-');
    plot([ones(n_channels,1).*(mean(n_stats)-std(n_stats)) ones(n_channels,1).*(mean(n_stats)+std(n_stats))], 'r--');
    
    str = {['Mean ' num2str(mean(n_stats))] ['Std  ' num2str(std(n_stats))]};dim = [.2 .5 .3 .3];
    annotation('textbox',dim,'String',str,'FitBoxToText','on')
end

%%

if nargin < 2
    dmtx = ismrm_calculate_noise_decorrelation_mtx(nt2); % figure,imagesc(abs(dmtx));
else
    % Required relative sampling time for accurate pre-whitening
%     n_scaling = 0.79 *data_samp_time/ (noise_test.head.sample_time_us(1)*1e-6);
    n_scaling = Siemens_rBW * data_samp_time / (noise_test.head.sample_time_us(1)*1e-6);
    dmtx = ismrm_calculate_noise_decorrelation_mtx(nt2, n_scaling ); % figure,imagesc(abs(dmtx));
end

%% 3D coil scaling 
% Don't know what to do with this... 
% encoded/recon dimensions available in xml

% noise_ind = noise_ind+1:length(noise_test.data);
% 
% n_samples = double(noise_test.head.number_of_samples(noise_ind(1)));
% n_channels = double(noise_test.head.active_channels(noise_ind)); 
% pe_tab_1 = 1+(noise_test.head.idx.kspace_encode_step_1(noise_ind));
% pe_tab_2 = 1+(noise_test.head.idx.kspace_encode_step_2(noise_ind));
% 
% nt3_AC = zeros(n_samples, max(pe_tab_1), max(pe_tab_2), n_channels(1)); 
% nt3_BC = zeros(n_samples, max(pe_tab_1), max(pe_tab_2), n_channels(2)); 
% 
% for i = noise_ind
%     if noise_test.head.active_channels(i) == 2
%         nt3_BC(:,noise_test.head.idx.kspace_encode_step_1(i)+1,noise_test.head.idx.kspace_encode_step_2(i)+1,:)= double(reshape(complex(noise_test.data{i}(1:2:end), noise_test.data{i}(2:2:end)), [n_samples 1 1 2]));
%     else
%         nt3_AC(:,noise_test.head.idx.kspace_encode_step_1(i)+1,noise_test.head.idx.kspace_encode_step_2(i)+1,:)= double(reshape(complex(noise_test.data{i}(1:2:end), noise_test.data{i}(2:2:end)), [n_samples 1 1 noise_test.head.active_channels(i)]));
%     end
% end
% 
% img_AC = ismrm_transform_kspace_to_image(nt3_AC);
% img_BC = ismrm_transform_kspace_to_image(nt3_BC);
% img_AC = sqrt(sum(img_AC.*conj(img_AC),4));
% img_BC = sqrt(sum(img_BC.*conj(img_BC),4));
% 
% % montage_RR(img_AC);montage_RR(img_BC);
% scaled_coil_image = img_AC./img_BC;
% mask = img_AC > std(img_AC(:))+mean(img_AC(:)); montage_RR(mask)
% montage_RR(scaled_coil_image.*mask);

end


