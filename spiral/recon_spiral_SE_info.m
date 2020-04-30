function [seheader] = recon_spiral_SE_info(dfile)
% [mag_RT, seheader, complex_out] = recon_spiral_SE_beta(dfile,  nfile, SpiDes)
% Spin-echo sequence reconstruction (uses user_float field to determine TSE
% factor, FOV, matrix etc).
% R Ramasawmy NHLBI Aug 2018

%% Set up
raw_data = h5read(dfile,'/dataset/data');
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
interleaves = double(max(raw_data.head.idx.kspace_encode_step_1))+1; 
reps = (1 + double(max(raw_data.head.idx.repetition)));
averages = (1 + double(max(raw_data.head.idx.average)));
slices = (1 + double(max(raw_data.head.idx.slice)));
% interleaves = gangles/averages;

dt = raw_data.head.sample_time_us(1)*1e-6;

disp(['Samples: ' num2str(samples)])
disp(['Interleaves: ' num2str(interleaves)])
disp(['Channels: ' num2str(channels)])
disp(['BW: ' num2str(dt)])
disp(['Reps: ' num2str(reps)])
disp(['Averages: ' num2str(averages)])
disp(['Slices: ' num2str(slices)])

seheader.samples = samples;
seheader.dt = dt;
seheader.number_aqs = length(raw_data.data);
seheader.averages = averages;
seheader.channels = channels;

end