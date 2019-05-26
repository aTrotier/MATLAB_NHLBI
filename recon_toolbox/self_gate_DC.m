function [sg_info] = self_gate_DC(kspace_pw, i_rd, phases)
% Stand-alone test function
% [sg_info] = self_gate_DC(dfile, nfile, phases [#cardiac <30> #resp <1>])
% This does a lot of smoothing for peak detection..
% PCA = 1 for motion 

if nargin < 3
    phases = [30 1];
end
cardiac_frames = phases(1);
resp_phases = phases(2);

% nasty switch for compatibility
if ischar(kspace_pw) && ischar(i_rd) % assume file names are being fed in
    raw_data = h5read(kspace_pw,'/dataset/data');
    dmtx = ismrm_dmtx_RR(i_rd);
    
    i_rd = read_h5_header(kspace_pw);
    TR = 1e-3 * i_rd.sequenceParameters.TR;
    
    
    channels = raw_data.head.available_channels(1);
    samples = raw_data.head.number_of_samples(1);
    
    clear data;
    for i = 1:length(raw_data.data)
        d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
        d = reshape(d, samples, 1, channels);
        data(i,:) = double(squeeze(d(1,1,:)));
    end
    data = ismrm_apply_noise_decorrelation_mtx(data, dmtx);
else
    channels = size(kspace_pw,ndims(kspace_pw));
    
    % take first data point
    S.subs = repmat({':'},1,ndims(kspace_pw)); % handle 2D/3D images
    S.subs{1} = 1;
    S.type = '()';
    data = subsref(kspace_pw, S);
    
end


%% Prep data and PCA (1-channel out)

data = abs(data);
% figure, plot(data)

cmp = jet(single(channels));
channel_legend = cellstr(num2str([1:channels]'));

% average two flow sets together for even phase-contrast binning
if i_rd.encoding.encodingLimits.set.maximum
    TR = (i_rd.encoding.encodingLimits.set.maximum + 1)*TR;
    data = squeeze(mean(reshape(data, [2 length(data)/2 channels]),1));
    % figure, plot(data)
end

% Take SVD
[U,S,V] = svd(data,'econ');
% figure, plot(diag(S))

% % Rank Contribution [ academic : if you want to see which rank of motion data]
% rank_matrix = zeros(size(U));
% for i = 1:length(S)
%     rank_matrix(:,i) = U(:,i)*S(i,i)*V(i,i)';
% end
% figure, colorful_plots(rank_matrix)
%
% figure, clear P_out D_out;
% for i = 1:channels
%     subplot(1,2,1); hold on;
%     plot(smooth(rank_matrix(:,i)), '-', 'Color', cmp(i,:))
%
%     subplot(1,2,2); hold on;
%     Y = fft(smooth( rank_matrix(:,i) ));
%     P2 = abs(Y/length(U));
%     P1 = P2(1:(length(U))/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     f = (1/TR)*(0:((length(U))/2))/(length(U));
%     plot(f, P1, '-', 'Color', cmp(i,:))
% %     P_out(:,i) = P1;
% end
% xlim([0.1 3.6])
% legend(channel_legend)


% PCA
channels_out = 1;
data_pca = U(:,1:channels_out)*S(1:channels_out,1:channels_out)*V(1:channels_out,1:channels_out)';
% figure, plot(data_pca)

smooth_data = 1;
if smooth_data
    for i = 1:channels_out
        data_pca(:,i) = smooth(data_pca(:,i));
    end
end

%% Fourier transform and visualise motion frequencies
figure, clear P_out;
for i = 1:channels_out
    subplot(1,2,1); hold on;
    plot(data_pca(:,i), '-', 'Color', cmp(i,:))
    
    subplot(1,2,2); hold on;
    Y = fft(data_pca);
    P2 = abs(Y/length(U));
    P1 = P2(1:(length(U))/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    freq_range = (1/TR)*(0:((length(U))/2))/(length(U));
    plot(freq_range, P1, '-', 'Color', cmp(i,:))
    P_out(:,i) = P1;
end
xlim([0.1 3.6]); xlabel('Freq (Hz)');
legend(channel_legend)

sg_info.dc_data = data_pca;
sg_info.dc_fft = P_out;

%% CARDIAC binning
tempf = find(0.5 < freq_range & freq_range < 3.5); % MIN/MAX HEART RATE = 30/210 bpm
P_temp = P_out(tempf, :);
[ty, coil_ecg] = find(P_temp == max(P_temp(:)));
freq_ecg = freq_range(ty+min(tempf)-1);

pass_band_ecg = P_out(:,coil_ecg);
px = (1:length(P_out));
px(tempf) = [];
dyn_cardiac = smooth(data_pca(:,coil_ecg));

figure, plot(dyn_cardiac, 'k-'); hold on;
d_inter = round((1/freq_ecg)/TR)+5;
if(~mod(d_inter,2))
    d_inter = d_inter+1;
end
ave_cardiac = smooth(data_pca(:,coil_ecg), d_inter);

plot(ave_cardiac, 'r-')
legend({'Smoothed ECG + RESP','Smoothed RESP'});

% Peaks
[pks, locs] = findpeaks(smooth(data_pca(:,coil_ecg)-smooth(data_pca(:,coil_ecg), d_inter)), 'MinPeakDistance',round(0.75* (1/freq_ecg)/TR) );

% remove steady-state descent
steady_state_inds = find(diff(ave_cardiac)  > median(diff(ave_cardiac)) - std(diff(ave_cardiac)));
steady_state_start = steady_state_inds(1);
steady_state_inds = find(locs > steady_state_start);
pks = pks(steady_state_inds);
locs = locs(steady_state_inds);

    %% Bin each cardiac cycle
rpd_y = pks; rpd_x = locs;

iRRvec = [[rpd_x(1:end-1)'];rpd_x(2:end)'];

bin_data1 = zeros(1,length(data_pca));

%     figure; hold on;
for iRR = 1:length(rpd_x)-1
    RRvec = iRRvec(1, iRR):iRRvec(2, iRR);
    
    temp_vec =  1:length(RRvec);
    bin_data2 = zeros(size(temp_vec)); bin_data2(1) = 1; bin_data2(end) = cardiac_frames;
    bin_data2(2:length(RRvec)-1) =   round( temp_vec(1:end-2)*( 1.001*cardiac_frames/(length(RRvec)-1)) );
    bin_data1(RRvec) = bin_data2';
%         plot(bin_data2);
end


    %% export

sg_info.ecg_binned = bin_data1;
sg_info.ecg_freq = freq_ecg;
sg_info.ecg_chan = coil_ecg;


%% RESPIRATORY binning

tempf = find(min(freq_range(freq_range>0)) < freq_range & freq_range < 0.5); % MIN/MAX RESP RATE = ~0/30 bpm

P_temp = P_out(tempf, :);
[ty, coil_resp] = find(P_temp == max(P_temp(:)));
freq_resp = freq_range(ty+min(tempf)-1);

pass_band_resp = P_out(:,coil_resp);
px = (1:length(P_out));
px(tempf) = [];

d_inter2 = round((1/freq_resp)/TR)+5;
if(~mod(d_inter2,2))
    d_inter2 = d_inter2+1;
end
dyn_resp = smooth(data_pca(:,coil_resp), d_inter);
ave_resp = smooth(smooth(data_pca(:,coil_resp),d_inter2));

figure, plot(dyn_resp, 'k-'); hold on;
plot(ave_resp, 'r-'); plot(dyn_resp - ave_resp + median(ave_resp), 'b--')
legend({'Smoothed RESP','Ave RESP'});

% Peaks detection
[pks, locs] = findpeaks(dyn_resp-ave_resp, 'MinPeakDistance',round(0.75* (1/freq_resp)/TR) );
steady_state_inds = find(locs > steady_state_start);
pks = pks(steady_state_inds);
locs = locs(steady_state_inds);


    %% Respiratory phases

    rpd_y = pks; rpd_x = locs;
    iRRvec = [[rpd_x(1:end-1)'];rpd_x(2:end)'];
    
    bin_data1 = zeros(1,length(data_pca));
    
    % figure, hold on;
    for iRR = 1:length(rpd_x)-1
        
        RRvec = iRRvec(1, iRR):iRRvec(2, iRR);
        temp_vec =  1:length(RRvec);
        temp_resp = dyn_resp - ave_resp;
        resp_binned = floor((temp_resp(RRvec)-min(temp_resp(RRvec)))/(1.001*range(temp_resp(RRvec))/resp_phases));
        bin_data1(RRvec) = resp_binned' + 1;
        
        %     plot(resp_binned);
    end


    %% export
    
    sg_info.resp_binned = bin_data1;
    sg_info.resp_freq = freq_resp;
    sg_info.resp_chan = coil_resp;
    
    figure, hold on, plot(sg_info.resp_binned, 'r-');   plot(sg_info.ecg_binned /(cardiac_frames/resp_phases) , 'b-');

%% Export pad output for 2-sets
if i_rd.encoding.encodingLimits.set.maximum
    temp = sg_info.resp_binned;
    sg_info.resp_binned = reshape(repmat(temp, [2 1]), [1 length(temp)*2]);
    
    temp = sg_info.ecg_binned;
    sg_info.ecg_binned = reshape(repmat(temp, [2 1]), [1 length(temp)*2]);
end



end

function [pks, locs] = self_gate_DC_fp(dyn_data, freq_data, TR)


d_inter = round((1/freq_data)/TR)+5;
if(~mod(d_inter,2))
    d_inter = d_inter+1;
end
dyn_trace = smooth(dyn_data);
ave_trace = smooth(dyn_data,d_inter); % smooth(smooth(x, d_inter))

figure, plot(dyn_trace, 'k-'); hold on;
plot(ave_trace, 'r-'); plot(dyn_trace - ave_trace + median(ave_trace), 'b--')
legend({'Smoothed','Ave'});

% Peaks detection
[pks, locs] = findpeaks(dyn_trace - ave_trace, 'MinPeakDistance',round(0.75* (1/freq_data)/TR) );

end
