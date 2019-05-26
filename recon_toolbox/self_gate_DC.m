function [sg_info] = self_gate_DC(dfile, nfile, phases, TR)
% [sg_info] = self_gate_DC(dfile, nfile, phases [#cardiac <30> #resp <1>])
if nargin < 3
   phases = [30 1]; 
%    TR = 0.0186;
end

raw_data = h5read(dfile,'/dataset/data');
i_rd = read_h5_header(dfile);
TR = 2 * 1e-3 * i_rd.sequenceParameters.TR;

cardiac_frames = phases(1);
resp_phases = phases(2);

dmtx = ismrm_dmtx_RR(nfile);

channels = raw_data.head.available_channels(1);
samples = raw_data.head.number_of_samples(1);

%%
clear data;
for i = 1:length(raw_data.data)
    d = complex(raw_data.data{i}(1:2:end), raw_data.data{i}(2:2:end));
    d = reshape(d, samples, 1, channels);
    data(i,:) = double(squeeze(d(1,1,:)));
end
data = ismrm_apply_noise_decorrelation_mtx(data, dmtx);
% figure, plot(real(data))

cmp = jet(single(channels));
for i = 1:channels
    leg_chan{i} = num2str(i);
end
% 
% figure, hold on,
% for i = 1:channels
%     temp = mean(reshape(real(data(:,i)), [2 length(data)/2]),1);
%     plot(smooth(real(temp)))
% end

figure, clear P_out D_out;
for i = 1:channels
    temp = mean(reshape(abs(data(:,i)), [2 length(data)/2]),1);
    D_out(:,i) = temp;
    subplot(1,2,1); hold on;
    plot(smooth(real(temp)), '-', 'Color', cmp(i,:))

    subplot(1,2,2); hold on;
    Y = fft(smooth(real(temp)));
    P2 = abs(Y/length(temp));
    P1 = P2(1:(length(temp))/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = (1/TR)*(0:((length(temp))/2))/(length(temp));
    plot(f, P1, '-', 'Color', cmp(i,:))
    P_out(:,i) = P1;
end
xlim([0.1 3.6])
legend(leg_chan)

sg_info.dc_data = D_out;
sg_info.dc_fft = P_out;

%% CARDIAC binning
tempf = find(0.5 < f & f < 3.5); % MIN/MAX HEART RATE = 30/210 bpm
P_temp = P_out(tempf, :);
[ty, coil_ecg] = find(P_temp == max(P_temp(:)));
freq_ecg = f(ty+min(tempf)-1);

pass_band_ecg = P_out(:,coil_ecg);
px = (1:length(P_out));
px(tempf) = [];
dyn_cardiac = smooth(D_out(:,coil_ecg));

figure, plot(dyn_cardiac, 'k-'); hold on;
d_inter = round((1/freq_ecg)/TR)+5;
if(~mod(d_inter,2))
    d_inter = d_inter+1;
end
ave_cardiac = smooth(D_out(:,coil_ecg), d_inter);

plot(ave_cardiac, 'r-')
legend({'Smoothed ECG + RESP','Smoothed RESP'});

% Peaks
[pks, locs] = findpeaks(smooth(D_out(:,coil_ecg)-smooth(D_out(:,coil_ecg), d_inter)), 'MinPeakDistance',round(0.75* (1/freq_ecg)/TR) );
% Rising edge
[pks2, locs2] = findpeaks(smooth(diff(smooth(D_out(:,coil_ecg)-smooth(D_out(:,coil_ecg), d_inter)))), 'MinPeakDistance',round(0.75* (1/freq_ecg)/TR) );
[mean(diff(locs)) mean(diff(locs2)); std(diff(locs)) std(diff(locs2)) ]

          %% compare rising edge and peaks
          
%           figure, hold on,
%           x1 = smooth(D_out(:,coil_ecg)-smooth(D_out(:,coil_ecg), d_inter));
%           plot(x1./max(x1(:)), 'b-');
%           plot(locs, pks/max(x1(:)), 'bo');
%           x2 = smooth(diff(smooth(D_out(:,coil_ecg)-smooth(D_out(:,coil_ecg), d_inter))));
%           plot(x2./max(x2(:)), 'r-');
%           plot(locs2, pks2/max(x2(:)), 'ro');
%           plot(smooth(D_out(:,coil_ecg))/max(D_out(:,coil_ecg)), 'k-')
           
    %% Bin each cardiac cycle
    rpd_y = pks; rpd_x = locs; 
% cardiac frames
    
    iRRvec = [[0 rpd_x(1:end-1)']+1;rpd_x'];

    bin_data1 = zeros(1,rpd_x(end));
%     figure;
    for iRR = 1:length(rpd_x)
% %         RRvec = iRRvec(1, iRR):iRRvec(2, iRR);
% %         temp = dyn_cardiac(RRvec);
% %         bin_step = range(temp)/(cardiac_frames);
% %
% %         bin_data2 = 1+floor((temp-min(temp))/(bin_step*1.001));
% %
% %         bin_data1(RRvec) = bin_data2';

% simpler method, ?? convert to "rpd"??
         RRvec = iRRvec(1, iRR):iRRvec(2, iRR);
%       RRvec2 = floor((1:(length(21:32)*100))/(1.001*(32 - 21)/(cardiac_frames-1)))
%          bin_data2 = floor( (1:length(RRvec)) / (1.001*(iRRvec(2, iRR) - iRRvec(1, iRR))/(cardiac_frames-1)) );
         temp_vec =  1:length(RRvec);
         bin_data2 = zeros(size(temp_vec)); bin_data2(1) = 1; bin_data2(end) = cardiac_frames;
         bin_data2(2:length(RRvec)-1) =   round( temp_vec(1:end-2)*( 1.001*cardiac_frames/(length(RRvec)-1)) );
      bin_data1(RRvec) = bin_data2';

%                 plot(bin_data2);
    end


    %% export

    sg_info.ecg_binned = bin_data1;
    sg_info.ecg_freq = freq_ecg;
    sg_info.ecg_chan = coil_ecg;


  

%% RESP binning
tempf = find(min(f(f>0)) < f & f < 0.5); % MIN/MAX RESP RATE = ~0/30 bpm

P_temp = P_out(tempf, :);
[ty, coil_resp] = find(P_temp == max(P_temp(:)));
freq_resp = f(ty+min(tempf)-1);

pass_band_resp = P_out(:,coil_resp);
px = (1:length(P_out));
px(tempf) = [];

d_inter2 = round((1/freq_resp)/TR)+5;
if(~mod(d_inter2,2))
    d_inter2 = d_inter2+1;
end
dyn_resp = smooth(D_out(:,coil_resp), d_inter);
ave_resp = smooth(smooth(D_out(:,coil_resp),d_inter2));

figure, plot(dyn_resp, 'k-'); hold on;
plot(ave_resp, 'r-'); plot(dyn_resp - ave_resp + median(ave_resp), 'b--')
legend({'Smoothed RESP','Ave RESP'});

    % Peaks detection
    [pks, locs] = findpeaks(dyn_resp-ave_resp, 'MinPeakDistance',round(0.75* (1/freq_resp)/TR) );
    % Rising edge
    [pks2, locs2] = findpeaks(smooth(diff( smooth( dyn_resp-ave_resp ))), 'MinPeakDistance',round(0.75* (1/freq_resp)/TR) );
%     [[mean(diff(locs)) mean(diff(locs2))]; [std(diff(locs)) std(diff(locs2))]]

  %%
% 
%     figure, hold on,
%     x1 = dyn_resp-ave_resp;
%     plot(x1./max(x1(:)), 'b-');
%     plot(locs, x1(locs), 'bo');
%     x2 = smooth(diff( smooth( dyn_resp-ave_resp )));
%     plot(x2./max(x2(:)), 'r-');
%     plot(locs2, x2(locs), 'ro');
%     plot(smooth(D_out(:,coil_resp) - min(D_out(:,coil_resp)))/(max(D_out(:,coil_resp)) - min(D_out(:,coil_resp))), 'k-')
% 

    %% Respiratory phases
    
    % remove non-steady state signal (0.5 = aggressive)
%     ss_start = find(ave_resp < 1/2*std(ave_resp) + median(ave_resp), 1, 'first');
%     dyn_resp = dyn_resp(ss_start:end);
%     ave_resp = ave_resp(ss_start:end);
%     figure, plot(dyn_resp)

 rpd_y = pks; rpd_x = locs; 
 iRRvec = [[0 rpd_x(1:end-1)']+1;rpd_x'];

    bin_data1 = zeros(1,rpd_x(end));
%     figure;
    for iRR = 1:length(rpd_x)

         RRvec = iRRvec(1, iRR):iRRvec(2, iRR);
%       RRvec2 = floor((1:(length(21:32)*100))/(1.001*(32 - 21)/(cardiac_frames-1)))
%          bin_data2 = floor( (1:length(RRvec)) / (1.001*(iRRvec(2, iRR) - iRRvec(1, iRR))/(cardiac_frames-1)) );
         temp_vec =  1:length(RRvec);
         temp_resp = dyn_resp - ave_resp;
         resp_binned = floor((temp_resp(RRvec)-min(temp_resp(RRvec)))/(1.001*range(temp_resp(RRvec))/resp_phases));
         bin_data1(RRvec) = resp_binned' + 1;

%                 plot(bin_data2);
    end

%     resp_binned = floor((dyn_resp-min(dyn_resp))/(range(dyn_resp)/resp_phases));
%     resp_binned(resp_binned==resp_phases) = resp_phases-1;


    %% export
sg_info.resp_binned = bin_data1;
sg_info.resp_freq = freq_resp;
sg_info.resp_chan = coil_resp;

figure, hold on, plot(sg_info.resp_binned, 'r-');   plot(sg_info.ecg_binned /(cardiac_frames/resp_phases) , 'b-');   


%% ?/??/???
% temp_data = data(201:600,:);
% figure, plot(abs(temp_data))
% [u, s, v] = svd(temp_data, 'econ');
% figure, plot(abs(u))

 
%% PCA - Liheng Guo <need to optimise implementation>

% Guo_flag  = 0;
% if Guo_flag
%     figure,
%     for i = 1:5
%         eigVecs = LG_eigNavTrace(D_out, i);
%         subplot(5,2,2*i -1), plot(smooth(eigVecs(:,1)));
%         subplot(5,2,2*i), plot(smooth(eigVecs));
%     end
%     
%     eigVecs = LG_eigNavTrace(D_out, 8);
% [LG_cardiac, cardiac_cycles] = LG_CardCurve_eig(eigVecs, TR*1000, [50 210]);
% % [LG_cardiac, cardiac_cycles] = LG_CardCurve_eig(eigVecs(200:end,:), TR*1000, [100 210]);
% [LG_resp, resp_cycles] = LG_CardCurve_eig(eigVecs, TR*1000, [min(f(f>0))*60 30]);
% cardiac_freq  = cardiac_cycles/(TR*length(LG_cardiac));
% resp_freq  = resp_cycles/(TR*length(LG_cardiac));
% figure, plot(smooth(LG_cardiac), 'r-'); hold on, plot(smooth(LG_resp), 'b-')
% figure, plot(smooth(D_out(:,coil_ecg)), 'r-'), hold on, plot(smooth(D_out(:,coil_resp)), 'b-')
% 
%     
% end
% 
% 
% [p1, l1] = self_gate_DC_fp(LG_cardiac, cardiac_freq, TR);
% [p2, l2] = self_gate_DC_fp( D_out(:,coil_ecg), freq_ecg, TR); close all
% 
% [pr1, lr1] = self_gate_DC_fp( smooth(LG_resp, d_inter), freq_resp, TR);
% [pr2, lr2] = self_gate_DC_fp( smooth(D_out(:,coil_resp), d_inter), freq_resp, TR);

%%
% %
% % figure, clear temp;
% % for i = [1 11]
% %     temp = mean(reshape(real(data(201:end,i)), [2 length(data(201:end,1))/2]),1);
% %     subplot(1,2,1); hold on;
% %     plot(smooth(real(temp)))
% %
% %     subplot(1,2,2); hold on;
% %     Y = fft(smooth(real(temp)));
% %     P2 = abs(Y/length(temp));
% %     P1 = P2(1:(length(temp))/2+1);
% %     P1(2:end-1) = 2*P1(2:end-1);
% %     f = (1/TR)*(0:((length(temp))/2))/(length(temp));
% %     plot(f, P1)
% % end
% %
% % find(f > 1.72,1, 'first')
% %
% % P2 = P1;
% % P2(29:end) =0;
% % figure, plot(abs(ifft(P2)))
% %
% % [frequencies, dampings, basis, ahat] = hsvd_new(temp,55);
% %
% % dt = 1 /55;
% % t = (0:dt:(length(temp)-1)*dt);
% % tempf = find(0.5 < frequencies & frequencies < 3.5);
% %
% % figure,plot(abs(basis(:,tempf)*ahat(tempf)))
% % a = findpeaks(abs(basis(:,tempf)*ahat(tempf)))
% % figure
% % [b,c] = findpeaks(smooth(temp));
%
%% Ipshitta

% figure, clear x_vec;
% for i = 1:channels
%     temp = mean(reshape(abs(data(:,i)), [2 length(data)/2]),1);
%     temp = temp./max(temp(:));
%     x_vec(:,i) = temp;
% end
% plot(x_vec);
% % Take SVD
% [U,S,V] = svd(x_vec,'econ');
% figure, plot(diag(S))
% 
% % Rank Contribution
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
%     P2 = abs(Y/length(temp));
%     P1 = P2(1:(length(temp))/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     f = (1/TR)*(0:((length(temp))/2))/(length(temp));
%     plot(f, P1, '-', 'Color', cmp(i,:))
% %     P_out(:,i) = P1;
% end
% xlim([0.1 3.6])
% legend(leg_chan)
% 
% 
% % Raj's PCA
% pca_chan = 4;
% figure, plot(U(:,1:pca_chan)*S(1:pca_chan,1:pca_chan)*V(1:pca_chan,1:pca_chan)')
% 
% % tructate SVD
% 
% rank_K = 5;
% 
% s_diag = diag(S);
% s_diag_trunc = zeros(size(s_diag));
% s_diag_trunc(1:rank_K) = s_diag(1:rank_K);
% 
% S_trunc = zeros(size(S));
% 
% S_trunc(1:size(x_vec,2),1:size(x_vec,2)) = diag(s_diag_trunc);
% 
% x_rec_vec = U*S_trunc*V';
% figure, colorful_plots(x_rec_vec)
% 
% 
% 
% % X_rec = reshape(x_rec_vec,[size(x,1),size(x,2),size(x,3)]);
% X_rec = x_rec_vec;
% %  plot recon
% figure;
% 
% slc = 1;
% 
% imagesc([x(:,:,slc),X_rec(:,:,slc)]);colormap(gray);
% 
% figure,plot(X_rec);
% 
% sp1 = floor(sqrt(double(channels)));sp2 = ceil(double(channels)/sp1);
% 
% figure,
% for i = 1:double(channels)
%     subplot(sp1,sp2,i); hold on;
%     plot(x_vec(:,i), 'k-'),
%     plot(x_rec_vec(:,i), 'r-')
% end
% figure, plot(U)

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
