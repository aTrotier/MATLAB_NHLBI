function [bin_data, NominalInterval]= physio_Binning(rpd, frames)
% [bin_data]= physio_Binning(raw_data, frames)
% raw_data = raw_data.head.physiology_time_stamp;
% Rudimentary missed beat detection

%%
if nargin < 2
    frames = 30;
end

% rpd = raw_data.head.physiology_time_stamp;
rpd = 2.5*double(rpd(1,:)); % Siemens unit of 2.5 ms.

%% find peaks for 'beat-to-beat' binning
[rpd_y, rpd_x] = findpeaks(rpd);% [y, x] = findpeaks(data, 'MinPeakDistance', minPeakWidth,'MinPeakHeight', minPeakHeight);
%   figure, plot(rpd, 'b-'); hold on, plot(rpd_x, rpd_y, 'ro');
if length(rpd_y) < 2
    bin_data = 0;
    disp('No physio trace!'); figure('Name', 'Physio trace'), plot(rpd, 'b-');
    NominalInterval = 0;
else
    iRRvec = [[0 rpd_x(1:end-1)]+1;rpd_x]';
    
    %% "Arrythmia" detection ### M + 1*STD to pick out missed beats
    
    mean_RR = mean(rpd_y);
    std_RR = std(rpd_y);
%     aRR = find(rpd_y > mean_RR + std_RR); % ### MAY NEED TUNING ###
    aRR = find(rpd_y > median(rpd_y)*1.5); % ### MAY NEED TUNING ###
    
    % % RR alternative:
    % aRR = find(abs(rpd_y - mean_RR) > 1.25*std_RR); % ### MAY NEED TUNING ###
    
    % % Hui method:
    % med_RR = median(rpd_y);
    % aRR = find(abs(rpd_y - med_RR) > 0.5*med_RR); % ### MAY NEED TUNING ###
    
    rpd_xt = rpd_x;  %iRRvect = iRRvec; % rpd_yt = rpd_y;
    for i = 1:length(aRR)
        iRRvec(rpd_xt==rpd_xt(aRR(i)),:) = 0;
        rpd_x(rpd_xt==rpd_xt(aRR(i))) = 0;
        rpd_y(rpd_xt==rpd_xt(aRR(i))) = 0;
    end
    
    iRRvec = reshape(iRRvec(iRRvec~=0),length(iRRvec(iRRvec~=0))/2,2);
    rpd_x = rpd_x(rpd_x~=0);
    rpd_y = rpd_y(rpd_y~=0);
    NominalInterval = mean(rpd_y);
    
    %% interpret to local max (fudge factor weighted to final bin)
    
    bin_data1 = zeros(1,rpd_x(end)); % cleanervec = [];
    for iRR = 1:length(iRRvec)
        max_rpd = (rpd_y(iRR));
        bin_step = max_rpd/(frames);
        RRvec = iRRvec(iRR,1):iRRvec(iRR,2);
        bin_data2 = 1+floor((rpd(RRvec)/(bin_step*1.001)));
        %     bin_data2 = round((rpd(RRvec)/bin_step)-0.49);
        bin_data1(RRvec) = bin_data2;
        %     cleanervec = [cleanervec RRvec];
    end
    
    % % clean up
    % % % RRvec =[];
    % % % for i = 1:length(aRR)
    % % %     RRvsec = [RRvec iRRvect(aRR(i),1):iRRvect(aRR(i),2)];
    % % % end
    
    % bin_data1 = bin_data1(cleanervec);
    
    % figure, subplot(1,2,1), plot(bin_data1),ylim([0 frames]);subplot(1,2,2),hist(bin_data1,0:1:frames); xlim([1 frames])
    
    % bin to cell array
    bin_data=cell(1,frames);
    for i= 1:frames
        bin_data{i} = find(bin_data1 == i);
    end
end

end

% % % % arrthymia model
% % % rpd_orig = rpd;
% % % rpd(594:661) = rpd(594:661) + rpd(593);
% % % rpd(1409:1483) = rpd(1409:1483) + rpd(1408);

% % % % interpret to max
% % % max_rpd = max(rpd);
% % % bin_step = max_rpd/frames;
% % % % bins = bin_step:bin_step:max_rpd;
% % % bin_data1 = round(rpd/bin_step); % figure, plot(round(rpd/bin_step))
% % %
% % % for i= 1:frames
% % % bin_data{i} = find(bin_data1 == i);
% % % end
