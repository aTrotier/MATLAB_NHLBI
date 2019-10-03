function [weights] = DCF_voronoi_RR(traj, inoutflag, showPlots)
% [weights] = DCF_voronoi_RR(traj,<inoutflag>,<showPlots>)
% size(traj) = [samples, interleaves, 2]
% set inoutflag = 1 for in-out estimation
% assumes a constant velocity gradient at.. some point (linear fit?)

if nargin < 2
    inoutflag=0;
    showPlots=0;
elseif nargin < 3
    showPlots=0;
end

interleaves = size(traj,2);
traj = reshape(traj, size(traj,1)*interleaves,2);

if showPlots % plot voronoi diagram
    figure, 
    voronoi(traj(:,1), traj(:,2)); axis image; zoom(5);
end

[x1, x2] = voronoin(traj);

% patch visualisation % figure
max_kr = max(abs(complex(traj(:,1), traj(:,2))));
test_area = zeros(size(traj,1),1);

for i = 1 : size(traj ,1)
    ind = x2{i}';
    
    vec1 = [x1(ind,1), x1(ind,2)];
    
    % kr-limit correction
    % combine in to complex to calc magnitude and angle 
    comp_vec =  complex(vec1(:,1),vec1(:,2)); 
    
    % remove any infinite values
    inf_i = find(isinf(comp_vec)); 
    vec1(inf_i,:) = []; 
    comp_vec(inf_i) = [];
    
    abs_v = abs(comp_vec);
    angle_v = angle(comp_vec);
    
    % loop through vector of vertices and limit points larger than kr
    for j = 1:length(abs_v)
        if abs_v(j) > max_kr
            [x,y] = pol2cart(angle_v(j),max_kr );
            vec1(j,1) = x;
            vec1(j,2) = y;
        end
    end
    
    % patch visualisation %     patch(vec1(:,1), vec1(:,2), 1);
    
    % calculate area of vertices - equal to density comp
    test_area(i,1) = polyarea( vec1(:,1), vec1(:,2));

end

%%
if inoutflag == 0
    %%
    disp('Performing spiral out and fitted end-DCWs');
    % # guessing the number of arms #
    % This was used initially for uneven sample sizes.. though I can't remember
    % in what case that happened..
    
    %     temp = abs(diff(traj(:,1)));
    %     [a] = find(temp > (mean(temp) + 2*std(temp)) );
    %     % if length(unique(diff(a))) > 1
    %     % uneven sample length, only dealing with two for the mo
    %     interleaves = length(a) + 1;
    %
    %     start_samp = [1 1+a'];
    %     end_samp = [a' size(traj,1)];
    
    samplength = length(test_area)/interleaves;
    start_samp = [1:samplength:length(test_area)];
    end_samp = [samplength:samplength:length(test_area)];
  
    if showPlots
    figure, % ### plot weights and fits ###clear 
    end
    
    weights = [];
    for i = 1:interleaves
        
        ta2 = test_area(start_samp(i):end_samp(i));
        
%         % detect rising edge of outer traj, and fudge it
%         pad = 10;
% %         outer_arm_i = find(diff(ta2) > nanmedian(diff(ta2)),1, 'first');
%         outer_arm_i = find(diff(ta2) > nanmean(diff(ta2)),1, 'first');
%         
%         pad2 = 20; x_samp_i = (outer_arm_i - pad2):(outer_arm_i - pad);
%         ta3 = ta2(x_samp_i);
%         ta3_fit = [ones(size(x_samp_i))' x_samp_i']\ta3;
%         ta3 = ta2;
%         ta3((outer_arm_i - pad):length(ta2)) = ta3_fit(1) + ta3_fit(2)*((outer_arm_i - pad):length(ta2));
        ta3 = ta2;
%         ta3(find(isnan(ta2))) = nanmedian(ta2);
        ta3(find(isnan(ta2))) = median(ta2,'omitnan');
        ta3(find(ta3 > median(ta3)*2 )) = median(ta3);
        
        weights = [weights; ta3];
        
        % ### plot weights and fits ###
        if showPlots
            subplot(ceil(sqrt(interleaves)), ceil(interleaves/(ceil(sqrt(interleaves)))),i);
            plot(ta2, 'k-') % voronoi weights
            hold on, plot(ta3, 'r-'); ylim([0 1.25*max(ta3)]) % fitted output
        end
    end
    
    %% Alternative/legacy methods to fit end
    %     samples =a(1); interleaves = size(traj,1)/samples;
    %
    %     % lazy-boy pad the end
    %     ta2 = test_area(1:samples); pad = 10;
    %     outer_arm_i = find(diff(ta2) > nanmean(diff(ta2)),1, 'first');
    %     outer_arm_weight = ta2(outer_arm_i - pad);
    %     ta2(outer_arm_i-pad : samples) = outer_arm_weight;
    %     figure, plot(test_area(1:samples), 'k-');
    %     hold on, plot(ta2, 'b-'); ylim([0 max(ta2)*1.25])
    %
    %     % linearly fit the end
    %     pad2 = 20; x_samp_i = (outer_arm_i - pad2):(outer_arm_i - pad);
    %     ta3 = test_area(x_samp_i);
    %     ta3_fit = [ones(size(x_samp_i))' x_samp_i']\ta3;
    %     ta3 = test_area(1:samples);
    %     ta3((outer_arm_i - pad):samples) = ta3_fit(1) + ta3_fit(2)*((outer_arm_i - pad):samples);
    %     % figure, plot(test_area(1:samples), 'k-');
    %     hold on, plot(ta3, 'r-'); ylim([0 max(ta3)*1.25]); % xlim([x_samp_i]);
    %
    %     test_area= reshape(test_area, samples, interleaves);
    %     test_area((outer_arm_i - pad):samples,:) = repmat(ta3_fit(1) + ta3_fit(2)*((outer_arm_i - pad):samples)',[1 interleaves]);
    %     weights= reshape(test_area, interleaves*samples, 1);
    %
    % end
    
    
else
    %%
    % Out-in test (with sampling of rewinders?)
    disp('Fudging spiral in-out');
    
    % Legacy code
    %     clear test_area;
    %     for i = 1 : length(x2)
    %         ind = x2{i}';
    %         test_area(i,1) = polyarea( x1(ind,1), x1(ind,2));
    %     end
    %
    %     interleaves = 16;
    %     start_samp  = [1:length(x2)/16:length(x2)];
    %     end_samp    = [length(x2)/16:length(x2)/16:length(x2)];
    
    samplength = length(test_area)/interleaves;
    start_samp = [1:samplength:length(test_area)];
    end_samp = [samplength:samplength:length(test_area)];
    if showPlots
    figure, % ### plot weights and fits ###
    end
    weights = [];
    for i = 1:interleaves
        
        ta2 = test_area(start_samp(i):end_samp(i));
        
        % attempt to smooth out irregularities and edge blow-ups
        temp_a = find(ta2 > median(smooth(ta2) ) + std( ta2(floor(samplength/8):ceil(samplength*(3/8)) ) ) ); % std range 256:800 for 16
        
        ta3 = ta2;
%         ta3(find(isnan(ta2))) = nanmedian(ta2);
        ta3(find(isnan(ta2))) = median(ta2,'omitnan');
        ta3(temp_a) = median(smooth(ta2));
        
        weights = [weights ta3'];
        
        % ### plot weights and fits ###
        if showPlots
        subplot(ceil(sqrt(interleaves)), ceil(interleaves/(ceil(sqrt(interleaves)))),i);
        plot(ta2, 'k-')
        hold on, plot(ta3, 'r-'); ylim([0 1.25*max(ta3)])
        end
    end
    
end

end

