function [Mag_data, Pha_data, info] = h5_read_RR(filename)
% Mag_data, Pha_data, info] = h5_read_RR(filename)
% Reads in Physio Interpolated Data, from Gadgetron Recon
% Assumes same number of frames for each heartbeat
if nargin < 1
   [a,b] = uigetfile('*.*');
   filename = [b filesep a]; clear a b;
end
a = h5info(filename);
b = a.Groups.Groups;

m_ind = 0;
p_ind = 0;
for i = 1 : length(b)
    temp = b(i).Name;
    disp(['Group ' num2str(i) ': ' temp]);
    ind = regexp(temp, 'image_');
    lnum(i) = length(temp(ind:end));
    
end

m_s_len = min(lnum);
p_s_len = max(lnum);

for i = 1 : length(b)
    temp = b(i).Name;
    
    if(length(unique(lnum)) > 1) % if phase images are present
        
        if lnum(i) == m_s_len
            % Magnitude Image
            m_ind = m_ind + 1;
            Mag_data(:,:,:,m_ind) = squeeze(double(h5read(filename, [temp '/data'])));
            
            if m_ind == 1
                info = h5read(filename, [temp '/header']);
            end
            
        elseif lnum(i) == p_s_len
            % Phase Image
            
            p_ind = p_ind + 1;
            Pha_data(:,:,:,p_ind) = squeeze(double(h5read(filename, [temp '/data'])));
            
        end
        
    else % Just magnitude
        m_ind = m_ind + 1;
        Mag_data(:,:,:,m_ind) = squeeze(double(h5read(filename, [temp '/data'])));
        
        if m_ind == 1
            info = h5read(filename, [temp '/header']);
        end
        
        Pha_data = [];
    end
end


end