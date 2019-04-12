function [outputStruct] = h5_read_multigroup(filename, dataSelect)
% Load all data
% [a] = h5_read_multigroup('csm_transfer\out.h5');
% Select data to load
% [b] = h5_read_multigroup('csm_transfer\out.h5', 0);
% Load this data (i.e. > 0)
% [c] = h5_read_multigroup('csm_transfer\out.h5', 9);

%%
if nargin < 1
    [a,b] = uigetfile('*.*');
    filename = [b '\' a]; clear a b;
end
a = h5info(filename);

numRecons = length(a.Groups);
if nargin > 1
    if dataSelect == 0
    disp(sprintf('************\nChoose Data \n************'))
    for nRc = 1:numRecons
        b = a.Groups(nRc).Groups;
        disp(['Image set: ' num2str(nRc) '   Recon ID: ' a.Groups(nRc).Name])
    end
    disp('************');
    dataSelect = input('Choose Image Set: ');
%     dataSelect = str2double(dataSelect);
    end
    loopVector = dataSelect;
else
    loopVector = 1:numRecons;
end

%%

for nRc = loopVector
    b = a.Groups(nRc).Groups;
    disp(['Loading Image set: ' num2str(nRc) '   Recon ID: ' a.Groups(nRc).Name])
    
    m_ind = 0;
    p_ind = 0;
    for i = 1 : length(b)
        temp = b(i).Name;
        ind = regexp(temp, 'image_');
        lnum(i) = length(temp(ind:end));
        
    end
    
    m_s_len = min(lnum);
    p_s_len = max(lnum);
    
    clear Mag_data Pha_data info
    
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
    
    outputStruct.(matlab.lang.makeValidName(['data' num2str(nRc)])).M = Mag_data;
    outputStruct.(matlab.lang.makeValidName(['data' num2str(nRc)])).P = Pha_data;
    outputStruct.(matlab.lang.makeValidName(['data' num2str(nRc)])).i = info;
end

end