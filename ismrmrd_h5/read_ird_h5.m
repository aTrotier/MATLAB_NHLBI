function [outputStruct] = read_ird_h5(filename)
% [output_struct] = read_ird_h5(filename)
% Reads in all images from an offline Gadgetron Recon
%
% R Ramasawmy, NHLBI, Oct 2019

% user UI if empty
if nargin < 1
    [a,b] = uigetfile('*.*');
    filename = [b filesep a]; clear a b;
end
a = h5info(filename);

%%
% check the number of reconstructions
numRecons = length(a.Groups); disp(' '); disp('===================================='); disp(' ');
for iNR = 1:numRecons

    disp(['Loading Image set: ' num2str(iNR) '   Recon ID: ' a.Groups(iNR).Name]); disp(' ');
    outputStruct.(matlab.lang.makeValidName(['group_' num2str(iNR)])).id = a.Groups(iNR).Name;

    % optional ismrmrd header attachment
    % info = h5read(filename, [temp '/header']);
    % outputStruct.(matlab.lang.makeValidName(['group_' num2str(iNR)])).info = info;

    % grab all the data in this recon group
    b = a.Groups(iNR).Groups;
    for i=1:length(b)
        temp = b(i).Name;
        ind = regexp(temp, 'image_'); % just take the image ID
        image_name = temp(ind:end); disp(temp)

         data_temp = squeeze(double(h5read(filename, [temp '/data']))); % simplify for image proc
         outputStruct.(matlab.lang.makeValidName(['group_' num2str(iNR)])).(matlab.lang.makeValidName(image_name)) = data_temp;

    end

    disp(' '); disp('===================================='); disp(' ');

end


end
