function ismrmrd_s = read_h5_header(dfile)
% ismrmrd_s = read_h5_header(dfile)
%
% Wrapper to read xml header on ISMRMRD data
% dfile = 'path_file_name.h5'; > converted ISMRMRD format

if ~isstruct(dfile)
    %% Extract header
    xml_string = h5read(dfile, '/dataset/xml');
    ismrmrd_s = xml2hdr(xml_string{1});
    
    % +++ example, can edit output for more common usage: +++
    % % new_output.FOV = [ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.x ...
    % %                   ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.y ...
    % %                   ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.z];
    % % new_output.matrix = [ismrmrd_s.encoding.reconSpace.matrixSize.x ...
    % %                      ismrmrd_s.encoding.reconSpace.matrixSize.y ...
    % %                      ismrmrd_s.encoding.reconSpace.matrixSize.z];
    % % new_output.TE = ismrmrd_s.sequenceParameters.TE;
    % % new_output.TR = ismrmrd_s.sequenceParameters.TR;
    % ++++++
    
else
    %% Print out parameters
    
    disp('########################################');
    disp('    Printing out iRD structure');
    disp('########################################');
    disp(' ');
    
    % printstruct(dfile); % sourced generic structure printing % https://github.com/spunt/printstruct/blob/master/printstruct.m
    
    % Simple MRI-specific print
    disp(['    ' dfile.measurementInformation.protocolName])
    disp(dfile.studyInformation)
    disp(['    Meas ID (iRD) : ' dfile.measurementInformation.measurementID]);
    
    disp(' ');disp('### Sequence parameters ###');
    disp(['TE/TR = ' num2str(dfile.sequenceParameters.TE) '/' num2str(dfile.sequenceParameters.TR) ' ms, Flip angle: ' num2str(dfile.sequenceParameters.flipAngle_deg)])
    
    disp(' ');disp('### Encoding info ###');
    disp(['Encoding: ' dfile.encoding.trajectory]);disp(' ');
    disp(['Encoded Res: ' num2str([dfile.encoding.encodedSpace.matrixSize.x dfile.encoding.encodedSpace.matrixSize.y dfile.encoding.encodedSpace.matrixSize.z])])
    disp(['Encoded FOV: ' num2str([dfile.encoding.encodedSpace.fieldOfView_mm.x  dfile.encoding.encodedSpace.fieldOfView_mm.y dfile.encoding.encodedSpace.fieldOfView_mm.z])])
    disp(['Recon Res: ' num2str([dfile.encoding.reconSpace.matrixSize.x dfile.encoding.reconSpace.matrixSize.y  dfile.encoding.reconSpace.matrixSize.z])])
    disp(['Recon FOV: ' num2str([dfile.encoding.reconSpace.fieldOfView_mm.x dfile.encoding.reconSpace.fieldOfView_mm.y dfile.encoding.reconSpace.fieldOfView_mm.z])])
    disp(['Channels: ' num2str(dfile.acquisitionSystemInformation.receiverChannels)])
    
    disp('Acceleration Factor');
    disp(dfile.encoding.comment.parallelImaging.accelerationFactor)
    
    disp(' ');disp(['System Field ' num2str(dfile.acquisitionSystemInformation.systemFieldStrength_T) '(T) Res Freq: ' num2str(dfile.experimentalConditions.H1resonanceFrequency_Hz/1e6) 'MHz'])
    
    % disp(' ');disp('### User Parameters ###');
    % fntemp = fieldnames(dfile.userParameters);
    % for i = 1:length(fntemp)
    %     disp(fntemp{i})
    %     disp(dfile.userParameters.(matlab.lang.makeValidName(fntemp{i})))
    % end
    
    
    
    
end

end