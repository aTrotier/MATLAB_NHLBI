function [R_DCS_PCS] = lookup_PCS_DCS(PCS_description)
% function [R_DCS_PCS] = lookup_PCS_DCS(PCS_description)
% Mapping Patient Coordinate System (PCS) to Device Coordinate System (DCS)/Physical Coordinate System
% 
% String definitions from ismrmrd
% https://github.com/ismrmrd/ismrmrd/ > schema/ismrmrd.xsd
% 
% R_DCS_PCS definitions from Siemens' Idea User Guide
% 
% PCS_description = "STRING" describing patient coordinate system
%   i.e. PCS_description = 'HFS' > head-first-supine
%   i.e. ismrmrd_header.measurementInformation.patientPosition
% R_DCS_PCS = [3x3] coordinate system
%  i.e. R_DCS_PCS = [0 1 0; 1 0 0; 0 0 1] > (P;R;S)
%
% This needs hands & pictures to understand..
% To get    R_DCS_GCS (Device, Gradient)
% you need  R_DCS_PCS (Device, Patient)
% and       R_PCS_GCS (Patient, Gradient) i.e. logical PRS (phase-read-slice) 
% so the physical axis (i.e. to GIRF)
% R_DCS_GCS = R_DCS_PCS*R_PCS_GCS
% "physical" = "patient" x "logical"
%
% There is an 80% chance this is incorrect. 
% 
% Ramasawmy NHLBI 05-2020

switch PCS_description
    case 'HFP'
        % Head first /  prone
        R_DCS_PCS = [   0  1  0
                       -1  0  0
                        0  0 -1  ];
    case 'HFS'
        
        % Head first /  supine
        R_DCS_PCS = [   0 -1  0
                        1  0  0
                        0  0 -1  ];
    case 'HFDR'
        % Head first /  right lateral
        R_DCS_PCS = [   1  0  0
                        0  1  0
                        0  0 -1];
    case 'HFDL'
        % Head first /  left lateral
        R_DCS_PCS = [  -1  0  0
                        0 -1  0
                        0  0 -1];
    case 'FFP'
        % Feet first /  prone
        R_DCS_PCS = [   0  1  0
                        1  0  0
                        0  0  1  ];
    case 'FFS'
        % Feet first /  supine 
        R_DCS_PCS = [   0 -1  0
                       -1  0  0
                        0  0  1  ];
    case 'FFDR'
        % Feet first /  right lateral
        R_DCS_PCS = [   1  0  0
                        0 -1  0
                        0  0  1  ];
    case 'FFDL'
        % Feet first /  left lateral
        R_DCS_PCS = [  -1  0  0
                        0  1  0
                        0  0  1  ];

end


end