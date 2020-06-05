function h5_hdr_structure
% example header structure (spiral, reduced)
% <strong>hdr_structure </strong>=  ...
%    |    
%    |--- <strong>acquisitionSystemInformation</strong>
%    |       |    
%    |       |--- coilLabel(1)
%    |       |       |    
%    |       |       |---- coilName : 'HeadNeck_16:1:HL_4'
%    |       |       |-- coilNumber : 2
%    |       |    
%    |       |--- coilLabel(2)
%    |       |       |    
%    |       |       |---- coilName : 'HeadNeck_16:1:HL_1'
%    |       |       |-- coilNumber : 1
%    |       |    
%    |       |--- coilLabel(3)
%    |       |       |    
%    |       |       |---- coilName : 'HeadNeck_16:1:HL_2'
%    |       |       |-- coilNumber : 3
%    |       |    
%    |       |--- coilLabel(4)
%    |       |       |    
%    |       |       |---- coilName : 'HeadNeck_16:1:HL_5'
%    |       |       |-- coilNumber : 4
%    |       |    
%    |       |--- coilLabel(5)
%    |       |       |    
%    |       |       |---- coilName : 'HeadNeck_16:1:HL_3'
%    |       |       |-- coilNumber : 6
%    |       |    
%    |       |--- coilLabel(6)
%    |       |       |    
%    |       |       |---- coilName : 'HeadNeck_16:1:HL_6'
%    |       |       |-- coilNumber : 5
%    |       |    
%    |       |----------------- institutionName : 'National'
%    |       |--------------------- systemModel : 'Investigational_Device_VE11S_S118'
%    |       |-------------------- systemVendor : 'SIEMENS'
%    |       |---------------- receiverChannels : 6
%    |       |-- relativeReceiverNoiseBandwidth : 0.793
%    |       |----------- systemFieldStrength_T : 0.55
%    |    
%    |--- <strong>encoding</strong>
%    |       |    
%    |       |--- encodedSpace
%    |       |       |    
%    |       |       |--- fieldOfView_mm
%    |       |       |       |    
%    |       |       |       |-- x : 208
%    |       |       |       |-- y : 208
%    |       |       |       |-- z : 8
%    |       |       |    
%    |       |       |--- matrixSize
%    |       |       |       |    
%    |       |       |       |-- x : 160
%    |       |       |       |-- y : 160
%    |       |       |       |-- z : 1
%    |       |       |    
%    |       |    
%    |       |--- encodingLimits
%    |       |       |    
%    |       |       |--- average
%    |       |       |       |    
%    |       |       |       |--- center : 0
%    |       |       |       |-- maximum : 0
%    |       |       |       |-- minimum : 0
%    |       |       |    
%    |       |       |--- contrast
%    |       |       |       |    
%    |       |       |       |--- center : 0
%    |       |       |       |-- maximum : 0
%    |       |       |       |-- minimum : 0
%    |       |       |    
%    |       |       |--- kspace_encoding_step_1
%    |       |       |       |    
%    |       |       |       |--- center : 64
%    |       |       |       |-- maximum : 7
%    |       |       |       |-- minimum : 0
%    |       |       |    
%    |       |       |--- kspace_encoding_step_2
%    |       |       |       |    
%    |       |       |       |--- center : 0
%    |       |       |       |-- maximum : 0
%    |       |       |       |-- minimum : 0
%    |       |       |    
%    |       |       |--- phase
%    |       |       |       |    
%    |       |       |       |--- center : 0
%    |       |       |       |-- maximum : 0
%    |       |       |       |-- minimum : 0
%    |       |       |    
%    |       |       |--- repetition
%    |       |       |       |    
%    |       |       |       |--- center : 0
%    |       |       |       |-- maximum : 19
%    |       |       |       |-- minimum : 0
%    |       |       |    
%    |       |       |--- segment
%    |       |       |       |    
%    |       |       |       |--- center : 0
%    |       |       |       |-- maximum : 0
%    |       |       |       |-- minimum : 0
%    |       |       |    
%    |       |       |--- set
%    |       |       |       |    
%    |       |       |       |--- center : 0
%    |       |       |       |-- maximum : 0
%    |       |       |       |-- minimum : 0
%    |       |       |    
%    |       |       |--- slice
%    |       |       |       |    
%    |       |       |       |--- center : 0
%    |       |       |       |-- maximum : 0
%    |       |       |       |-- minimum : 0
%    |       |       |    
%    |       |    
%    |       |--- parallelImaging
%    |       |       |    
%    |       |       |--- accelerationFactor
%    |       |       |       |    
%    |       |       |       |-- kspace_encoding_step_1 : 1
%    |       |       |       |-- kspace_encoding_step_2 : 1
%    |       |       |    
%    |       |       |----- calibrationMode : 'other'
%    |       |    
%    |       |--- reconSpace
%    |       |       |    
%    |       |       |--- fieldOfView_mm
%    |       |       |       |    
%    |       |       |       |-- x : 208
%    |       |       |       |-- y : 208
%    |       |       |       |-- z : 8
%    |       |       |    
%    |       |       |--- matrixSize
%    |       |       |       |    
%    |       |       |       |-- x : 160
%    |       |       |       |-- y : 160
%    |       |       |       |-- z : 1
%    |       |       |    
%    |       |    
%    |       |--- trajectoryDescription
%    |       |       |    
%    |       |       |--- userParameterDouble(1)
%    |       |       |       |    
%    |       |       |       |--- name : 'MaxGradient_G_per_cm'
%    |       |       |       |-- value : 2.4
%    |       |       |    
%    |       |       |--- userParameterDouble(2)
%    |       |       |       |    
%    |       |       |       |--- name : 'MaxSlewRate_G_per_cm_per_s'
%    |       |       |       |-- value : 14400
%    |       |       |    
%    |       |       |--- userParameterDouble(3)
%    |       |       |       |    
%    |       |       |       |--- name : 'FOVCoeff_1_cm'
%    |       |       |       |-- value : 20.8
%    |       |       |    
%    |       |       |--- userParameterDouble(4)
%    |       |       |       |    
%    |       |       |       |--- name : 'krmax_per_cm'
%    |       |       |       |-- value : 3.84615
%    |       |       |    
%    |       |       |--- userParameterLong(1)
%    |       |       |       |    
%    |       |       |       |--- name : 'interleaves'
%    |       |       |       |-- value : 8
%    |       |       |    
%    |       |       |--- userParameterLong(2)
%    |       |       |       |    
%    |       |       |       |--- name : 'fov_coefficients'
%    |       |       |       |-- value : 1
%    |       |       |    
%    |       |       |--- userParameterLong(3)
%    |       |       |       |    
%    |       |       |       |--- name : 'SamplingTime_ns'
%    |       |       |       |-- value : 4000
%    |       |       |    
%    |       |       |----------- identifier : 'HargreavesVDS2000'
%    |       |    
%    |       |------------- trajectory : 'spiral'
%    |       |-------- echoTrainLength : [ ]
%    |    
%    |--- <strong>experimentalConditions</strong>
%    |       |    
%    |       |-- H1resonanceFrequency_Hz : 2.36283e+07
%    |    
%    |---<strong>measurementInformation</strong>
%    |       |    
%    |       |--- measurementDependency(1)
%    |       |       |    
%    |       |       |-- dependencyType : 'SenMap'
%    |       |       |--- measurementID : '41624_33051224_62116429_138'
%    |       |    
%    |       |--- measurementDependency(2)
%    |       |       |    
%    |       |       |-- dependencyType : 'Noise'
%    |       |       |--- measurementID : '41624_33051224_62116429_138'
%    |       |    
%    |       |---- frameOfReferenceUID : '1.3.12.2.1107.5.2.18.41624.1.20200413142118295.0.0.4960'
%    |       |---------- measurementID : '41624_33051224_62116429_150'
%    |       |-------- patientPosition : 'HFS'
%    |       |----------- protocolName : 'ax_0mm_spiral_gmax24'
%    |    
%    |--- <strong>sequenceParameters</strong>
%    |       |    
%    |       |------------- TE : 1
%    |       |------------- TR : 40
%    |       |-- flipAngle_deg : 20
%    |    
%    |--- <strong>studyInformation</strong>
%    |       |    
%    |       |-- studyTime : '18:49:53'
%    |    
%    |--- <strong>userParameters</strong>
%    |       |    
%    |       |--- userParameterDouble(1)
%    |       |       |    
%    |       |       |--- name : 'MaxwellCoefficient_0'
%    |       |       |-- value : 0
%    |       |    
%    |       |--- userParameterDouble(2)
%    |       |       |    
%    |       |       |--- name : 'MaxwellCoefficient_1'
%    |       |       |-- value : 0
%    |       |    
%    |       |--- userParameterDouble(3)
%    |       |       |    
%    |       |       |--- name : 'MaxwellCoefficient_2'
%    |       |       |-- value : 0
%    |       |    
%    |       |--- userParameterDouble(4)
%    |       |       |    
%    |       |       |--- name : 'MaxwellCoefficient_3'
%    |       |       |-- value : 0
%    |       |    
%    |   
%  
help h5_hdr_structure
end