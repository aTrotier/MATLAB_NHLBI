function ismrmrd_s = xml2hdr(xml_string)
% ismrmrd_s = xml2hdr(xml_string)
%
% Method to extract scan parameters
% xml_string = h5read(h5FILE, '/dataset/xml'); // ismrmrd_s = read_h5_header(dfile)
% Non-intelligent way of parsing duplicate labelling, just pads it with
% _NUM (i.e. coilLabel).
% The extract method can probably be vastly improved..
%
% % ismrmrdHeader \\
%     studyInformation
%     measurementInformation
%     acquisitionSystemInformation
%     experimetnatlConditions
%     encoding
%     sequenceParameters
%     userParameters
%
% +++ Quick index: +++
% % DIMENSIONS
% ismrmrd_s.encoding.reconSpace.fieldOfView_mm.<x y z>
% ismrmrd_s.encoding.reconSpace.matrixSize.<x y z>
% ismrmrd_s.encoding.encodedSpace.fieldOfView_mm.<x y z>
% ismrmrd_s.encoding.encodedSpace.matrixSize.<x y z>
%
% % CONTRAST
% ismrmrd_s.sequenceParameters.TR
% ismrmrd_s.sequenceParameters.TE
% ismrmrd_s.sequenceParameters.TI
% ++++++
%
% R Ramasawmy NHLBI Nov 2018


%%
if iscell(xml_string)
    xml_string = xml_string{1};
end

%     subfields = { ...
%         'studyInformation',  ...
%         'measurementInformation', ...
%         'acquisitionSystemInformation', ...
%         'experimetnatlConditions', ...
%         'encoding', ...
%         'sequenceParameters', ...
%         'userParameters'}; % can potentially use this with regexp for future implementations

run_legacy_method = 0; % Older version, until Adrienne broke it! New method is theoretically more stable (N=2 tests)

if run_legacy_method
    a = strsplit(xml_string);
    start_val = 1;
else
    a = strsplit(xml_string, '\n'); % using new line to better separate parameters
    start_val = 3;
end

    ismrmrd_s = struct();
    current_subfield = [];

    for i = start_val:length(a)-2 % {</ismrmrdHeader>, ''}, can alternatively look for /ismrmrdHeader to finish..

        if run_legacy_method
            [xml_code,value]=xml_interp_RR_old( a{i} );
        else
            [xml_code,value]=xml_interp_RR( a{i} );
        end

        switch xml_code
            case 1
                % start section
                if isempty(current_subfield) % parent section ("subfield" level)
                    current_subfield{1} = value;
                    fields{1} = value;
                else
                    current_field = value;
                    if strcmp(current_field, last_field) % check for dubplicate label fields
                        repeat_counter = repeat_counter+1;
                        fields = {current_subfield{:} [value '_' num2str(repeat_counter)]};
                    else
                        repeat_counter = 0;
                        fields = {fields{:} value};
                    end
                end
                last_field = value;

            case 2
                % parameter value
                fields_write = {fields{:} value.name};
                ismrmrd_s = setfield(ismrmrd_s, fields_write{:}, value.value);

            case 3
                % end section

                if strcmp(current_subfield, value)
                    % close subfield
                    current_subfield = [];
                    fields = {};
                else
                    % step back
                    fields = fields(1,1:length(fields)-1);
                end

        end

    end

end
function [xml_code, value] = xml_interp_RR(str)
% xml_code = 1 (start section), 2 (parameter), 3 (end section)
% 0s for header information, which has current format:
test1 = regexp(str,'>');
test2 = regexp(str,'<');

if isempty(test1) || isempty(test2)
    xml_code = 0; value = [];
else

    if length(test1) > 1
        % +++ parameter +++
        xml_code = 2;
        value1 = str(( test1(1)+1 ):( test2(2)-1 ));
        test3 = isstrprop(value1, 'alpha');
        test4 = isstrprop(value1, 'punct');
        % determine if punctuation is "decimal point"
        if length(find(test4)) == 1
            if regexp(value1(test4), '.')
                test4 = zeros(size(test4));
            end
        end
        % convert to value if necessary
        if sum(test3 + test4) == 0
            value1 = str2double(value1);
        end
        value.name  = str(test2(1)+1:test1(1)-1);
        value.value = value1;
    else
        test3 = regexp(str, '</');
        if test3
            % +++ end section +++
            xml_code = 3;
            value = str((test3+2):test1-1);
        else
            % +++ start section +++
            xml_code = 1;
            value = str((test2+1):test1-1);
        end
    end
end
end
function [xml_code, value] = xml_interp_RR_old(str)
% xml_code = 1 (start section), 2 (parameter), 3 (end section)
% 0s for header information, which has current format:
% %
% %     '<?xml'    'version="1.0"?>'    '<ismrmrdHeader'    'xmlns="http://www�'
% %     'xmlns:xsi="http:/�'    'xmlns:xs="http://�'    'xsi:schemaLocatio�'
% %      'ismrmrd.xsd">'

test1 = regexp(str,'>');
test2 = regexp(str,'<');

if isempty(test1) || isempty(test2)
    xml_code = 0; value = [];
else

    if length(test1) > 1
        % +++ parameter +++
        xml_code = 2;
        value1 = str(( test1(1)+1 ):( test2(2)-1 ));
        test3 = isstrprop(value1, 'alpha');
        test4 = isstrprop(value1, 'punct');
        % determine if punctuation is "decimal point"
        if length(find(test4)) == 1
            if regexp(value1(test4), '.')
                test4 = zeros(size(test4));
            end
        end
        % convert to value if necessary
        if sum(test3 + test4) == 0
            value1 = str2double(value1);
        end
        value.name  = str(test2(1)+1:test1(1)-1);
        value.value = value1;
    else
        test2 = regexp(str(2), '/');
        if test2
            % +++ end section +++
            xml_code = 3;
            value = str(3:test1-1);
        else
            % +++ start section +++
            xml_code = 1;
            value = str(2:test1-1);
        end
    end
end

end
