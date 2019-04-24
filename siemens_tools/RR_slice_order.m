function slice_order = RR_slice_order(numSlices, subset)
% slice_order = RR_slice_order(numSlices, <subset>)
% slice_order [vec] or with optional subset - single number or 3D-stack
% numSlices [number or 3D-stack], subset [number] which correctly ordered slice you
% want
% Returns the acquired order of interleaves slices. 
% subset [range] can further choose a selection of slices
% Empirically derived from Siemens interleaved slice ordering
% R Ramasawmy NHLBI Nov 2018

if length(size(numSlices)) == 3
    % rearrange 3D stack
    
    if nargin < 2
        slice_order = numSlices(:,:,RR_slice_order(size(numSlices,3) ));
    else
        slice_order = numSlices(:,:,RR_slice_order(size(numSlices,3), subset ));
    end
else
    % empiricial determination of order 
    
    slice_order = zeros(1, numSlices);
    
    if mod(numSlices,2)
        % odd number of slices
        vec = 1:round(numSlices/2);
        slice_order(1:2:end) = vec ;
        slice_order(2:2:end) = vec(1:end-1)+ round(numSlices/2);
    else
        % even number of slices
        vec = 1:numSlices/2;
        slice_order(1:2:end) = vec + numSlices/2;
        slice_order(2:2:end) = vec;
    end
    
    if nargin > 1
        slice_order = slice_order(subset);
    end
end

end