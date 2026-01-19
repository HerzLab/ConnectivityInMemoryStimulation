function [Network_value] = mean_by_atlas_4d(Volume, Atlas, atlas_order)
% mean_by_atlas_4d: Extracts average time series from each ROI defined by Atlas
%
% Inputs:
%   Volume:      4D NIfTI file path or 4D numeric array [x,y,z,t]
%   Atlas:       3D NIfTI file path or 3D numeric array [x,y,z]
%   atlas_order: vector of ROI labels to extract
%
% Output:
%   Network_value: [timepoints x number of ROIs]

    % --- Load Volume (4D) ---
    if ischar(Volume) || isstring(Volume)
        V_info = niftiinfo(Volume);
        V = niftiread(V_info);  % [X Y Z T]
    elseif isnumeric(Volume) && ndims(Volume) == 4
        V = Volume;
    else
        error('Volume must be a 4D numeric array or a valid NIfTI file path.');
    end

    [nx, ny, nz, nt] = size(V);
    V_2d = reshape(V, nx*ny*nz, nt);  % [voxels x time]

    % --- Load Atlas (3D) ---
    if ischar(Atlas) || isstring(Atlas)
        A = niftiread(Atlas);  % [X Y Z]
    elseif isnumeric(Atlas) && isequal(size(Atlas), [nx, ny, nz])
        A = Atlas;
    else
        error('Atlas must be a 3D numeric array matching Volume size, or a valid NIfTI path.');
    end

    A_2d = reshape(A, nx*ny*nz, 1);  % [voxels x 1]

    % --- Preallocate output ---
    n_rois = length(atlas_order);
    Network_value = nan(nt, n_rois);  % [time x region]

    % --- Extract time series for each ROI ---
    for i = 1:n_rois
        roi_label = atlas_order(i);
        roi_index = (A_2d == roi_label);

        if any(roi_index)
            roi_data = V_2d(roi_index, :);  % [n_voxels x time]
            Network_value(:, i) = mean(roi_data, 1, 'omitnan');
        else
            Network_value(:, i) = NaN;
        end
    end
end
