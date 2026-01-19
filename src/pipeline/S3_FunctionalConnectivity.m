% Description:
%   This script identifies functional features that drive memory performance using
%   Resting-state and Task-based (Word List) connectivity.
%   It starts by creating normative templates (averaging FC maps across healthy controls)
%   for both Rest and WL conditions.
%   It then correlates individual patient Functional Connectivity (FC) patterns with
%   recall performance (delta_recall).
%   Finally, it performs feature selection to identify the top predictive edges/nodes
%   in the 'Close' stimulation group for downstream analysis.
% Inputs:
%   - Atlas: cocommp438_dil1_allsubj_mode_2mm.nii in data/
%   - Brain Mask: cocommp438_dil3_allsubj_mode_2mm.nii in data/
%   - MNI ROIs: MNI_ROI in data/
%   - Clinical Data: subjects.xlsx in data/
%   - Pre-computed ROI Signals: HC_Parcels_ROI_signal.mat（Pre-computed ROI timeseries for HC parcels）
%                               HC_electrode_signal.mat (Pre-computed electrode timeseries for HC electrodes)
% Outputs:
%   - FC Maps: Mean_Rest_FC.tiff, Mean_WL_FC.tiff (Normative FC maps)
%   - Correlation Maps: [Condition]_[Group]_r.tiff (Correlation maps)
%   - Top Features: Func_feature.mat (Top 10 weighted features/parcels from Close group)
% External toolboxes documentation:
% (Please ensure these toolboxes are added to your MATLAB path using addpath)
%   - Enigma Toolbox: https://enigma-toolbox.readthedocs.io/en/latest/
%   - Permutools: https://github.com/mickcrosse/PERMUTOOLS
% Author: Qirui Zhang, Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
% Email: qirui.zhang@jefferson.edu; fmrizhangqr@126.com
% Date: 12/24/2025

clear; clc;

%% Setup Paths
% Determine project root (assumes src/pipeline/S3_FunctionalConnectivity.m)
this_file = mfilename('fullpath');

if ~isempty(this_file)
    % Case 1: script/function is being executed from a .m file
    project_root = fileparts(fileparts(fileparts(this_file)));
else
    % Case 2: executed from Command Window or interactive context
    project_root = fileparts(fileparts(pwd));
end

rootDir = project_root;
addpath(genpath(fullfile(rootDir, 'src', 'utils')));

%% Define Directories
% Using relative paths based on the project structure
dataDir = fullfile(rootDir, 'data');
outputDir = fullfile(rootDir, 'outputs', 'FunctionalConnectivity');
nemoToHCPPath = fullfile(dataDir, 'Nemo2HCP360.xlsx');

% (Paths for XCP and HC list removed as signals are loaded directly)

% ROI Directory (Requested to be dataDir/MNI_ROI)
roiDir = fullfile(dataDir, 'MNI_ROI');

% Check if output directory exists, if not create it
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Load Atlas and Mask
atlasPath = fullfile(dataDir, 'cocommp438_dil1_allsubj_mode_2mm.nii');
brainMaskPath = fullfile(dataDir, 'cocommp438_dil3_allsubj_mode_2mm.nii');



%% Load Pre-computed Signals
fprintf('Loading HC Parcels ROI Signal and Electrode Signal...\n');
if exist(fullfile(dataDir, 'HC_Parcels_ROI_signal.mat'), 'file') && ...
        exist(fullfile(dataDir, 'HC_electrode_signal.mat'), 'file')
    load(fullfile(dataDir, 'HC_Parcels_ROI_signal.mat'), 'WL_signal', 'Rest_signal');
    load(fullfile(dataDir, 'HC_electrode_signal.mat'), 'ROI_WL_signal', 'ROI_Rest_signal');
else
    error('Pre-computed signal files not found in data directory.');
end

% Infer dimensions from loaded data
nHC = size(WL_signal, 3);
nSubjects = size(ROI_WL_signal, 3);

%% Load Clinical Data
fprintf('Loading Clinical Data...\n');
subjectsFile = fullfile(dataDir, 'subjects.xlsx');
subjectsTable = readtable(subjectsFile);

loopVar = subjectsTable.randomstim;
recallVar = subjectsTable.delta_recall;

%% 3. FC Calculation
fprintf('Calculating Functional Connectivity...\n');
% Signals are loaded above

% Pre-allocate
rRest = zeros(438, nHC, nSubjects);
pRest = zeros(438, nHC, nSubjects);
rWL = zeros(438, nHC, nSubjects);
pWL = zeros(438, nHC, nSubjects);

for i = 1:nSubjects
    for k = 1:nHC
        % Electrode Signal [Time x 1]
        roiRestIndi = ROI_Rest_signal(:, k, i);
        roiWlIndi = ROI_WL_signal(:, k, i);

        % Atlas Signal [Time x 438]
        restAtlasIndi = Rest_signal(:, :, k);
        wlAtlasIndi = WL_signal(:, :, k);

        [rR, pR] = corr(roiRestIndi, restAtlasIndi);
        [rW, pW] = corr(roiWlIndi, wlAtlasIndi);

        rRest(:, k, i) = rR';
        pRest(:, k, i) = pR';

        rWL(:, k, i) = rW';
        pWL(:, k, i) = pW';
    end
end

% Average across HCs (dim 2)
% Result: [438 x nSubjects]
rRestMean = squeeze(mean(rRest, 2, 'omitnan'));
rWlMean = squeeze(mean(rWL, 2, 'omitnan'));

% Save Mean FC Maps (averaged across subjects)
parcelValue = mean(rRestMean, 2);
outputName = fullfile(outputDir, 'Mean_Rest_FC.tiff');
colorRange = [-0.6 0.6];
cmap = 'fake_parula';
parcel_to_HCPsurface(parcelValue, outputName, colorRange, cmap, nemoToHCPPath);

parcelValue = mean(rWlMean, 2);
outputName = fullfile(outputDir, 'Mean_WL_FC.tiff');
parcel_to_HCPsurface(parcelValue, outputName, colorRange, cmap, nemoToHCPPath);

% 02R1195E (Index 2 subject)
parcelValue = rWlMean(:, 2);
outputName = fullfile(outputDir, '02R1195E.tiff');
parcel_to_HCPsurface(parcelValue, outputName, colorRange, cmap, nemoToHCPPath);

%% 4. Correlation Analysis
fprintf('Running Correlation Analysis with Recall...\n');

groups = struct('name', {}, 'idx', {}, 'recall', {});
groups(1).name = 'close';
groups(1).idx = find(loopVar == 0);
groups(1).recall = recallVar(groups(1).idx);

groups(2).name = 'all';
groups(2).idx = (1:length(loopVar))';
groups(2).recall = recallVar;

groups(3).name = 'random';
groups(3).idx = find(loopVar == 1);
groups(3).recall = recallVar(groups(3).idx);

means = struct('label', {'WL', 'Rest'}, 'mat', {rWlMean, rRestMean});

colorRange = [-0.5 0.5];
cmap = 'RdBu_r0';
alphaThr = 0.05;

% Initialize Stats structure
Stats = struct();

for i = 1:numel(means)
    meanLabel = means(i).label;
    fcMatrix = means(i).mat; % [438 x nSub]

    for j = 1:numel(groups)
        grpLabel = groups(j).name;
        sel = groups(j).idx;
        recallV = groups(j).recall;

        % Data for Permutation
        x = repmat(recallV, 1, 438);
        y = fcMatrix(:, sel)'; % [nSel x 438]

        % Permutation Correlation
        [rVec, pPerm, ~, ~, ~] = permucorr2(x, y, 'correct', 1, 'tail', 'both', 'type', 'pearson');

        % Uncorrected Correlation
        [rUncFull, pUncFull] = corr(recallV, y, 'Type', 'Pearson', 'Rows', 'pairwise');
        pUnc = pUncFull;

        % Store in Stats structure
        Stats.(meanLabel).(grpLabel).r = rVec;
        Stats.(meanLabel).(grpLabel).p_unc = pUnc;
        Stats.(meanLabel).(grpLabel).p_perm = pPerm;

        % Filename Base
        base = sprintf('%s_%s', meanLabel, grpLabel);

        % Save Unthresholded
        parcel_to_volume(rVec, atlasPath, fullfile(outputDir, [base, '_r.nii']), 1:438);
        parcel_to_HCPsurface(rVec, fullfile(outputDir, [base, '_r.tiff']), colorRange, cmap, nemoToHCPPath);

        % Save Uncorrected Significant (p < 0.05)
        rUncorrMask = rVec;
        rUncorrMask(pUnc > 0.05) = 0;
        parcel_to_volume(rUncorrMask, atlasPath, fullfile(outputDir, [base, '_r_p0p05_uncorr.nii']), 1:438);
        parcel_to_HCPsurface(rUncorrMask, fullfile(outputDir, [base, '_r_p0p05_uncorr.tiff']), colorRange, cmap, nemoToHCPPath);

        % Save Corrected Significant (p_perm < 0.05)
        rPermMask = rVec;
        rPermMask(pPerm > alphaThr) = 0;
        parcel_to_volume(rPermMask, atlasPath, fullfile(outputDir, [base, '_r_permu.nii']), 1:438);
        parcel_to_HCPsurface(rPermMask, fullfile(outputDir, [base, '_r_permu.tiff']), colorRange, cmap, nemoToHCPPath);

        % Save Stats Table
        T = table((1:438)', rVec(:), pUnc(:), pPerm(:), ...
            'VariableNames', {'Parcel', 'r', 'p_uncorrected', 'p_perm_corrected'});
        writetable(T, fullfile(outputDir, [base, '_r_p_table.csv']));

        fprintf('[Done] %s-%s: nSub=%d -> %s\n', meanLabel, grpLabel, length(recallV), base);
    end
end
close all
% Save Stats structure
save(fullfile(outputDir, 'FunctionalConnectivity_Stats.mat'), 'Stats');

%% Feature Selection
fprintf('Selecting Top Features from Close group results...\n');

% WL - Close
rWlClose = Stats.WL.close.r;

% Rest - Close
rRestClose = Stats.Rest.close.r;

% Top 10 WL (Signed sort as per original logic)
[~, sortedIdxWl] = sort(rWlClose, 'descend');
topIdxWl = sortedIdxWl(1:10);
wlTop10Feature = rWlMean(topIdxWl, :)';

% Top 10 Rest (Absolute sort as per original logic)
[~, sortedIdxRest] = sort(abs(rRestClose), 'descend');
topIdxRest = sortedIdxRest(1:10);
restTop10Feature = rRestMean(topIdxRest, :)';

% Save to outputDir
save(fullfile(outputDir, 'Func_feature.mat'), 'wlTop10Feature', 'restTop10Feature');
