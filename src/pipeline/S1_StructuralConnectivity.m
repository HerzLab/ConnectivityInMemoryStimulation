% Description:
%   This script analyzes structural connectivity (SC) data as part of the study:
%   "Structural and Functional Connectivity Predict the Effects of Direct Brain Stimulation on Memory".
%   It characterizes the underlying structural scaffold of the brain by calculating
%   Mean Structural Connectivity (SC) and Coefficient of Variation (CV).
%   Crucially, it correlates nodal structural strength (degree) with memory recall
%   performance (delta_recall) across different patient groups (Close, All, Random).
% Inputs:
%   - Connectivity Data: *ifod2act_chacoconn_cocommp438subj_mean.csv in data/NeMo/
%   - Atlas: cocommp438_dil2_allsubj_mode.nii.gz in data/
%   - Clinical Data: subjects.xlsx (containing 'delta_recall', 'randomstim') in root/data/
% Outputs:
%   - Volume Maps (*.nii): Mean_StructuralConnectivity, CV_StructuralConnectivity,
%                          StructuralConnectivity_close_r (Correlation map for Close group).
%   - Surface Maps (*.tiff): Visualizations of degree and correlations.
%   - Statistics: Rs.mat containing correlation coefficients and p-values.
% External toolboxes documentation:
% (Please ensure these toolboxes are added to your MATLAB path using addpath)
%   - Enigma Toolbox: https://enigma-toolbox.readthedocs.io/en/latest/
%   - Permutools: https://github.com/mickcrosse/PERMUTOOLS
% Author: Qirui Zhang, Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
% Email: qirui.zhang@jefferson.edu; fmrizhangqr@126.com
% Date: 12/24/2025

clear; clc;

%% Setup Paths
% Determine project root (assumes src/pipeline/S1_GLM_WhiteMatter.m)
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
dataDir = fullfile(rootDir, 'data', 'NeMo');
outputDir = fullfile(rootDir, 'outputs', 'StructuralConnectivity');
nemoToHCPPath = fullfile(rootDir, 'data', 'Nemo2HCP360.xlsx');


% Check if output directory exists, if not create it
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Load Data
% Note: Atlas path seems to be in a sibling directory to the project root 'github'
% Adjusted to point to the likely location based on previous context
atlasPath = fullfile(rootDir, 'data', 'cocommp438_dil2_allsubj_mode.nii.gz');

fileList = dir(fullfile(dataDir, '*ifod2act_chacoconn_cocommp438subj_mean.csv'));

nSubjects = length(fileList);
% Initialize data containers
% We will initialize inside the loop or pre-allocate if dims are known.

for i = 1:nSubjects
    filePath = fullfile(fileList(i).folder, fileList(i).name);
    individualData = load(filePath);

    % Initialize arrays on first iteration to handle dynamic size
    if i == 1
        [nRegions, ~] = size(individualData);
        connectivityMatrixAll = zeros(nRegions, nRegions, nSubjects);
        connectivityDataAll = zeros(nRegions, nSubjects);
    end

    individualStructuralConnectivity = sum(individualData, 1);

    connectivityMatrixAll(:,:,i) = individualData;
    connectivityDataAll(:,i) = individualStructuralConnectivity'; % Transpose to match (regions x subjects)
end

% Handle zero values
positiveVals = connectivityDataAll(connectivityDataAll > 0);
minVal = min(positiveVals);

connectivityDataClean = connectivityDataAll;
connectivityDataClean(connectivityDataAll == 0) = 1.9e-07;
connectivityDataLog = log(connectivityDataClean);
connectivityDataLog10 = log10(connectivityDataClean);

%% Save Mean StructuralConnectivity and CV Maps
meanStructuralConnectivity = mean(connectivityDataAll, 2); % Mean across subjects (dim 2)
cvStructuralConnectivity = std(connectivityDataAll, 0, 2) ./ meanStructuralConnectivity;

parcel_to_volume(meanStructuralConnectivity, atlasPath, fullfile(outputDir, 'Mean_StructuralConnectivity.nii'), 1:438);
parcel_to_volume(cvStructuralConnectivity, atlasPath, fullfile(outputDir, 'CV_StructuralConnectivity.nii'), 1:438);
parcel_to_volume(connectivityDataAll(:,2), atlasPath, fullfile(outputDir, '02R1195E.nii'), 1:438);

meanStructuralConnectivityLog = log10(meanStructuralConnectivity);
parcel_to_volume(meanStructuralConnectivityLog, atlasPath, fullfile(outputDir, 'Mean_StructuralConnectivity_log.nii'), 1:438);

% Surface Plots
parcelValue = meanStructuralConnectivityLog;
outputName = fullfile(outputDir, 'Mean_StructuralConnectivity_log.tiff');
colorRange = [-2 2];
cmap = 'fake_parula';
parcel_to_HCPsurface(parcelValue, outputName, colorRange, cmap, nemoToHCPPath);

parcelValue = log(connectivityDataAll(:,2));
outputName = fullfile(outputDir, '02R1195E_log.tiff');
colorRange = [-4 5];
cmap = 'fake_parula';
parcel_to_HCPsurface(parcelValue, outputName, colorRange, cmap, nemoToHCPPath);


%% Correlation Analysis
% Load Clinical Data
subjectsFile = fullfile(rootDir, 'data', 'subjects.xlsx');
subjectsTable = readtable(subjectsFile);

% Extract variables using column names
% User specified: 'randomstim' for loop, 'delta_recall' for recall
loopVar = subjectsTable.randomstim;
recallVar = subjectsTable.delta_recall;

% Define Groups
groups = struct('name', {}, 'indices', {});
% All
groups(1).name = 'all';
groups(1).indices = true(nSubjects, 1);
% Close (randomstim == 0)
groups(2).name = 'close';
groups(2).indices = (loopVar == 0);
% Random (randomstim == 1)
groups(3).name = 'random';
groups(3).indices = (loopVar == 1);

% Initialize structure to store results for later use
analysisResults = struct();

for i = 1:length(groups)
    groupName = groups(i).name;
    idx = groups(i).indices;

    fprintf('Running correlation analysis for group: %s\n', groupName);

    % Prepare data for correlation
    % x: repeated recall scores for each region
    % Filter recall and connectivity data by group indices
    currentRecall = recallVar(idx);

    % Note: x and y must have same number of rows (nGroupSubjects)
    x = repmat(currentRecall, 1, 438);
    y = connectivityDataAll(:, idx)';

    % Perform Permutation Correlation
    [r, p, ~, ~, ~] = permucorr2(x, y, 'correct', 1, 'tail', 'both', 'type', 'kendall', 'verbose', 1);

    % Store in results structure
    analysisResults.(groupName).r = r;
    analysisResults.(groupName).p = p;

    % Filter significant results
    r_sig = r;
    r_sig(p > 0.05) = 0;

    % Construct filenames
    volNameR = sprintf('StructuralConnectivity_%s_r.nii', groupName);
    volNameSig = sprintf('Node_%s_r_permu.nii', groupName);
    surfNameR = sprintf('StructuralConnectivity_%s_r.tiff', groupName);
    surfNameSig = sprintf('StructuralConnectivity_%s_r_permu.tiff', groupName);

    % Save Volume Maps
    parcel_to_volume(r, atlasPath, fullfile(outputDir, volNameR), 1:438);
    parcel_to_volume(r_sig, atlasPath, fullfile(outputDir, volNameSig), 1:438);

    % Save Surface Maps
    colorRange = [-0.5 0.5];
    cmap = 'RdBu_r0';

    parcel_to_HCPsurface(r, fullfile(outputDir, surfNameR), colorRange, cmap, nemoToHCPPath);
    parcel_to_HCPsurface(r_sig, fullfile(outputDir, surfNameSig), colorRange, cmap, nemoToHCPPath);
end



%% Save Results
save(fullfile(outputDir, 'analysisResults.mat'), 'analysisResults');

%% Feature Selection
% Top All Features
r_all(isnan(r_all)) = 0;
[~, sortedIdx] = sort(r_all, 'descend');
topIdx = sortedIdx(1:10);

scTop10FeatureLog = connectivityDataLog(topIdx, :)';
scTop10Feature = connectivityDataAll(topIdx, :)';

% Top Close Features
r_close(isnan(r_close)) = 0;
[~, sortedIdx] = sort(r_close, 'descend');
topIdx = sortedIdx(1:10);

scTop10FeatureLog_Close = connectivityDataLog(topIdx, :)';
scTop10Feature_Close = connectivityDataAll(topIdx, :)';
