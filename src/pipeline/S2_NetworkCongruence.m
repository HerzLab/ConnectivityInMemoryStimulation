% Description:
%   This script calculates the Dice coefficient between the Word List (WL) functional
%   activation map and structural connectivity networks. It serves to quantify
%   the spatial "match" (congruence) between function and structure.
%   Hypothesis: Higher congruence between the functional network and the structural
%   scaffold predicts better memory recall performance.
%   Analysis:
%       1. Define WL Network from activation map (Top K regions).
%       2. Define SC Network seeded from WL nodes (Top connected regions).
%       3. Compute Dice overlap and correlate with recall.
% Inputs:
%   - activation Map: WL_Active-rest_HC_1mm.nii (Average functional activation map for WL task) in data/
%   - Atlas: cocommp438_dil2_allsubj_mode.nii.gz in data/
%   - Clinical Data: subjects.xlsx (with 'randomstim', 'delta_recall') in data/
%   - Region Values: Region_values.xlsx (with 'heatmap') in data/
%   - Connectivity Data: *ifod2act_chacoconn_cocommp438subj_mean.csv in data/NeMo/
% Outputs:
%   - Network Maps: WL_node_map.nii (Functional Network), WL_node_map.tiff
%   - Correlation Plots: Dice_r_close.tiff (Matrix plot of Congruence vs Recall correlations)

% Author: Qirui Zhang, Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
% Email: qirui.zhang@jefferson.edu; fmrizhangqr@126.com
% Date: 12/24/2025

clear; clc;

%% Setup Paths
% Determine project root (assumes src/pipeline/S2_NetworkCongruence.m)
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
nemoDir = fullfile(dataDir, 'NeMo');
outputDir = fullfile(rootDir, 'outputs', 'NetworkCongruence');
nemoToHCPPath = fullfile(dataDir, 'Nemo2HCP360.xlsx');

% Check if output directory exists, if not create it
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Load Input Files
volumeFile = fullfile(dataDir, 'WL_Active-rest_HC_1mm.nii');
atlasPath = fullfile(dataDir, 'cocommp438_dil2_allsubj_mode.nii.gz'); % Matches S1 path
subjectsFile = fullfile(dataDir, 'subjects.xlsx');
regionValuesFile = fullfile(dataDir, 'Region_values.xlsx');

%% Load Clinical Data
subjectsTable = readtable(subjectsFile);
loopVar = subjectsTable.randomstim;
recallVar = subjectsTable.delta_recall;

% Define Groups (Indices)
allIdx = true(length(loopVar), 1);
closeIdx = (loopVar == 0);
randomIdx = (loopVar == 1);

recallAll = recallVar;
recallClose = recallVar(closeIdx);
recallRandom = recallVar(randomIdx);

%% Load Region Values (Heatmap)
% User requested readtable with header 'heatmap'
regionTable = readtable(regionValuesFile);
heatmapValues = regionTable.heatmap;

% Determine excluded indices (where heatmap > 0)
outHeatmapIndex = setdiff(1:438, find(heatmapValues > 0));

%% Process Volume and Atlas
atlasOrder = 1:438;
% Calculate mean activation per region
% Calculate mean activation per region
networkValue = mean_by_atlas(volumeFile, atlasPath, atlasOrder);

% Save Network Map Volume
parcel_to_volume(networkValue, atlasPath, fullfile(outputDir, 'WL_node_map.nii'), 1:438);

% Save Network Map Surface
parcelValue = networkValue;
outputName = fullfile(outputDir, 'WL_node_map.tiff');
colorRange = [0 6];
cmap = 'Reds';
parcel_to_HCPsurface(parcelValue, outputName, colorRange, cmap, nemoToHCPPath);

% Thresholded Map
networkValue05 = networkValue;
networkValue05(networkValue < 1.6) = 0;
outputName = fullfile(outputDir, 'WL_node_map05.tiff');
parcel_to_HCPsurface(networkValue05, outputName, colorRange, cmap, nemoToHCPPath);

%% Load Connectivity Data
fileList = dir(fullfile(nemoDir, '*ifod2act_chacoconn_cocommp438subj_mean.csv'));
nSubjects = length(fileList);

% Pre-allocate
nRegions = 438; % Known dimension
connectivityDataAll = zeros(nRegions, nSubjects);
connectivityMatrixAll = zeros(nRegions, nRegions, nSubjects);

for i = 1:nSubjects
    filePath = fullfile(fileList(i).folder, fileList(i).name);
    individualData = load(filePath);

    individualDegree = sum(individualData, 1); % Sum columns

    connectivityMatrixAll(:,:,i) = individualData;
    connectivityDataAll(:,i) = individualDegree';
end

%% Plot Example Subject (02R1195E - Index 2)
exampleSubjectIdx = 2;
if nSubjects >= exampleSubjectIdx
    parcelValue = connectivityDataAll(:, exampleSubjectIdx);
    outputName = fullfile(outputDir, '02R1195E.tiff');
    colorRange = [0 50];
    cmap = 'Reds';
    parcel_to_HCPsurface(parcelValue, outputName, colorRange, cmap, nemoToHCPPath);

    % Binary Plot
    val = parcelValue;
    val(val <= 1) = 0;
    val(val > 1) = 1;
    outputName = fullfile(outputDir, '02R1195E_T1.tiff');
    parcel_to_HCPsurface(val, outputName, [0 1.3], 'fake_parula', nemoToHCPPath);
end

%% Dice Analysis
positiveRegionNumber = sum(networkValue > 0);

wlLevels = 0.1:0.1:1;
wlLevelsName = {'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'};
connLevels = [0.01, 0.1, 1, 10, 50];
connLevelsName = {'0.01','0.1','1','10','50'};

% Pre-allocate Result Matrices
rCloseAll = zeros(length(connLevels), length(wlLevels));
pCloseAll = zeros(length(connLevels), length(wlLevels));
rRandomAll = zeros(length(connLevels), length(wlLevels));
pRandomAll = zeros(length(connLevels), length(wlLevels));
rAllAll = zeros(length(connLevels), length(wlLevels));
pAllAll = zeros(length(connLevels), length(wlLevels));

% Matching Matrix: (Subjects x Thresholds x Levels)
matchingMatrix = zeros(nSubjects, length(connLevels), length(wlLevels));

for k = 1:length(wlLevels)
    % Select Top K regions from Activation Map
    topK = round(positiveRegionNumber * wlLevels(k));
    [maxVals, ~] = maxk(networkValue, topK);
    minTopVal = min(maxVals);
    wlMask = (networkValue >= minTopVal); % Regions in activation network

    for j = 1:length(connLevels)
        for i = 1:nSubjects
            indi = connectivityDataAll(:, i);
            % Threshold connectivity degree
            wmMask = (indi > connLevels(j));

            % Compute Dice Coefficient
            matchingMatrix(i, j, k) = dice(wlMask', wmMask);
        end
    end

    % Correlate Dice with Recall for this WL Threshold (k), across all Conn Thresholds (:, k)
    % Results for THIS wlLevel (k), across all levels (rows of result)

    [r, p] = corr(matchingMatrix(closeIdx, :, k), recallClose, 'Type', 'Spearman');
    rCloseAll(:, k) = r'; % Store as column (nConnLevels x 1) for this k
    pCloseAll(:, k) = p';

    [r1, p1] = corr(matchingMatrix(randomIdx, :, k), recallRandom, 'Type', 'Spearman');
    rRandomAll(:, k) = r1';
    pRandomAll(:, k) = p1';

    [r2, p2] = corr(matchingMatrix(allIdx, :, k), recallAll, 'Type', 'Spearman');
    rAllAll(:, k) = r2';
    pAllAll(:, k) = p2';
end

%% Plotting Results (FDR Correction)
figureFiles = {'Dice_r_random.tiff', 'Dice_r_close.tiff', 'Dice_r_all.tiff'};
rMatrices = {rRandomAll, rCloseAll, rAllAll};
pMatrices = {pRandomAll, pCloseAll, pAllAll};
alpha = 0.05;

for k = 1:3
    rFig = rMatrices{k};
    pFig = pMatrices{k};

    rFig(isnan(rFig)) = 0;
    pFig(isnan(pFig)) = 1;

    % Benjamini-Hochberg FDR
    pVec = pFig(:);
    qVec = mafdr(pVec, 'BHFDR', true);
    qFig = reshape(qVec, size(pFig));

    % Plot
    f = figure('Visible', 'off'); % Don't display popups
    imagesc(rFig, [-0.6 0.6]);
    c = RdBu_r;
    colormap(c);
    colorbar('southoutside');

    % Axis Labels
    yticks(1:length(connLevels));
    xticks(1:length(wlLevels));
    xticklabels(wlLevelsName);
    yticklabels(connLevelsName);

    set(gca, 'XAxisLocation', 'top');
    set(gca, 'FontSize', 18);
    set(gcf, 'Color', 'w', 'Position', [200, 200, 700, 600]);

    % Add Significance Markers
    hold on;
    for i = 1:length(connLevels)
        for j = 1:length(wlLevels)
            if qFig(i,j) < alpha
                % ** = FDR Significant
                text(j, i, '**', 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'FontSize', 18, ...
                    'Color', 'w', 'FontWeight', 'bold');
            elseif pFig(i,j) < 0.05
                % * = Uncorrected Significant
                text(j, i, '*', 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'FontSize', 18, ...
                    'Color', 'w', 'FontWeight', 'bold');
            end
        end
    end
    hold off;

    print(gcf, '-dtiff', '-r300', fullfile(outputDir, figureFiles{k}));
    close(f);
end
