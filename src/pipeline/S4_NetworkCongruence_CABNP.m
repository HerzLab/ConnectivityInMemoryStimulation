% S4_NetworkCongruence_CABNP.m
% Description:
%   This script calculates the Dice coefficient between Cole-Anticevic (CAB-NP)
%   networks and structural connectivity (SC) networks.
%   Network definitions are loaded from Nemo2HCP360.xlsx.
%
% Author: Qirui Zhang, Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
% Date: 05/02/2026

clear; clc;

%% Setup Paths
this_file = mfilename('fullpath');
if ~isempty(this_file)
    project_root = fileparts(fileparts(fileparts(this_file)));
else
    project_root = fileparts(fileparts(pwd));
end
rootDir = project_root;
addpath(genpath(fullfile(rootDir, 'src', 'utils')));

%% Define Directories
dataDir   = fullfile(rootDir, 'data');
nemoDir   = fullfile(dataDir, 'NeMo');
outputDir = fullfile(rootDir, 'outputs', 'NetworkCongruence_CABNP');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Define Input Files
volumeFile    = fullfile(dataDir, 'WL_Active-rest_HC_1mm.nii');
atlasPath     = fullfile(dataDir, 'cocommp438_dil2_allsubj_mode.nii.gz');
subjectsFile  = fullfile(dataDir, 'subjects.xlsx');
nemoToHCPPath = fullfile(dataDir, 'Nemo2HCP360.xlsx');

% Surface plots subdirectory
surfaceDir = fullfile(outputDir, 'surfaces');
if ~exist(surfaceDir, 'dir')
    mkdir(surfaceDir);
end

%% Load Clinical Data
subjectsTable = readtable(subjectsFile);
loopVar   = subjectsTable.randomstim;
recallVar = subjectsTable.delta_recall;

% Define Groups
allIdx    = true(length(loopVar), 1);
closeIdx  = (loopVar == 0);
randomIdx = (loopVar == 1);

recallAll    = recallVar;
recallClose  = recallVar(closeIdx);
recallRandom = recallVar(randomIdx);

%% Generate Reference Network (Top 30% Verbal Memory)
fprintf('Generating Reference Network (Verbal Memory)... \n');
atlasOrder = 1:438;
refValue = mean_by_atlas(volumeFile, atlasPath, atlasOrder);
posCount = sum(refValue > 0);
topK = round(posCount * 0.3);
[maxVals, ~] = maxk(refValue, topK);
refMask = (refValue >= min(maxVals));
fprintf('  Reference mask created with %d regions.\n', sum(refMask));

%% Load Structural Connectivity Data
fprintf('Loading Structural Connectivity (SC) data...\n');
fileList   = dir(fullfile(nemoDir, '*ifod2act_chacoconn_cocommp438subj_mean.csv'));
nSubjects  = length(fileList);
nRegions   = 438;
connectivityDataAll = zeros(nRegions, nSubjects);

for i = 1:nSubjects
    filePath         = fullfile(fileList(i).folder, fileList(i).name);
    individualData   = load(filePath);
    individualDegree = sum(individualData, 1);
    connectivityDataAll(:, i) = individualDegree';
end
fprintf('  Loaded %d subjects.\n', nSubjects);

%% Load Network Definitions from Excel
% Col 1: ROI ID, Col 5: Network Name
opts = detectImportOptions(nemoToHCPPath);
opts.VariableNamingRule = 'preserve';
mappingTable = readtable(nemoToHCPPath, opts);

% Identify columns by position if names are missing
roiCol = 1;
netCol = 5;

roiIds = mappingTable{:, roiCol};
netNamesRaw = mappingTable{:, netCol};

% Find unique networks (ignoring NaNs)
uniqueNets = unique(netNamesRaw);
% If it's a cell array of strings, remove empty or 'NaN'
if iscell(uniqueNets)
    uniqueNets(cellfun(@isempty, uniqueNets)) = [];
    uniqueNets(strcmp(uniqueNets, 'NaN')) = [];
end
nNetworks = length(uniqueNets);

%% Loop Over CAB-NP Networks
rCloseVec  = nan(nNetworks, 1);
pCloseVec  = nan(nNetworks, 1);
rRandomVec = nan(nNetworks, 1);
pRandomVec = nan(nNetworks, 1);
rAllVec    = nan(nNetworks, 1);
pAllVec    = nan(nNetworks, 1);
diceWithRef = nan(nNetworks, 1); % Overlap with verbal memory network
nPositiveRegions = zeros(nNetworks, 1);
connThreshold = 1;

for n = 1:nNetworks
    netName = uniqueNets{n};
    fprintf('[%d/%d] Processing: %s ...\n', n, nNetworks, netName);
    
    % Step 1: Create mask for this network
    % Find all ROI IDs belonging to this network
    targetRoiIds = roiIds(strcmp(netNamesRaw, netName));
    
    nsMask = false(nRegions, 1);
    for r = 1:length(targetRoiIds)
        id = targetRoiIds(r);
        if ~isnan(id) && id > 0 && id <= nRegions
            nsMask(id) = true;
        end
    end
    nPositiveRegions(n) = sum(nsMask);
    
    if nPositiveRegions(n) == 0
        fprintf('  WARNING: No regions for %s. Skipping.\n', netName);
        continue;
    end
    
    % Step 1b: Dice with Reference (Verbal Memory Network)
    diceWithRef(n) = dice(nsMask(:), refMask(:));
    
    % Step 1c: Surface visualization
    safeNetName = regexprep(netName, '[^a-zA-Z0-9]', '_');
    surfOutputName = fullfile(surfaceDir, [safeNetName, '.tiff']);
    try
        parcel_to_HCPsurface(double(nsMask), surfOutputName, [0 1.5], 'Reds', nemoToHCPPath);
        close(gcf);
    catch ME
        fprintf('  Surface plot error: %s\n', ME.message);
    end
    
    % Step 2: Dice per subject
    diceScores = zeros(nSubjects, 1);
    for i = 1:nSubjects
        wmMask = (connectivityDataAll(:, i) > connThreshold);
        diceScores(i) = dice(nsMask(:), wmMask(:));
    end
    
    % Step 3: Correlate with delta_recall
    [rC, pC] = corr(diceScores(closeIdx), recallClose, 'Type', 'Spearman');
    [rR, pR] = corr(diceScores(randomIdx), recallRandom, 'Type', 'Spearman');
    [rA, pA] = corr(diceScores(allIdx), recallAll, 'Type', 'Spearman');
    
    rCloseVec(n) = rC;  pCloseVec(n) = pC;
    rRandomVec(n) = rR;  pRandomVec(n) = pR;
    rAllVec(n) = rA;  pAllVec(n) = pA;
    
    fprintf('  N_regions=%d | R_Close=%.3f (p=%.3f)\n', nPositiveRegions(n), rC, pC);
    
    % Scatter Plot
    fScatter = figure('Visible', 'off');
    set(fScatter, 'Color', 'w', 'Position', [200, 200, 500, 450]);
    diceClose = diceScores(closeIdx);
    scatter(diceClose, recallClose, 60, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.7);
    hold on;
    validMask = ~isnan(diceClose) & ~isnan(recallClose);
    if sum(validMask) > 2
        pCoeffs = polyfit(diceClose(validMask), recallClose(validMask), 1);
        xRange = linspace(min(diceClose), max(diceClose), 100);
        plot(xRange, polyval(pCoeffs, xRange), '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
    end
    hold off;
    xlabel('Dice Coefficient', 'FontSize', 13);
    ylabel('\Delta Recall', 'FontSize', 13);
    title(sprintf('CAB-NP: %s\nSpearman r = %.3f, p = %.3f', netName, rC, pC), 'FontSize', 12);
    grid on; set(gca, 'GridAlpha', 0.3, 'FontSize', 11);
    print(fScatter, '-dtiff', '-r300', fullfile(surfaceDir, [safeNetName, '_scatter.tiff']));
    close(fScatter);
end

%% Add Benchmark Result (Verbal Encoding Memory)
uniqueNets{end+1, 1} = 'Verbal encoding memory';
nPositiveRegions(end+1, 1) = NaN;
diceWithRef(end+1, 1) = NaN;
rCloseVec(end+1, 1) = 0.5848;
pCloseVec(end+1, 1) = 0.00005;
rRandomVec(end+1, 1) = NaN;
pRandomVec(end+1, 1) = NaN;
rAllVec(end+1, 1) = NaN;
pAllVec(end+1, 1) = NaN;

%% Save Results
resultTable = table(uniqueNets, nPositiveRegions, diceWithRef, ...
    rCloseVec, pCloseVec, rRandomVec, pRandomVec, rAllVec, pAllVec, ...
    'VariableNames', {'Network', 'N_Regions', 'Dice_with_Ref', ...
    'R_Close', 'P_Close', 'R_Random', 'P_Random', 'R_All', 'P_All'});

% FDR Correction for P_Close
validIdx = ~isnan(resultTable.P_Close);
if any(validIdx)
    pValues = resultTable.P_Close(validIdx);
    % Benjamini-Hochberg FDR
    qValues = mafdr(pValues, 'BHFDR', true);
    resultTable.Q_Close = nan(height(resultTable), 1);
    resultTable.Q_Close(validIdx) = qValues;
end

resultTable = sortrows(resultTable, 'R_Close', 'descend');
csvPath = fullfile(outputDir, 'CABNP_Network_Rankings.csv');
writetable(resultTable, csvPath);
fprintf('\nResults saved to: %s\n', csvPath);
disp(resultTable);

%% Bar Plot
f = figure('Visible', 'off');
set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 450]);

barValues = resultTable.R_Close;
pVals     = resultTable.P_Close;
qVals     = resultTable.Q_Close;
diceVals  = resultTable.Dice_with_Ref;
nTotal    = height(resultTable);
barNames  = cell(nTotal, 1);

for i = 1:nTotal
    if isnan(diceVals(i))
        barNames{i} = resultTable.Network{i};
    else
        barNames{i} = sprintf('%s (Dice=%.2f)', resultTable.Network{i}, diceVals(i));
    end
end

% Sky blue
bh = barh(1:nTotal, barValues, 'FaceColor', [0.53 0.81 0.98], 'EdgeColor', 'none');
set(gca, 'YTick', 1:nTotal, 'YTickLabel', barNames);
set(gca, 'YDir', 'reverse');
xlabel('Spearman R (Close Group)', 'FontSize', 14);
title('CAB-NP Network Congruence vs Recall', 'FontSize', 15);
set(gca, 'FontSize', 11);
xlim([min(0, min(barValues)-0.05), 0.7]);

% Significance markers (FDR-corrected)
hold on;
for i = 1:nTotal
    if ~isnan(qVals(i)) && qVals(i) < 0.05
        % ** = FDR Significant
        text(barValues(i) + (sign(barValues(i))*0.01), i, '**', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'VerticalAlignment', 'middle');
    elseif ~isnan(pVals(i)) && pVals(i) < 0.05
        % * = Uncorrected Significant
        text(barValues(i) + (sign(barValues(i))*0.01), i, '*', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'VerticalAlignment', 'middle');
    end
end
hold off;
grid on; set(gca, 'GridAlpha', 0.3);

barFile = fullfile(outputDir, 'CABNP_R_rankings.tiff');
print(gcf, '-dtiff', '-r300', barFile);
close(f);
fprintf('Bar plot saved: %s\n', barFile);

fprintf('\nDone.\n');
