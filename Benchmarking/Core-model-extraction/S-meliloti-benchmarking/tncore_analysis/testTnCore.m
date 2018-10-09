%% Set up workspave

clear all
addpath(genpath('../../software/cobratoolbox/'));
rmpath(genpath('../../software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../../software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../../software/tiger/'));
addpath(genpath('../../software/TnCore/'));
addpath(genpath('../../software/FastCore/'));
initCobraToolbox;
changeCobraSolver('gurobi', 'all');
rmpath(genpath('../../software/FastCore/'));
addpath(genpath('../../software/FastCore/'));

%% Load the data

% Get the model
load('iGD1575.mat');

% Get omics data
tnseq = table2cell(readtable('TnSeqData.txt', 'ReadVariableNames', false, ...
    'ReadRowNames', false, 'Delimiter', '\t'));
rnaseq = table2cell(readtable('RnaSeqData.txt', 'ReadVariableNames', false, ...
    'ReadRowNames', false, 'Delimiter', '\t'));

%% Make models

% Without RNA-seq data
[coreModelA, reducedModelA] = tncore_core(iGD1575, tnseq, [], ...
    [], [], [], [], [], [], 1);

% With RNA-seq data without keeping high
[coreModelB, reducedModelB] = tncore_core(iGD1575, tnseq, [], ...
    [], [], [], rnaseq, [], 0, 1);

% With RNA-seq data with keeping high
[coreModelC, reducedModelC] = tncore_core(iGD1575, tnseq, [], ...
    [], [], [], rnaseq, [], 1, 1);

%% Save and clear

save('allWorkspace.mat');
save('coreModelA.mat', 'coreModelA');
save('coreModelB.mat', 'coreModelB');
save('coreModelC.mat', 'coreModelC');

clear;

