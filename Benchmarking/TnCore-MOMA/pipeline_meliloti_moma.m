%% Prepare matlab

initCobraToolbox;
clear;

%% iGD1575 core model, MOMA, without RNA-seq

% Load the input
load('iGD1575.mat');
load('inputVariables_iGD1575.mat');

% Set growth threshold
solutionOriginal = optimizeCbModel(iGD1575,'max');
growthThresh = 0.5 * solutionOriginal.f;

% Run the analysis
[~, ~, ~, coreModelC, coreModelFastC, reducedModelC, ~, ~, ~] = ...
    tncore_main(iGD1575, 50000, growthThresh, tnseq, {1}, [], [], 2);

% Save the data
clearvars -except *C
save('outputC.mat');
save('coreModelC.mat', 'coreModelC');
clear;

%% iGD1575 core model, MOMA, with RNA-seq

% Load the input
load('iGD1575.mat');
load('inputVariables_iGD1575.mat');

% Set growth threshold
solutionOriginal = optimizeCbModel(iGD1575,'max');
growthThresh = 0.5 * solutionOriginal.f;

% Run the analysis
[~, ~, ~, coreModelD, coreModelFastD, reducedModelD, ~, ~, ~] = ...
    tncore_main(iGD1575, 50000, growthThresh, tnseq, {1}, rnaseq, [], 2);

% Save the data
clearvars -except *D
save('outputD.mat');
save('coreModelD.mat', 'coreModelD');
clear;

%% Exit matlab

quit
