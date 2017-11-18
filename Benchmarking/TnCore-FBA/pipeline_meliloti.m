%% Prepare models

initCobraToolbox;
clear;
prepareModels_meliloti;

%% iGD1575 core model, FBA, no RNA-seq

% Load the input
load('iGD1575.mat');
load('inputVariables_iGD1575.mat');

% Set growth threshold
solutionOriginal = optimizeCbModel(iGD1575,'max');
growthThresh = 0.5 * solutionOriginal.f;

% Run the analysis
[~, ~, ~, coreModelA, coreModelFastA, reducedModelA, ~, ~, ~] = ...
    tncore_main(iGD1575, 50000, growthThresh, tnseq, {1}, [], []);

% Save the data
clearvars -except *A
save('outputA.mat');
save('coreModelA.mat', 'coreModelA');
clear;

%% iGD1575 core model, FBA, with RNA-seq

% Load the input
load('iGD1575.mat');
load('inputVariables_iGD1575.mat');

% Set growth threshold
solutionOriginal = optimizeCbModel(iGD1575,'max');
growthThresh = 0.5 * solutionOriginal.f;

% Run the analysis
[~, ~, ~, coreModelB, coreModelFastB, reducedModelB, ~, ~, ~] = ...
    tncore_main(iGD1575, 50000, growthThresh, tnseq, {1}, rnaseq, []);

% Save the data
clearvars -except *B
save('outputB.mat');
save('coreModelB.mat', 'coreModelB');
clear;

%% Refinement of iGD1575

% Load data
load('iGD1575.mat');
load('inputVariables_iGD1575.mat');

% Set growth threshold
solutionOriginal = optimizeCbModel(iGD1575,'max');
growthThresh = 0.5 * solutionOriginal.f;

% Perform the refinement
[refinedModel1] = tncore_refine(iGD1575, tnseq, 25000, growthThresh, [], 0);

% Save the data
clearvars -except *1
save('output_refined_iGD1575.mat');
clear;

%% Refinement of draft S. meliloti model

% Load data
load('draftMeliloti.mat');
load('inputVariables_iGD1575.mat');

% Set growth threshold
solutionOriginal = optimizeCbModel(draftMeliloti,'max');
growthThresh = 0.5 * solutionOriginal.f;

% Perform the refinement
[refinedModel2] = tncore_refine(draftMeliloti, tnseq, 25000, growthThresh, [], 0);

% Save the data
clearvars -except *2
save('output_refined_draftMeliloti.mat');
clear;

quit
