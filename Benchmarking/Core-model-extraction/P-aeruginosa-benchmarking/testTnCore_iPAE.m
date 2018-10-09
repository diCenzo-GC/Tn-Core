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
load('iPAE1146.mat');

% Get omics data
load('inputVariables_iPAE1146.mat');

%% Make Tn-Core models

% Tn-Core Without RNA-seq data
[coreModel_tncore_1, reducedModel_1] = tncore_core(iPAE1146, tnseq, [], ...
    [], [], [], [], [], [], 1);

% Tn-Core with RNA-seq data
[coreModel_tncore_2, reducedModel_2] = tncore_core(iPAE1146, tnseq, [], ...
    [], [], [], rnaseq, [], [], 1);

%% Prepare model for other methods

% Refine the model to remove dead ends
model = iPAE1146;
deadEndMetabolites = detectDeadEnds(model);
if ~isempty(deadEndMetabolites)
    test = 1;
    while test > 0
        deadEndMetabolites = detectDeadEnds(model);
        if isempty(deadEndMetabolites)
            test = 0;
        else
            deadEndMetabolites3 = cell(length(deadEndMetabolites),1);
            for n = 1:length(deadEndMetabolites);
                deadEndMetabolites2 = num2cell(deadEndMetabolites);
                deadEndMetabolites3{n,1} = ...
                    model.mets{deadEndMetabolites2{n,1},1};
            end
            [deadEndReactions] = findRxnsFromMets(model,deadEndMetabolites3);
            model = removeRxns(model,deadEndReactions);
            if length(model.mets) < size(model.S, 1)
                model.S(deadEndMetabolites,:) = [];
                model.b(deadEndMetabolites,:) = [];
            end
            clear deadEndMetabolites
            clear deadEndMetabolites2
        end
    end
    model = tncore_remove(model);
end

% Remove unnecessary unknowns from the model
unknownID = findGeneIDs(model,'Unknown');
A = [' | x(' num2str(unknownID) ')'];
B = ['x(' num2str(unknownID) ') | '];
model.grRules = strrep(model.grRules,' or Unknown','');
model.grRules = strrep(model.grRules,'Unknown or ','');
model.rxnNotes = strrep(model.rxnNotes,' or Unknown','');
model.rxnNotes = strrep(model.rxnNotes,'Unknown or ','');
model.rules = strrep(model.rules,A,'');
model.rules = strrep(model.rules,B,'');
rxnGeneMatNewFull = cell(length(model.rxns),length(model.genes));
for n = 1:length(model.genes)
    [~,ListResults] = findRxnsFromGenes(model, model.genes{n},0,1);
    rxnIDs = num2cell(findRxnIDs(model,ListResults(:,1)));
    for m = 1:length(rxnIDs)
        rxnGeneMatNewFull{rxnIDs{m},n} = 1;
    end
end
temp = cellfun('isempty',rxnGeneMatNewFull);
rxnGeneMatNewFull(temp) = {0};
rxnGeneMatNewFull = cell2mat(rxnGeneMatNewFull);
rxnGeneMatNew = sparse(double(rxnGeneMatNewFull));
model.rxnGeneMat = rxnGeneMatNew;
iPAE1146 = model;

% Save updated model
save('iPAE1146_modified.mat', 'iPAE1146');

%% Prepare Tn-seq data for GIMME

% Convert tnseq data to log scale
for n = 1:length(tnseq)
    if strmatch(0,tnseq{n,1},'exact')
        tnseq{n,1} = [];
    end
end
is0 = cellfun('isempty',tnseq);
tnseq(is0) = num2cell(min(cell2mat(tnseq(:,1))));
tnseq(:,1) = num2cell(log10(cell2mat(tnseq(:,1))));

% Calculate statistics of tnseq data
allTnseq = num2cell(sort(cell2mat(tnseq(:,1))));
totalTnseq = length(allTnseq);
prc25 = round(0.25 * totalTnseq);
prc75 = round(0.75 * totalTnseq);
IQrange = allTnseq{prc75} - allTnseq{prc25};
outThresh = allTnseq{prc25} - IQrange;
isOutlier = cell2mat(allTnseq) < outThresh;
temp = allTnseq(isOutlier);
nonOutliers = setdiff(cell2mat(allTnseq),cell2mat(temp));
medianValue = median(nonOutliers);
stdev = std(nonOutliers);
expressThresh = medianValue - (3.5 * stdev);

% Pull out just model genes and data
modelTnseq = cell(length(iPAE1146.genes),2);
for n = 1:length(iPAE1146.genes)
    ID = strmatch(iPAE1146.genes{n},tnseq(:,2));
    if isempty(ID)
        modelTnseq{n,2} = iPAE1146.genes{n};
        modelTnseq{n,1} = max(cell2mat(tnseq(:,1)));
    else
        modelTnseq{n,2} = iPAE1146.genes{n};
        modelTnseq{n,1} = tnseq{ID,1};
    end
end

% Convert tnseq data format
iPAE1146_express = modelTnseq(:,1);
iPAE1146_genes = modelTnseq(:,2);
iPAE1146_express = num2cell(cell2mat(iPAE1146_express) * -1);
minValue = min(cell2mat(iPAE1146_express));
for n = 1:length(iPAE1146_express)
    iPAE1146_express{n} = iPAE1146_express{n} - minValue;
end
expressThresh = (expressThresh * -1) - minValue;

%% GIMME with Tn-seq data

% Initiate tiger
start_tiger('gurobi');
iPAE1146_tiger = cobra_to_tiger(iPAE1146);

% Prepare expression data
iPAE1146_express = cell2mat(iPAE1146_express);

% Run GIMME
[GENE_STATES_iPAE1146, GENES_iPAE1146, SOL_iPAE1146, TIGER_iPAE1146] = ...
    tncore_gimme(iPAE1146_tiger, iPAE1146_express, expressThresh, ...
    iPAE1146_genes, 0.5);

% Rebuild the core model in COBRA format

isPresent = logical(GENE_STATES_iPAE1146);
genesToDelete = GENES_iPAE1146(~isPresent);
[coreModel,~,constrRxnNames] = deleteModelGenes(iPAE1146, genesToDelete);
testSol = optimizeCbModel(coreModel);
origSol = optimizeCbModel(iPAE1146);
if testSol.f < 0.5 * origSol.f
    for n = 1:length(genesToDelete)
        geneID = findGeneIDs(iPAE1146, genesToDelete{n,1});
        genesToDelete{n,2} = iPAE1146_express(geneID);
    end
end
genesToDelete = sortrows(genesToDelete,2);
tempModel = iPAE1146;
allRxnsToRemove = {};
allNewGenesToDelete = {};
for n = 1:length(genesToDelete)
    [~,~,constrRxnNames] = deleteModelGenes(tempModel, genesToDelete{n,1});
    testModel = removeRxns(tempModel, constrRxnNames);
    testSol = optimizeCbModel(testModel);
    if testSol.f >= 0.5 * origSol.f
        [tempModel,~,constrRxnNames] = deleteModelGenes(tempModel, genesToDelete{n,1});
        allRxnsToRemove = vertcat(allRxnsToRemove, constrRxnNames);
        allNewGenesToDelete = vertcat(allNewGenesToDelete, genesToDelete{n,1});
    end
end
[tempModel,~,constrRxnNames] = deleteModelGenes(iPAE1146, allNewGenesToDelete);
tempModel = removeRxns(tempModel, constrRxnNames);
coreModel_gimme = tncore_delete(tempModel);

clearvars -except coreModel_gimme coreModel_tncore_1 coreModel_tncore_2

%% Identify core reactions

% Load inputs
load('iPAE1146_modified.mat');
load('inputVariables_iPAE1146.mat');

% Convert tnseq data to log scale
for n = 1:length(tnseq)
    if strmatch(0,tnseq{n,1},'exact')
        tnseq{n,1} = [];
    end
end
is0 = cellfun('isempty',tnseq);
tnseq(is0) = num2cell(min(cell2mat(tnseq(:,1))));
tnseq(:,1) = num2cell(log10(cell2mat(tnseq(:,1))));

% Pull out tnseq data for model genes
modelTnseq = cell(length(iPAE1146.genes),2);
for n = 1:length(iPAE1146.genes)
    ID = strmatch(iPAE1146.genes{n},tnseq(:,2));
    if isempty(ID)
        modelTnseq{n,2} = iPAE1146.genes{n};
        modelTnseq{n,1} = max(cell2mat(tnseq(:,1)));
    else
        modelTnseq{n,2} = iPAE1146.genes{n};
        modelTnseq{n,1} = tnseq{ID,1};
    end
end

% Calculate statistics of tnseq data
allTnseq = num2cell(sort(cell2mat(tnseq(:,1))));
totalTnseq = length(allTnseq);
prc25 = round(0.25 * totalTnseq);
prc75 = round(0.75 * totalTnseq);
IQrange = allTnseq{prc75} - allTnseq{prc25};
outThresh = allTnseq{prc25} - IQrange;
isOutlier = cell2mat(allTnseq) < outThresh;
temp = allTnseq(isOutlier);
nonOutliers = setdiff(cell2mat(allTnseq),cell2mat(temp));
medianValue = median(nonOutliers);
stdev = std(nonOutliers);

% Find essential genes
modelTnseqGimme = modelTnseq;
essentialGenes = {};
for n = 1:length(modelTnseqGimme)
    if modelTnseqGimme{n,1} < medianValue - (3.5 * stdev)
        essentialGenes = vertcat(essentialGenes, modelTnseqGimme{n,2});
    end
end

% Get core reactions
[~, ~, coreRxns] = deleteModelGenes(iPAE1146, essentialGenes);

clearvars -except coreModel_gimme coreModel_tncore_1 coreRxns iPAE1146 ...
    essentialGenes coreModel_tncore_2

%% Make fastcore model

% Run fastcc
A = fastcc(iPAE1146, 1.1e-6);

% Make the consistent model
A = iPAE1146.rxns(A);
B = setdiff(iPAE1146.rxns, A);
iPAE1146_consistent = removeRxns(iPAE1146, B);

% Identify core reactions
coreRxns2 = intersect(coreRxns, iPAE1146_consistent.rxns);
coreReactions = findRxnIDs(iPAE1146_consistent, coreRxns2);

% Run fastcore
finalRxns = fastcore(coreReactions, iPAE1146_consistent, 1.1e-6);

% Prepare core model
finalRxns = iPAE1146_consistent.rxns(finalRxns);
rxnsToRemove = setdiff(iPAE1146_consistent.rxns, finalRxns);
coreModel_fastcore = removeRxns(iPAE1146_consistent, rxnsToRemove);
coreModel_fastcore = tncore_remove(coreModel_fastcore);

clearvars -except coreModel_gimme coreModel_tncore_1 coreModel_fastcore ...
    coreRxns iPAE1146 iPAE1146_consistent coreReactions coreModel_tncore_2 ...
    essentialGenes

save('temp.mat');
save('iPAE1146_consistent.mat', 'iPAE1146_consistent');
clear;

%% Make minnw model

% Run minnw
coreModel_minnw_intermediate = FaMoRe_iPAE1146();

% Load previous stuff
load('temp.mat');

% Produce core model with genes included
rxnsToRemove = setdiff(iPAE1146_consistent.rxns, coreModel_minnw_intermediate.rxns);
coreModel_minnw = removeRxns(iPAE1146_consistent, rxnsToRemove);
coreModel_minnw = tncore_remove(coreModel_minnw);

clearvars -except coreModel_gimme coreModel_tncore_1 coreModel_fastcore ...
    coreRxns iPAE1146 iPAE1146_consistent coreReactions coreModel_minnw ...
    essentialGenes coreModel_tncore_2

%% Make GIMME model with RNA-seq

% Load inputs
load('inputVariables_iPAE1146.mat');

% Pull out just model genes and data
genes = iPAE1146.genes;
modelExpress = cell(length(genes),1);
for n = 1:length(genes)
    for m = 1:length(rnaseq)
        if strmatch(genes{n},rnaseq{m,2},'exact')
            modelExpress{n} = rnaseq{m};
        end
    end
    if isempty(modelExpress{n})
        modelExpress{n} = 0;
    end
end
iPAE1146_express = modelExpress;
iPAE1146_genes = genes;

% Determine expression threshold
expressThresh = sum(cell2mat(rnaseq(:,1))) / 5000;

% Initiate tiger
start_tiger('gurobi');
iPAE1146_tiger = cobra_to_tiger(iPAE1146);

% Prepare expression data
iPAE1146_express = cell2mat(iPAE1146_express);

% Run GIMME
[GENE_STATES_iPAE1146, GENES_iPAE1146, SOL_iPAE1146, TIGER_iPAE1146] = ...
    tncore_gimme(iPAE1146_tiger, iPAE1146_express, expressThresh, ...
    iPAE1146_genes, 0.5);

% Rebuild the core model in COBRA format
isPresent = logical(GENE_STATES_iPAE1146);
genesToDelete = GENES_iPAE1146(~isPresent);
[coreModel,~,constrRxnNames] = deleteModelGenes(iPAE1146, genesToDelete);
testSol = optimizeCbModel(coreModel);
origSol = optimizeCbModel(iPAE1146);
if testSol.f < 0.5 * origSol.f
    for n = 1:length(genesToDelete)
        geneID = findGeneIDs(iPAE1146, genesToDelete{n,1});
        genesToDelete{n,2} = iPAE1146_express(geneID);
    end
end
genesToDelete = sortrows(genesToDelete,2);
tempModel = iPAE1146;
allRxnsToRemove = {};
allNewGenesToDelete = {};
for n = 1:length(genesToDelete)
    [~,~,constrRxnNames] = deleteModelGenes(tempModel, genesToDelete{n,1});
    testModel = removeRxns(tempModel, constrRxnNames);
    testSol = optimizeCbModel(testModel);
    if testSol.f >= 0.5 * origSol.f
        [tempModel,~,constrRxnNames] = deleteModelGenes(tempModel, genesToDelete{n,1});
        allRxnsToRemove = vertcat(allRxnsToRemove, constrRxnNames);
        allNewGenesToDelete = vertcat(allNewGenesToDelete, genesToDelete{n,1});
    end
end
[tempModel,~,constrRxnNames] = deleteModelGenes(iPAE1146, allNewGenesToDelete);
tempModel = removeRxns(tempModel, constrRxnNames);
coreModel_gimme_rna = tncore_delete(tempModel);

clearvars -except coreModel_gimme coreModel_tncore_1 coreModel_fastcore ...
    coreRxns iPAE1146 iPAE1146_consistent coreReactions coreModel_minnw ...
    coreModel_gimme_rna coreModel_tncore_2 essentialGenes

%% Compare gene essentialities

% Single gene deletion analysis
iPAE1146_grRatio = singleGeneDeletion(iPAE1146, 'MOMA');
coreModel_tncore_1_grRatio = singleGeneDeletion(coreModel_tncore_1, 'MOMA');
coreModel_tncore_2_grRatio = singleGeneDeletion(coreModel_tncore_2, 'MOMA');
coreModel_fastcore_grRatio = singleGeneDeletion(coreModel_fastcore, 'MOMA');
coreModel_minnw_grRatio = singleGeneDeletion(coreModel_minnw, 'MOMA');
coreModel_gimme_grRatio = singleGeneDeletion(coreModel_gimme, 'MOMA');
coreModel_gimme_rna_grRatio = singleGeneDeletion(coreModel_gimme_rna, 'MOMA');

% List of genes of interest
carbonGenes = {'PA3183'; 'PA3182'; 'PA3181'; 'PA0330'; ...
    'PA0548'; 'PA2796'; 'PA3001'; 'PA0552'; 'PA3635'; 'PA4329'; ...
    'PA3416'; 'PA1580'; 'PA1562'; 'PA2623'; 'PA1585'; 'PA1583'; ...
    'PA4470'; 'PA3452'; 'PA5555'; 'PA5192'; 'PA3471'};

% Make output variable
outputCompared = cell(length(carbonGenes), 8);

% Pull out gene numbers
iPAE1146_IDs = findGeneIDs(iPAE1146, carbonGenes);
coreModel_tncore_1_IDs = findGeneIDs(coreModel_tncore_1, carbonGenes);
coreModel_tncore_2_IDs = findGeneIDs(coreModel_tncore_2, carbonGenes);
coreModel_fastcore_IDs = findGeneIDs(coreModel_fastcore, carbonGenes);
coreModel_minnw_IDs = findGeneIDs(coreModel_minnw, carbonGenes);
coreModel_gimme_IDs = findGeneIDs(coreModel_gimme, carbonGenes);
coreModel_gimme_rna_IDs = findGeneIDs(coreModel_gimme_rna, carbonGenes);

% Get the data
for n = 1:length(carbonGenes)
    
    if iPAE1146_IDs(n) ~= 0
        outputCompared{n,1} = round(iPAE1146_grRatio(iPAE1146_IDs(n)),3);
    else
        outputCompared{n,1} = -1;
    end
    
    if coreModel_tncore_1_IDs(n) ~= 0
        outputCompared{n,2} = round(coreModel_tncore_1_grRatio(coreModel_tncore_1_IDs(n)),3);
    else
        outputCompared{n,2} = -1;
    end

     if coreModel_tncore_2_IDs(n) ~= 0
         outputCompared{n,3} = round(coreModel_tncore_2_grRatio(coreModel_tncore_2_IDs(n)),3);
     else
         outputCompared{n,3} = -1;
     end

    if coreModel_fastcore_IDs(n) ~= 0
        outputCompared{n,4} = round(coreModel_fastcore_grRatio(coreModel_fastcore_IDs(n)),3);
    else
        outputCompared{n,4} = -1;
    end

    if coreModel_minnw_IDs(n) ~= 0
        outputCompared{n,5} = round(coreModel_minnw_grRatio(coreModel_minnw_IDs(n)),3);
    else
        outputCompared{n,5} = -1;
    end
    
     if coreModel_gimme_IDs(n) ~= 0
         outputCompared{n,6} = round(coreModel_gimme_grRatio(coreModel_gimme_IDs(n)),3);
     else
         outputCompared{n,6} = -1;
     end
     
     if coreModel_gimme_rna_IDs(n) ~= 0
         outputCompared{n,7} = round(coreModel_gimme_rna_grRatio(coreModel_gimme_rna_IDs(n)),3);
     else
         outputCompared{n,7} = -1;
     end
    
end

% Record Tn-seq data
load('inputVariables_iPAE1146.mat');
for n = 1:length(carbonGenes)
    pos = strmatch(carbonGenes(n), tnseq(:,2), 'exact');
    outputCompared{n,8} = tnseq{pos, 1};
end

% Add labels
headers = {'Gene', 'iPAE1146', 'tncore_1', 'tncore_2', 'fastcore', ...
    'minnw', 'gimme', 'gimme_rna', 'TnSeq'};
outputCompared = horzcat(carbonGenes, outputCompared);
outputCompared = vertcat(headers, outputCompared);
outputCompared = transpose(outputCompared);

%% Save and clear

% Save model and workspace
save('allWorkspace_iPAE.mat');

% Export comparison as an Excell table
toExport = cell2table(outputCompared);
writetable(toExport, 'iPAE1146_carbonGenesCompared.xlsx', 'WriteVariableNames', false);

clear;

