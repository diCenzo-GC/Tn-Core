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

% Get the iRsp1140 model
iRsp1140 = readCbModel('iRsp1140.xml');

% Get omics data
load('inputVariables_iRsp1140.mat');

% Set objective
iRsp1140 = changeObjective(iRsp1140, 'RXN1391');

% Set medium
EX_list = iRsp1140.rxns(findExcRxns(iRsp1140));
iRsp1140.rev(findExcRxns(iRsp1140)) = 1;
iRsp_medium = {'RXN0217';'RXN0223';'RXN0213';'RXN0222';'RXN0219';'RXN0191';...
    'RXN0190';'RXN0224';'RXN0196';'RXN0188';'RXN0221';'RXN1395';'RXN1326';...
    'RXN0193';'RXN1167';'RXN1329';'RXN1158';'RXN1354';'RXN0220';'RXN0192'};
iRsp1140 = changeRxnBounds(iRsp1140, EX_list, 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, EX_list, 1000, 'u');
iRsp1140 = changeRxnBounds(iRsp1140, iRsp_medium, -100, 'l');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0224', -5, 'l');

% Set additional model constraints
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0109', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0631', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN1827', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN1843', 0, 'b');

% Save the model
save('iRsp1140.mat', 'iRsp1140');

%% Make Tn-Core model

% Tn-Core with custom bounds
[coreModel_tncore_1, reducedModel_iRsp] = tncore_core(iRsp1140, tnseq, [], ...
    [], [], [], [], [], [], 1);

%% Prepare model for other methods

% Refine the model to remove dead ends
model = iRsp1140;
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
iRsp1140 = model;

% Save updated model
save('iRsp1140_modified.mat', 'iRsp1140');

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
modelTnseq = cell(length(iRsp1140.genes),2);
for n = 1:length(iRsp1140.genes)
    ID = strmatch(iRsp1140.genes{n},tnseq(:,2));
    if isempty(ID)
        modelTnseq{n,2} = iRsp1140.genes{n};
        modelTnseq{n,1} = max(cell2mat(tnseq(:,1)));
    else
        modelTnseq{n,2} = iRsp1140.genes{n};
        modelTnseq{n,1} = tnseq{ID,1};
    end
end

% Convert tnseq data format
iRsp1140_express = modelTnseq(:,1);
iRsp1140_genes = modelTnseq(:,2);
iRsp1140_express = num2cell(cell2mat(iRsp1140_express) * -1);
minValue = min(cell2mat(iRsp1140_express));
for n = 1:length(iRsp1140_express)
    iRsp1140_express{n} = iRsp1140_express{n} - minValue;
end
expressThresh = (expressThresh * -1) - minValue;

%% GIMME with Tn-seq data

% Initiate tiger
start_tiger('gurobi');
iRsp1140_tiger = cobra_to_tiger(iRsp1140);

% Prepare expression data
iRsp1140_express = cell2mat(iRsp1140_express);

% Run GIMME
[GENE_STATES_iRsp1140, GENES_iRsp1140, SOL_iRsp1140, TIGER_iRsp1140] = ...
    tncore_gimme(iRsp1140_tiger, iRsp1140_express, expressThresh, ...
    iRsp1140_genes, 0.5);

% Rebuild the core model in COBRA format

isPresent = logical(GENE_STATES_iRsp1140);
genesToDelete = GENES_iRsp1140(~isPresent);
[coreModel,~,constrRxnNames] = deleteModelGenes(iRsp1140, genesToDelete);
testSol = optimizeCbModel(coreModel);
origSol = optimizeCbModel(iRsp1140);
if testSol.f < 0.5 * origSol.f
    for n = 1:length(genesToDelete)
        geneID = findGeneIDs(iRsp1140, genesToDelete{n,1});
        genesToDelete{n,2} = iRsp1140_express(geneID);
    end
end
genesToDelete = sortrows(genesToDelete,2);
tempModel = iRsp1140;
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
[tempModel,~,constrRxnNames] = deleteModelGenes(iRsp1140, allNewGenesToDelete);
tempModel = removeRxns(tempModel, constrRxnNames);
coreModel_gimme = tncore_delete(tempModel);

clearvars -except coreModel_gimme coreModel_tncore_1

%% Identify core reactions

% Load inputs
load('iRsp1140_modified.mat');
load('inputVariables_iRsp1140.mat');

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
modelTnseq = cell(length(iRsp1140.genes),2);
for n = 1:length(iRsp1140.genes)
    ID = strmatch(iRsp1140.genes{n},tnseq(:,2));
    if isempty(ID)
        modelTnseq{n,2} = iRsp1140.genes{n};
        modelTnseq{n,1} = max(cell2mat(tnseq(:,1)));
    else
        modelTnseq{n,2} = iRsp1140.genes{n};
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
[~, ~, coreRxns] = deleteModelGenes(iRsp1140, essentialGenes);

clearvars -except coreModel_gimme coreModel_tncore_1 ...
    coreRxns iRsp1140 essentialGenes

%% Make fastcore model

% Run fastcc
A = fastcc(iRsp1140, 1.1e-6);

% Make the consistent model
A = iRsp1140.rxns(A);
B = setdiff(iRsp1140.rxns, A);
iRsp1140_consistent = removeRxns(iRsp1140, B);

% Identify core reactions
coreRxns2 = intersect(coreRxns, iRsp1140_consistent.rxns);
coreReactions = findRxnIDs(iRsp1140_consistent, coreRxns2);

% Run fastcore
finalRxns = fastcore(coreReactions, iRsp1140_consistent, 1.1e-6);

% Prepare core model
finalRxns = iRsp1140_consistent.rxns(finalRxns);
rxnsToRemove = setdiff(iRsp1140_consistent.rxns, finalRxns);
coreModel_fastcore = removeRxns(iRsp1140_consistent, rxnsToRemove);
coreModel_fastcore = tncore_remove(coreModel_fastcore);

clearvars -except coreModel_gimme coreModel_tncore_1 coreModel_fastcore ...
    coreRxns iRsp1140 iRsp1140_consistent coreReactions essentialGenes

save('temp.mat');
save('iRsp1140_consistent.mat', 'iRsp1140_consistent');
clear;

%% Make minnw model

% Run minnw
coreModel_minnw_intermediate = FaMoRe_iRsp1140();

% Load previous stuff
load('temp.mat');

% Produce core model with genes included
rxnsToRemove = setdiff(iRsp1140_consistent.rxns, coreModel_minnw_intermediate.rxns);
coreModel_minnw = removeRxns(iRsp1140_consistent, rxnsToRemove);
coreModel_minnw = tncore_remove(coreModel_minnw);

clearvars -except coreModel_gimme coreModel_tncore_1 ...
    coreModel_fastcore coreRxns iRsp1140 iRsp1140_consistent coreReactions ...
    coreModel_minnw essentialGenes 

%% Compare gene essentialities

% Single gene deletion analysis
iRsp1140_grRatio = singleGeneDeletion(iRsp1140, 'MOMA');
coreModel_tncore_grRatio = singleGeneDeletion(coreModel_tncore_1, 'MOMA');
coreModel_fastcore_grRatio = singleGeneDeletion(coreModel_fastcore, 'MOMA');
coreModel_minnw_grRatio = singleGeneDeletion(coreModel_minnw, 'MOMA');
coreModel_gimme_grRatio = singleGeneDeletion(coreModel_gimme, 'MOMA');
iRsp1140_grRatio_FBA = singleGeneDeletion(iRsp1140, 'FBA');
coreModel_tncore_grRatio_FBA = singleGeneDeletion(coreModel_tncore_1, 'FBA');
coreModel_fastcore_grRatio_FBA = singleGeneDeletion(coreModel_fastcore, 'FBA');
coreModel_minnw_grRatio_FBA = singleGeneDeletion(coreModel_minnw, 'FBA');
coreModel_gimme_grRatio_FBA = singleGeneDeletion(coreModel_gimme, 'FBA');

% List of genes of interest
carbonGenes = {'RSP_2734'; 'RSP_2735'; 'RSP_2645'; 'RSP_0359'; ...
    'RSP_2956'; 'RSP_0561'; 'RSP_2959'; 'RSP_4044'; 'RSP_2491'; 'RSP_1848'; ...
    'RSP_4047'; 'RSP_1994'; 'RSP_1806'; 'RSP_1559'; 'RSP_0965'; 'RSP_0976'; ...
    'RSP_2138'; 'RSP_0968'; 'RSP_2090'; 'RSP_1680'; 'RSP_1593'; 'RSP_2298'};

% Make output variable
outputCompared_MOMA = cell(length(carbonGenes), 6);

% Pull out gene numbers
iRsp1140_IDs = findGeneIDs(iRsp1140, carbonGenes);
coreModel_tncore_1_IDs = findGeneIDs(coreModel_tncore_1, carbonGenes);
coreModel_fastcore_IDs = findGeneIDs(coreModel_fastcore, carbonGenes);
coreModel_minnw_IDs = findGeneIDs(coreModel_minnw, carbonGenes);
coreModel_gimme_IDs = findGeneIDs(coreModel_gimme, carbonGenes);

% Get the data
for n = 1:length(carbonGenes)
    
    if iRsp1140_IDs(n) ~= 0
        outputCompared_MOMA{n,1} = round(iRsp1140_grRatio(iRsp1140_IDs(n)),3);
    else
        outputCompared_MOMA{n,1} = -1;
    end
    
    if coreModel_tncore_1_IDs(n) ~= 0
        outputCompared_MOMA{n,2} = round(coreModel_tncore_grRatio(coreModel_tncore_1_IDs(n)),3);
    else
        outputCompared_MOMA{n,2} = -1;
    end

    if coreModel_fastcore_IDs(n) ~= 0
        outputCompared_MOMA{n,3} = round(coreModel_fastcore_grRatio(coreModel_fastcore_IDs(n)),3);
    else
        outputCompared_MOMA{n,3} = -1;
    end

    if coreModel_minnw_IDs(n) ~= 0
        outputCompared_MOMA{n,4} = round(coreModel_minnw_grRatio(coreModel_minnw_IDs(n)),3);
    else
        outputCompared_MOMA{n,4} = -1;
    end
    
     if coreModel_gimme_IDs(n) ~= 0
         outputCompared_MOMA{n,5} = round(coreModel_gimme_grRatio(coreModel_gimme_IDs(n)),3);
     else
         outputCompared_MOMA{n,5} = -1;
     end
     
end

% Record Tn-seq data
load('inputVariables_iRsp1140.mat');
for n = 1:length(carbonGenes)
    pos = strmatch(carbonGenes(n), tnseq(:,2), 'exact');
    if isempty(pos)
        outputCompared_MOMA{n,6} = -1;
    else
        outputCompared_MOMA{n,6} = tnseq{pos, 1};
    end
end

% Add labels
headers = {'Gene', 'iRsp1140', 'tncore_1', 'fastcore', 'minnw', 'gimme', 'TnSeq'};
outputCompared_MOMA = horzcat(carbonGenes, outputCompared_MOMA);
outputCompared_MOMA = vertcat(headers, outputCompared_MOMA);
outputCompared_MOMA = transpose(outputCompared_MOMA);

% Make output variable
outputCompared_FBA = cell(length(carbonGenes), 6);

% Pull out gene numbers
iRsp1140_IDs = findGeneIDs(iRsp1140, carbonGenes);
coreModel_tncore_1_IDs = findGeneIDs(coreModel_tncore_1, carbonGenes);
coreModel_fastcore_IDs = findGeneIDs(coreModel_fastcore, carbonGenes);
coreModel_minnw_IDs = findGeneIDs(coreModel_minnw, carbonGenes);
coreModel_gimme_IDs = findGeneIDs(coreModel_gimme, carbonGenes);

% Get the data
for n = 1:length(carbonGenes)
    
    if iRsp1140_IDs(n) ~= 0
        outputCompared_FBA{n,1} = round(iRsp1140_grRatio_FBA(iRsp1140_IDs(n)),3);
    else
        outputCompared_FBA{n,1} = -1;
    end
    
    if coreModel_tncore_1_IDs(n) ~= 0
        outputCompared_FBA{n,2} = round(coreModel_tncore_grRatio_FBA(coreModel_tncore_1_IDs(n)),3);
    else
        outputCompared_FBA{n,2} = -1;
    end

    if coreModel_fastcore_IDs(n) ~= 0
        outputCompared_FBA{n,3} = round(coreModel_fastcore_grRatio_FBA(coreModel_fastcore_IDs(n)),3);
    else
        outputCompared_FBA{n,3} = -1;
    end

    if coreModel_minnw_IDs(n) ~= 0
        outputCompared_FBA{n,4} = round(coreModel_minnw_grRatio_FBA(coreModel_minnw_IDs(n)),3);
    else
        outputCompared_FBA{n,4} = -1;
    end
    
     if coreModel_gimme_IDs(n) ~= 0
         outputCompared_FBA{n,5} = round(coreModel_gimme_grRatio_FBA(coreModel_gimme_IDs(n)),3);
     else
         outputCompared_FBA{n,5} = -1;
     end
     
end

% Record Tn-seq data
load('inputVariables_iRsp1140.mat');
for n = 1:length(carbonGenes)
    pos = strmatch(carbonGenes(n), tnseq(:,2), 'exact');
    if isempty(pos)
        outputCompared_FBA{n,6} = -1;
    else
        outputCompared_FBA{n,6} = tnseq{pos, 1};
    end
end

% Add labels
headers = {'Gene', 'iRsp1140', 'tncore_1', 'fastcore', 'minnw', 'gimme', 'TnSeq'};
outputCompared_FBA = horzcat(carbonGenes, outputCompared_FBA);
outputCompared_FBA = vertcat(headers, outputCompared_FBA);
outputCompared_FBA = transpose(outputCompared_FBA);

%% Save and clear

% Save model and workspace
save('allWorkspace_iRsp.mat');

% Export comparison as an Excell table
toExport = cell2table(outputCompared_MOMA);
writetable(toExport, 'iRsp1140_carbonGenesCompared_MOMA.xlsx', 'WriteVariableNames', false);
toExport = cell2table(outputCompared_FBA);
writetable(toExport, 'iRsp1140_carbonGenesCompared_FBA.xlsx', 'WriteVariableNames', false);

clear;

