function [contextModel, reducedModel] = tncore_core(model, tnseq, epsilon, ...
    binThresh, growthFrac, rnaseq, expressThresh, keepHigh, deadends)

%
% Uses FASTCORE and the gene-centric GIMME implementation of the TIGER 
% Toolbox to extract a context-specific core metabolic model from a 
% genome-scale metabolic model based on the Tn-seq data and optionally with
% RNA-seq data as well.
%
% USAGE
%   [contextModel, reducedModel] = tncore_core(model, tnseq, epsilon, ...
%       binThresh, growthFrac, rnaseq, expressThresh, keepHigh, deadends)
%
% INPUTS
%   model           COBRA model structure including gene-reaction
%                   associations
%   tnseq           Tn-seq data for all genes in the organism. First 
%                   column contains the Tn-seq scores, the second column 
%                   contains the gene names (Default = {})
%
% OPTIONAL INPUTS
%   epsilon         The epsilon value (the minimum flux rate that is
%                   considered non-zero) to be used during the running of
%                   FASTCC and FASTCORE. (Default = 1.1e-6)
%   binThresh       An array containing three values (number of standard
%                   deviations away from the mean of the log of the Tn-seq 
%                   data) to use for setting the limits of each bin. 
%                   (Default = {'3.5'; '2.5'; '1.5'}) 
%   growthFrac      The minimum allowable objective flux in the generated
%                   core model, as a fraction of growth of the input model
%                   (Default = 0.5) 
%   rnaseq          RNA-seq data for all genes in the organism. First
%                   column contains the RNA-seq data, the second column
%                   contains the gene names (Default = {})
%   expressThresh   The threshold for a gene to be considered expressed
%                   (Default = 0.02% of the sum of all expression values)
%   keepHigh        Indicates if highly expressed genes should be added
%                   back into the context-specific model, after the
%                   context-specific model has otherwise been prepared. Use 
%                   1 for yes, and 0 for no. (Default = 0)
%   deadends        Indicates if reactions producing dead-end metabolites
%                   should be removed from the core model. Use 1 for yes,
%                   and 0 for no. (Default = 1)
%
% OUTPUTS
%   contextModel    The context-specific core model based on the Tn-seq
%                   data (considering also the RNA-seq data, if provided)
%   reducedModel    The input model with dead-ends removed and with
%                   unnecessary 'Unknown' GPRs removed
%
% AUTHORS
%   George diCenzo and Marco Fondi - 06/04/2018
%   George diCenzo and Marco Fondi - updated - 12/11/2018
%

%% Check input variables

% Ensure there are enough inputs
assert(nargin >= 2,'This function requires atleast two inputs - a model and tnseq data');

% Set default epsilon value
if nargin < 3
    epsilon = 1.1e-6;
elseif isempty(epsilon)
    epsilon = 1.1e-6;
end

% Set default binning thresholds
if nargin < 4
    binThresh = {'3.5';'2.5';'1.5'};
elseif isempty(binThresh)
    binThresh = {'3.5';'2.5';'1.5'};
end

% Set growth threshold
if nargin < 5
    growthFrac = 0.5;
elseif isempty(growthFrac)
    growthFrac = 0.5;
end
solutionOriginal = optimizeCbModel(model,'max');
growthThresh = growthFrac * solutionOriginal.f;

% Is there RNA-seq data
if nargin < 6
    rnaPresent = 0;
elseif isempty(rnaseq)
    rnaPresent = 0;
else
    rnaPresent = 1;
end

% Check if there is RNA-seq data
if ~isempty(expressThresh)
    assert(~isempty(rnaseq), ...
        'RNA-seq data must be provided if expressThresh is set');
end

% Set default RNA-seq exprssion threshold
if ~isempty(rnaseq)
    if nargin < 7
        expressThresh = sum(cell2mat(rnaseq(:,1))) / 5000;
    elseif isempty(expressThresh)
        expressThresh = sum(cell2mat(rnaseq(:,1))) / 5000;
    end
end

% Set default to keep highly expressed genes
if nargin < 8
    keepHigh = 0;
elseif isempty(keepHigh)
    keepHigh = 0;
end

% Set default for deadends
if nargin < 9
    deadends = 1;
elseif isempty(deadends)
    deadends = 1;
end

%% Fix the model as necessary

if ~isfield(model, 'rxnGeneMat')
    model = tncore_fix(model);
elseif ~isfield(model, 'rules')
    model = tncore_fix(model);
elseif ~isfield(model, 'grRules')
    model = tncore_fix(model);
end

%% Refine the model to remove dead ends

deadEndMetabolites = detectDeadEnds(model, true);

if ~isempty(deadEndMetabolites)
    test = 1;
    while test > 0
        deadEndMetabolites = detectDeadEnds(model, true);
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
            model = tncore_remove_reactions(model,deadEndReactions);
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

%% Remove unnecessary unknowns from the model

% Find gene index number for Unknown
unknownID = findGeneIDs(model,'Unknown');

% Replace the unknowns that are not needed
A = [' | x(' num2str(unknownID) ')'];
B = ['x(' num2str(unknownID) ') | '];
model.grRules = strrep(model.grRules,' or Unknown','');
model.grRules = strrep(model.grRules,'Unknown or ','');
if isfield(model, 'rxnNotes')
    model.rxnNotes = strrep(model.rxnNotes,' or Unknown','');
    model.rxnNotes = strrep(model.rxnNotes,'Unknown or ','');
end
model.rules = strrep(model.rules,A,'');
model.rules = strrep(model.rules,B,'');

% Fix the rxnGeneMat
rxnGeneMatNewFull = cell(length(model.rxns),length(model.genes));

for n = 1:length(model.rxns)
    if ~isempty(model.rules{n})
        rulesTemp = strrep(model.rules{n}, '&', '|');
        rules = strsplit(rulesTemp, '|');
        for m = 1:length(rules)
            rules{m} = strrep(rules{m}, '(', '');
            rules{m} = strrep(rules{m}, ')', '');
            rules{m} = strrep(rules{m}, 'x', '');
            rules{m} = strrep(rules{m}, ' ', '');
        end
        for m = 1:length(rules)
            rxnGeneMatNewFull{n,str2num(rules{m})} = 1;
        end
    end
end

temp = cellfun('isempty',rxnGeneMatNewFull);
rxnGeneMatNewFull(temp) = {0};
rxnGeneMatNewFull = cell2mat(rxnGeneMatNewFull);
rxnGeneMatNew = sparse(double(rxnGeneMatNewFull));
model.rxnGeneMat = rxnGeneMatNew;

% Save the reduced full model
reducedModel = model;

%% Parse the tnseq input

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
modelTnseq = cell(length(model.genes),2);
for n = 1:length(model.genes)
    ID = strmatch(model.genes{n},tnseq(:,2));
    if isempty(ID)
        modelTnseq{n,2} = model.genes{n};
        modelTnseq{n,1} = max(cell2mat(tnseq(:,1)));
    else
        modelTnseq{n,2} = model.genes{n};
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

% Bin the data
modelTnseqGimme = modelTnseq;
for n = 1:length(modelTnseqGimme)
    if modelTnseqGimme{n,1} < medianValue - (str2num(binThresh{1}) * stdev)
        modelTnseqGimme{n,1} = 1001;
    elseif modelTnseqGimme{n,1} < medianValue - (str2num(binThresh{2}) * stdev)
        modelTnseqGimme{n,1} = 999.9;
    elseif modelTnseqGimme{n,1} < medianValue - (str2num(binThresh{3}) * stdev)
        modelTnseqGimme{n,1} = 990;
    elseif modelTnseqGimme{n,1} >= medianValue - (str2num(binThresh{3}) * stdev)
        modelTnseqGimme{n,1} = 0;
    end
end
essThresh = 1000;

%% Protect model genes essential for growth above the threshold

% Find essential genes
[~,grRateKO] = singleGeneDeletion(model);

% Protect the essential genes
lethal = grRateKO < growthThresh;
modelTnseqGimmeTemp = modelTnseqGimme;
modelTnseqGimmeTemp(lethal, 1) = num2cell(1001);
for n = 1:length(modelTnseqGimme)
    if strmatch(modelTnseqGimme{n,2}, 'Unknown', 'exact')
    elseif strmatch(modelTnseqGimme{n,2}, 'Spontaneous', 'exact')
    else
        modelTnseqGimme{n,1} = modelTnseqGimmeTemp{n,1};
    end
end

%% Protect model genes expressed above the threshold

if rnaPresent == 1

    % Pull out rnaseq data for model genes
    modelRNAseq = cell(length(model.genes),2);
    for n = 1:length(model.genes)
        ID = strmatch(model.genes{n},rnaseq(:,2));
        if isempty(ID)
            modelRNAseq{n,2} = model.genes{n};
            modelRNAseq{n,1} = min(cell2mat(rnaseq(:,1)));
        else
            modelRNAseq{n,2} = model.genes{n};
            modelRNAseq{n,1} = rnaseq{ID,1};
        end
    end
    
    % Determine expressed genes
    for n = 1:length(modelTnseqGimme)
        if modelRNAseq{n,1} >= expressThresh;
            if modelTnseqGimme{n,1} == 0
                modelTnseqGimme{n,1} = 1;
            end
        end
    end
    
end

%% Iterative reduction in the model

% Force objective to be full
sol = optimizeCbModel(model);
objFun = findRxnIDs(model, model.rxns(logical(model.c)));
origLB = model.lb(objFun);
model.lb(objFun) = sol.f * 0.99;

% Run fastcc
consistentRxns = fastcc(model, epsilon);
consistentRxns = model.rxns(consistentRxns);
nonConsistentRxns = setdiff(model.rxns, consistentRxns);
model = tncore_remove_reactions(model, nonConsistentRxns);
model = tncore_remove(model);

% Get genes from the least essential group and the inverse gene set
genesLow = {};
for n = 1:length(modelTnseqGimme)
    if modelTnseqGimme{n,1} == 0
        genesLow = vertcat(genesLow, modelTnseqGimme{n,2});
    end
end
genesHigh = setdiff(model.genes, genesLow);
genesHigh = intersect(genesHigh, model.genes);

% Fastcore to limit non-essential genes
[~, rxnList] = findRxnsFromGenes(model, genesHigh, [], 1);
rxnList = unique(rxnList(:,1));
rxnIDs = findRxnIDs(model, rxnList);
coreRxns = fastcore(rxnIDs, model, epsilon);
coreRxns = model.rxns(coreRxns);
nonCoreRxns = setdiff(model.rxns, coreRxns);
model = tncore_remove_reactions(model, nonCoreRxns);
model = tncore_remove(model);

if rnaPresent == 1
    % Get genes from the highly expressed least essential group
    genesLow = {};
    for n = 1:length(modelTnseqGimme)
        if modelTnseqGimme{n,1} == 1 || modelTnseqGimme{n,1} == 0
            genesLow = vertcat(genesLow, modelTnseqGimme{n,2});
        end
    end
    genesHigh = setdiff(model.genes, genesLow);
    genesHigh = intersect(genesHigh, model.genes);
    
    % Fastcore to limit non-essential genes
    [~, rxnList] = findRxnsFromGenes(model, genesHigh, [], 1);
    rxnList = unique(rxnList(:,1));
    rxnIDs = findRxnIDs(model, rxnList);
    coreRxns = fastcore(rxnIDs, model, epsilon);
    coreRxns = model.rxns(coreRxns);
    nonCoreRxns = setdiff(model.rxns, coreRxns);
    model = tncore_remove_reactions(model, nonCoreRxns);
    model = tncore_remove(model);
end

% Get genes from the second least essential group and the inverse gene set
genesLow = {};
for n = 1:length(modelTnseqGimme)
    if modelTnseqGimme{n,1} == 0 || modelTnseqGimme{n,1} == 990 || ...
            modelTnseqGimme{n,1} == 1
        genesLow = vertcat(genesLow, modelTnseqGimme{n,2});
    end
end
genesHigh = setdiff(model.genes, genesLow);
genesHigh = intersect(genesHigh, model.genes);

% Fastcore to limit non-essential genes
[~, rxnList] = findRxnsFromGenes(model, genesHigh, [], 1);
rxnList = unique(rxnList(:,1));
rxnIDs = findRxnIDs(model, rxnList);
coreRxns = fastcore(rxnIDs, model, epsilon);
coreRxns = model.rxns(coreRxns);
nonCoreRxns = setdiff(model.rxns, coreRxns);
model = tncore_remove_reactions(model, nonCoreRxns);
model = tncore_remove(model);

% Get genes from the third least essential group and the inverse gene set
genesLow = {};
for n = 1:length(modelTnseqGimme)
    if modelTnseqGimme{n,1} == 0 || modelTnseqGimme{n,1} == 990 || ...
            modelTnseqGimme{n,1} == 999.9 || modelTnseqGimme{n,1} == 1
        genesLow = vertcat(genesLow, modelTnseqGimme{n,2});
    end
end
genesHigh = setdiff(model.genes, genesLow);
genesHigh = intersect(genesHigh, model.genes);

% Fastcore to limit non-essential genes
[~, rxnList] = findRxnsFromGenes(model, genesHigh, [], 1);
rxnList = unique(rxnList(:,1));
rxnIDs = findRxnIDs(model, rxnList);
coreRxns = fastcore(rxnIDs, model, epsilon);
coreRxns = model.rxns(coreRxns);
nonCoreRxns = setdiff(model.rxns, coreRxns);
model = tncore_remove_reactions(model, nonCoreRxns);
model = tncore_remove(model);

% No longer force flux through objective
objFun = findRxnIDs(model, model.rxns(logical(model.c)));
model.lb(objFun) = origLB;

%% Remove unnecessary, low essential complexes

% Remove gene associations that are superseded by a more essential one
for n = 1:length(model.rxns)
    genesToKeep = {};
    genesToKeepTemp = {};
    genesToRemove = {};
    if strfind(model.grRules{n}, ' or ')
        if strfind(model.grRules{n}, ' and ')
            tempGrRules = strrep(model.grRules{n}, ' and ', ' or ');
            grRulesSplit = transpose(strsplit(tempGrRules, ' or '));
            for m = 1:length(grRulesSplit)
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, '(', '');
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, ')', '');
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, ' ', '');
                grRulesSplit{m,2} = modelTnseqGimme{...
                    strmatch(grRulesSplit{m,1}, modelTnseqGimme(:,2), 'exact'), 1};
            end
            grRulesSplit = sortrows(grRulesSplit, -2);
            if max(cell2mat(grRulesSplit(:,2))) > 0
                maxGrValue = max(cell2mat(grRulesSplit(:,2)));
                genesToKeep = {};
                genesToKeep_inv = {};
                for m = 1:length(grRulesSplit)
                    if grRulesSplit{m,2} == maxGrValue
                        genesToKeep = vertcat(genesToKeep, grRulesSplit{m,1});
                    else
                        genesToKeep_inv = vertcat(genesToKeep_inv, grRulesSplit{m,1});
                    end
                end
                genesToKeep_invTemporary = {};
                for m = 1:length(genesToKeep_inv)
                    [~, ~, constrRxnNames] = deleteModelGenes(model, genesToKeep_inv{m,1});
                    if strmatch(model.rxns{n}, constrRxnNames, 'exact')
                    else
                        genesToKeep_invTemporary = vertcat(genesToKeep_invTemporary, genesToKeep_inv(m,:));
                    end
                end
                genesToKeep_inv = genesToKeep_invTemporary;
                [~,~,constrRxns] = deleteModelGenes(model, genesToKeep_inv);
                if ~isempty(constrRxns)
                    if strmatch(model.rxns{n}, constrRxns, 'exact');
                        genesToKeepTemp = genesToKeep;
                        for m = 1:length(genesToKeep_inv)
                            genesToKeepTemp = vertcat(genesToKeepTemp, genesToKeep_inv{m,1});
                            genesToKeep_invTemp = setdiff(genesToKeep_inv(:,1), genesToKeepTemp);
                            [~,~,constrRxns] = deleteModelGenes(model, genesToKeep_invTemp);
                            if isempty(constrRxns)
                                tempGenes = vertcat(genesToKeep_inv{m,1}, genesToKeep_invTemp);
                                [~,~,constrRxns] = deleteModelGenes(model, tempGenes);
                                if isempty(constrRxns)
                                    genesToKeepTemp = genesToKeepTemp(1:end-1);
                                elseif ~strmatch(model.rxns{n}, constrRxns, 'exact')
                                    genesToKeepTemp = genesToKeepTemp(1:end-1);
                                end
                            elseif ~strmatch(model.rxns{n}, constrRxns, 'exact')
                                tempGenes = vertcat(genesToKeep, genesToKeep_invTemp);
                                [~,~,constrRxns] = deleteModelGenes(model, tempGenes);
                                if isempty(constrRxns)
                                    genesToKeepTemp = genesToKeepTemp(1:end-1);
                                elseif ~strmatch(model.rxns{n}, constrRxns, 'exact')
                                    genesToKeepTemp = genesToKeepTemp(1:end-1);
                                end
                            elseif strmatch(model.rxns{n}, constrRxns, 'exact')
                                tempGenes = vertcat(genesToKeep, genesToKeep_invTemp);
                                [~,~,constrRxns] = deleteModelGenes(model, tempGenes);
                                if isempty(constrRxns)
                                    genesToKeepTemp = genesToKeepTemp(1:end-1);
                                elseif ~strmatch(model.rxns{n}, constrRxns, 'exact')
                                    genesToKeepTemp = genesToKeepTemp(1:end-1);
                                end
                            end
                        end
                        genesToKeep = genesToKeepTemp;
                        genesToRemove = setdiff(genesToKeep_inv, genesToKeep);
                        if ~isempty(genesToRemove)
                            model.rxnGeneMat = full(model.rxnGeneMat);
                            for m = 1:length(genesToRemove)
                                model.grRules{n} = strrep(model.grRules{n}, ...
                                    genesToRemove{m}, 'toDelete');
                                genePos = strmatch(genesToRemove{m}, model.genes, 'exact');
                                ruleToRemove = ['x(' num2str(genePos) ')'];
                                if strmatch('toDelete', model.genes, 'exact')
                                    ruleToAdd = ['x(' num2str(length(model.genes)) ')'];
                                else
                                    ruleToAdd = ['x(' num2str(length(model.genes)+1) ')'];
                                end
                                model.rules{n} = strrep(model.rules{n}, ...
                                    ruleToRemove, ruleToAdd);
                                model.rxnGeneMat(n,genePos) = 0;
                            end
                            if strmatch('toDelete', model.genes, 'exact')
                                model.rxnGeneMat(n,length(model.genes)) = 1;
                            else
                                model.genes{end+1,1} = 'toDelete';
                                model.rxnGeneMat(n,length(model.genes)) = 1;
                            end
                        end
                    end
                end
            end
        end
    end
end
model.rxnGeneMat = num2cell(model.rxnGeneMat);
temp = cellfun('isempty',model.rxnGeneMat);
model.rxnGeneMat(temp) = {0};
model.rxnGeneMat = cell2mat(model.rxnGeneMat);
model.rxnGeneMat = sparse(double(model.rxnGeneMat));

% Fix gene list
model = tncore_remove(model);

%% Update the tnseq table

% Pull out tnseq data for new smaller model genes
modelTnseqGimmeNew = cell(length(model.genes),2);
for n = 1:length(model.genes)
    if strmatch('toDelete', model.genes{n}, 'exact')
        modelTnseqGimmeNew{n,2} = model.genes{n};
        modelTnseqGimmeNew{n,1} = -99000;
    else        
        ID = strmatch(model.genes{n},modelTnseqGimme(:,2));
        modelTnseqGimmeNew{n,2} = model.genes{n};
        modelTnseqGimmeNew{n,1} = modelTnseqGimme{ID,1};
    end
end
modelTnseqGimme = modelTnseqGimmeNew;

%% Run GIMME

% Get solver
global CBT_LP_SOLVER
solver = CBT_LP_SOLVER;
if isempty(solver)
    global CBTLPSOLVER
    solver = CBTLPSOLVER;
end

% Start tiger
start_tiger(solver);
tigerModelIn = cobra_to_tiger(model);

% Prepare expression data
modelRNAseqExp = cell2mat(modelTnseqGimme(:,1));

% Run GIMME
[geneStates, genesOut, tigerSol, tigerModelOut] = tncore_gimme(tigerModelIn, ...
    modelRNAseqExp, essThresh, modelTnseqGimme(:,2), growthFrac);

%% Build context specific model in COBRA

% Identify genes off in GIMME and delete from the reduced model
isPresent = logical(geneStates);
genesToDelete = genesOut(~isPresent);
[testModel, ~, constrRxnNames] = deleteModelGenes(model, genesToDelete);

% Test if the resulting model grows
testSol = optimizeCbModel(testModel);

if testSol.f < growthThresh
    
    % If model does not grow, add back in additional genes
    for n = 1:length(genesToDelete)
        geneID = findGeneIDs(model, genesToDelete{n,1});
        genesToDelete{n,2} = modelTnseqGimme{geneID,1};
    end
    
    genesToDelete = sortrows(genesToDelete,2);
    tempModel = model;
    rxnsToRemove = {};
    genesToDeleteNew = {};

    for n = 1:length(genesToDelete)
        [~,~,constrRxnNames] = deleteModelGenes(tempModel, genesToDelete{n,1});
        testModel = tncore_remove_reactions(tempModel, constrRxnNames);
        testSol = optimizeCbModel(testModel);
        if testSol.f >= growthThresh
            [tempModel,~,constrRxnNames] = deleteModelGenes(tempModel, genesToDelete{n,1});
            rxnsToRemove = vertcat(rxnsToRemove, constrRxnNames);
            genesToDeleteNew = vertcat(genesToDeleteNew, genesToDelete{n,1});
        end
    end
    
    % Save the expressed genes
    if rnaPresent == 1
        if keepHigh == 1
            expressed = cell2mat(modelRNAseq(:,1)) >= expressThresh;
            expressedGenes = modelRNAseq(expressed,2);
            expressedGenes = intersect(expressedGenes, model.genes);
            genesToDeleteNew = setdiff(genesToDeleteNew, expressedGenes);
        end
    end
    
    [tempModel,~,constrRxnNames] = deleteModelGenes(model, genesToDeleteNew);
    tempModel = tncore_remove_reactions(tempModel, constrRxnNames);
    contextModel = tncore_delete(tempModel);

else
    
    % If model growns
    
    % Save the expressed genes
    if rnaPresent == 1
        if keepHigh == 1
            expressed = cell2mat(modelRNAseq(:,1)) >= expressThresh;
            expressedGenes = modelRNAseq(expressed,2);
            expressedGenes = intersect(expressedGenes, model.genes);
            genesToDeleteNew = setdiff(genesToDeleteNew, expressedGenes);
        end
    end
    
    [tempModel,~,constrRxnNames] = deleteModelGenes(model, genesToDelete);
    tempModel = tncore_remove_reactions(tempModel, constrRxnNames);
    contextModel = tncore_delete(tempModel);
    
end

%% Refine the model to remove dead ends

if deadends == 1
    
    model = contextModel;
    deadEndMetabolites = detectDeadEnds(model, true);
    if ~isempty(deadEndMetabolites)
        test = 1;
        while test > 0
            deadEndMetabolites = detectDeadEnds(model, true);
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
                model = tncore_remove_reactions(model,deadEndReactions);
                if length(model.mets) < size(model.S, 1)
                    model.S(deadEndMetabolites,:) = [];
                    model.b(deadEndMetabolites,:) = [];
                end
                clear deadEndMetabolites
                clear deadEndMetabolites2
            end
        end
        contextModel = tncore_remove(model);
    end
    
end
