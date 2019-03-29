function [contextModel, reducedModel] = tncore_rna(model, rnaseq, ...
    expressThresh, growthFrac, coreGenes, tnseq, essThresh, forceHigh, epsilon)

%
% An adaptation of the tncore_core script for production of a context-
% specific model using RNA-seq data and optionally also integrating a set 
% of core genes determined via Tn-seq.
%
% USAGE
%   [contextModel, reducedModel] = tncore_rna(model, rnaseq, expressThresh, ...
%       growthFrac, coreGenes, tnseq, essThresh, forceHigh, epsilon)
%
% INPUTS
%   model           COBRA model structure including gene-reaction
%                   associations
%   rnaseq          RNA-seq data for all genes in the organism. First
%                   column contains the RNA-seq data, the second column
%                   contains the gene names
%
% OPTIONAL INPUTS
%   expressThresh   The threshold for a gene to be considered expressed
%                   (Default = 0.02% of the sum of all expression values)
%   growthFrac      The minimum allowable objective flux in the generated
%                   core model, as a fraction of growth of the input model
%                   (Default = 0.5) 
%   coreGenes       A list of all genes to be protected during context-
%                   specific model generation. Can either be empty, a cell 
%                   array of genes to protect, or {1} to indicate 
%                   to automatically determine the core as the essential 
%                   genes from the Tn-seq data (Default = {})
%   tnseq           Tn-seq data for all genes in the organism. First 
%                   column contains the Tn-seq scores, the second column 
%                   contains the gene names. Tn-seq data is necessary if 
%                   coreGenes is set to {1}, and Tn-seq data will be 
%                   ignored if coreGenes is set differently. (Default = {}) 
%   essThresh       The threshold (in number of standard deviations from
%                   the mean of the log of the Tn-seq data) for a gene to 
%                   be considered essential. (Default = 3.5)
%   forceHigh       An option to choose a workflow that favours inclusion of
%                   a pathway/complex with a highly expressed gene even if
%                   an alternative pathway/complex with no highly expressed
%                   gene, but average overall expression, is present. Use 1
%                   for yes, and 0 for no. (Default = 0)
%   epsilon         The epsilon value (the minimum flux rate that is
%                   considered non-zero) to be used during the running of
%                   FASTCC and FASTCORE when using forceHigh.
%                   (Default = 1.1e-6)
%
% OUTPUTS
%   contextModel    The context-specific  model based on the RNA-seq data
%                   (considering also the Tn-seq data, if provided)
%   reducedModel    The input model with dead-ends removed and with
%                   unnecessary 'Unknown' GPRs removed
%
% AUTHORS
%   George diCenzo and Marco Fondi - 10/12/2018
%

%% Check input variables

% Ensure there are enough inputs
assert(nargin >= 2,'This function requires atleast two inputs');

% Set default RNA-seq exprssion threshold
if nargin < 3
    expressThresh = sum(cell2mat(rnaseq(:,1))) / 5000;
elseif isempty(expressThresh)
    expressThresh = sum(cell2mat(rnaseq(:,1))) / 5000;
end

% Set growth threshold
if nargin < 4
    growthFrac = 0.5;
elseif isempty(growthFrac)
    growthFrac = 0.5;
end
solutionOriginal = optimizeCbModel(model,'max');
growthThresh = growthFrac * solutionOriginal.f;

% Set default coreGenes
if nargin < 5
    coreGenes = {};
elseif isempty(coreGenes)
    coreGenes = {};
end

% Is there Tn-seq data
if ~isempty(coreGenes)
    if coreGenes{1} == 1
        assert(nargin >= 6,'Tn-seq data must be provided if coreGenes = {1}');
        assert(~isempty(tnseq),'Tn-seq data must be provided if coreGenes = {1}');
    end
end

% Set default essThresh
if nargin < 7
    essThresh = 3.5;
elseif isempty(essThresh)
    essThresh = 3.5;
end

% Set default forceHigh
if nargin < 8
    forceHigh = 0;
elseif isempty(essThresh)
    forceHigh = 0;
end

% Set default epsilon value
if nargin < 9
    epsilon = 1.1e-6;
elseif isempty(epsilon)
    epsilon = 1.1e-6;
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

%% Parse the RNA-seq input

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

%% Turn on essential model genes

% Find essential genes
[~,grRateKO] = singleGeneDeletion(model);

% Expand the protected gene set with the essential genes
lethal = grRateKO < growthThresh;
expressTemp = modelRNAseq(1:end,1);
expressTemp(lethal) = {2 * expressThresh};
modelRNAseq(1:end,1) = expressTemp;

%% Turn on core genes / essential Tn-seq genes

% Convert tnseq data to log scale
if cell2mat(coreGenes) == 1
    
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

end

% Pull out tnseq data for model genes
if ~isempty(coreGenes)
    
    if cell2mat(coreGenes) == 1
        
        modelTnseq = cell(length(model.genes),1);
        for n = 1:length(model.genes)
            ID = strmatch(model.genes{n},tnseq(:,2));
            if isempty(ID)
                modelTnseq{n,1} = max(cell2mat(tnseq(:,1)));
            else
                modelTnseq{n,1} = tnseq{ID,1};
            end
        end
        isCore = cell2mat(modelTnseq) < medianValue - (essThresh * stdev);
        expressTemp = modelRNAseq(1:end,1);
        expressTemp(isCore) = {2 * expressThresh};
        modelRNAseq(1:end,1) = expressTemp;
        
    else
        
        geneIDs = num2cell(findGeneIDs(model, coreGenes));
        for n = 1:length(coreGenes)
            if geneIDs{n} ~= 0
                modelRNAseq{geneIDs{n},1} = 2 * expressThresh;
            end
        end  
        
    end
    
end

%% Iterative reduction in the model

if forceHigh == 1
    
    % Force objective to be full
    sol = optimizeCbModel(model);
    objFun = findRxnIDs(model, model.rxns(logical(model.c)));
    origLB = model.lb(objFun);
    model.lb(objFun) = sol.f;
    
    % Run fastcc
    consistentRxns = fastcc(model, epsilon);
    consistentRxns = model.rxns(consistentRxns);
    nonConsistentRxns = setdiff(model.rxns, consistentRxns);
    model = tncore_remove_reactions(model, nonConsistentRxns);
    model = tncore_remove(model);
    
    % Get genes from the low expression group and the inverse gene set
    genesLow = {};
    for n = 1:length(modelRNAseq)
        if modelRNAseq{n,1} < expressThresh
            genesLow = vertcat(genesLow, modelRNAseq{n,2});
        end
    end
    genesHigh = setdiff(model.genes, genesLow);
    genesHigh = intersect(genesHigh, model.genes);
    
    % Fastcore to limit low expression group
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
    
end

%% Remove unnecessary, low expressed complexes

if forceHigh == 1
    
    % Remove gene associations that are superseded by a highly expressed one
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
                    grRulesSplit{m,2} = modelRNAseq{...
                        strmatch(grRulesSplit{m,1}, modelRNAseq(:,2), 'exact'), 1};
                end
                grRulesSplit = sortrows(grRulesSplit, -2);
                if max(cell2mat(grRulesSplit(:,2))) >= expressThresh
                    genesToKeep = {};
                    genesToKeep_inv = {};
                    for m = 1:length(grRulesSplit)
                        if grRulesSplit{m,2} >= expressThresh
                            genesToKeep = vertcat(genesToKeep, grRulesSplit{m,1});
                        else
                            genesToKeep_inv = vertcat(genesToKeep_inv, grRulesSplit{m,1});
                        end
                    end
                    [~,~,constrRxns] = deleteModelGenes(model, genesToKeep_inv);
                    if ~isempty(constrRxns)
                        if strmatch(model.rxns{n}, constrRxns, 'exact');
                            genesToKeepTemp = genesToKeep;
                            for m = length(genesToKeep)+1:length(grRulesSplit)
                                genesToKeepTemp = vertcat(genesToKeepTemp, grRulesSplit{m,1});
                                genesToKeep_invTemp = setdiff(grRulesSplit(:,1), genesToKeepTemp);
                                [~,~,constrRxns] = deleteModelGenes(model, genesToKeep_invTemp);
                                if isempty(constrRxns)
                                    tempGenes = vertcat(genesToKeep, genesToKeep_invTemp);
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
                                    elseif ~strmach(model.rxns{n}, constrRxns, 'exact')
                                        genesToKeepTemp = genesToKeepTemp(1:end-1);
                                    end
                                end
                            end
                            genesToKeep = genesToKeepTemp;
                            genesToRemove = setdiff(grRulesSplit(:,1), genesToKeep);
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
    
end

%% Update the rnaseq table

if forceHigh == 1
    
    % Pull out tnseq data for new smaller model genes
    modelRNAseqNew = cell(length(model.genes),2);
    for n = 1:length(model.genes)
        if strmatch('toDelete', model.genes{n}, 'exact')
            modelRNAseqNew{n,2} = model.genes{n};
            modelRNAseqNew{n,1} = -99000;
        else
            ID = strmatch(model.genes{n},modelRNAseq(:,2));
            modelRNAseqNew{n,2} = model.genes{n};
            modelRNAseqNew{n,1} = modelRNAseq{ID,1};
        end
    end
    modelRNAseq = modelRNAseqNew;
    
end

%% Run GIMME in TIGER

% Get solver
global CBT_LP_SOLVER
solver = CBT_LP_SOLVER;
if isempty(solver)
    global CBTLPSOLVER
    solver = CBTLPSOLVER;
end

% Start tiger
start_tiger('gurobi');
tigerModelIn = cobra_to_tiger(model);

% Prepare expression data
modelRNAseqExp = cell2mat(modelRNAseq(:,1));

% Run GIMME
[geneStates, genesOut, tigerSol, tigerModelOut] = tncore_gimme(tigerModelIn, ...
    modelRNAseqExp, expressThresh, modelRNAseq(:,2), growthFrac);

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
        genesToDelete{n,2} = modelRNAseq{geneID,1};
    end
    
    genesToDelete = sortrows(genesToDelete,2);
    tempModel = model;
    rxnsToRemove = {};
    geneToDeleteNew = {};

    for n = 1:length(genesToDelete)
        [~,~,constrRxnNames] = deleteModelGenes(tempModel, genesToDelete{n,1});
        testModel = tncore_remove_reactions(tempModel, constrRxnNames);
        testSol = optimizeCbModel(testModel);
        if testSol.f >= growthThresh
            [tempModel,~,constrRxnNames] = deleteModelGenes(tempModel, genesToDelete{n,1});
            rxnsToRemove = vertcat(rxnsToRemove, constrRxnNames);
            geneToDeleteNew = vertcat(geneToDeleteNew, genesToDelete{n,1});
        end
    end
    
    [tempModel,~,constrRxnNames] = deleteModelGenes(model, geneToDeleteNew);
    tempModel = tncore_remove_reactions(tempModel, constrRxnNames);
    contextModel = tncore_delete(tempModel);

else
    
    % If model grows, process is done
    contextModel = tncore_remove_reactions(testModel, constrRxnNames);
    contextModel = tncore_delete(contextModel);
    
end
