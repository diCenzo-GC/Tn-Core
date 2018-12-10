function [coreGeneVar, coreRxnVar, coreGrowthVar, coreModel, ... 
            coreModelFast, reducedModel, coreGeneCat, genePresence, ...
            rxnPresence] = tncore_redundancy(model, iters, growthFrac, ...
            tnseq, coreGenes, essThresh, binThresh, rnaseq, expressThresh, ...
            method)

%
% Generates randomized core models; identifies the core model most 
% consistent with provided Tn-seq data and the model with the highest
% objective flux; can optionally include both Tn-seq and RNA-seq data for 
% generation of context specific metabolic models.
%
% USAGE
%   [coreGeneVar, coreRxnVar, coreGrowthVar, coreModel, coreModelFast,...
%       reducedModel, coreGeneCat, genePresence, rxnPresence] = ...
%       tncore_redundancy(model, iters, growthFrac, tnseq, coreGenes, ...
%       essThresh, binThresh, rnaseq, expressThresh, method)
%
% INPUTS
%   model           COBRA model structure including gene-reaction
%                   associations
%
% OPTIONAL INPUTS
%   iters           Number of random core models to generate 
%                   (Default = 1000)
%   growthFrac      The minimum allowable objective flux in the generated
%                   core model, as a fraction of growth of the input model
%                   (Default = 0.5) 
%   tnseq           Tn-seq data for all genes in the organism. First 
%                   column contains the Tn-seq scores, the second column 
%                   contains the gene names (Default = {})
%   coreGenes       A list of all essential genes to be protected during 
%                   random model generation. Can either be empty, a cell 
%                   array of genes to protect, or {1} to indicate 
%                   to automatically determine the core as the essential 
%                   genes from the Tn-seq data (Default = {})
%   essThresh       The threshold (in number of standard deviations from
%                   the mean of the log of the Tn-seq data) for a gene to 
%                   be considered essential. Used if coreGenes is set to 
%                   {1}. (Default = 3.5)
%   binThresh       An array containing three values (number of standard
%                   deviations away from the mean of the log of the Tn-seq 
%                   data) to use for setting the limits of each bin. 
%                   (Default = {'3.5'; '2.5'; '1.5'}) 
%   rnaseq          RNA-seq data for all genes in the organism. First
%                   column contains the RNA-seq data, the second column
%                   contains the gene names (Default = {})
%   expressThresh   The threshold for a gene to be considered expressed
%                   (Default = 0.02% of the sum of all expression values)
%   method          Indicates if FBA or MOMA should be used during core
%                   model generation. Use 1 for FBA and 2 for MOMA (Default
%                   = 1)
%
% OUTPUTS
%   coreGeneVar     A gene list indicating the genes that are present (gene
%                   name given) or deleted (gene name followed by
%                   _deleted) in each core model. Genes are in rows, 
%                   models are in columns.
%   coreRxnVar      A reaction list indicating the reactions that are
%                   present(reaction name given) or deleted (reaction name
%                   followed by _deleted) in each core model. Reactions
%                   are in rows, models are in columns.
%   coreGrowthVar   The objective flux through each of the core models.
%                   Order corresponds to the same order of models in the
%                   coreGeneVar and coreRxnVar.
%   coreModel       The core model most consistent with the Tn-seq data.
%   coreModelFast   The core model with the highest objective flux.
%   reducedModel    The input model with dead-ends removed and with
%                   unnecessary 'Unknown' GPRs removed.
%   coreGeneCat     The number of core model genes grouped into each
%                   essentiality category.
%   genePresence    A binary matrix of gene presence/absence in each core
%                   model. Genes are in rows and models in columns.
%   rxnPresence     A binary matrix of reaction presence/absence in each
%                   core model. Reactions are in rows and models in
%                   columns.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 01/11/2017
%   George diCenzo and Marco Fondi - updated - 06/04/2018
%   George diCenzo and Marco Fondi - updated - 12/11/2018
%

%% Check input variables

% Ensure there are enough inputs
assert(nargin >= 1,'This function requires atleast one input');

% Set default iterations
if nargin < 2
    iters = 1000;
elseif isempty(iters)
    iters = 1000;
end

% Set default growth threshold
if nargin < 3
    growthFrac = 0.5;
elseif isempty(growthFrac)
    growthFrac = 0.5;
end
solutionOriginal = optimizeCbModel(model,'max');
growthThresh = growthFrac * solutionOriginal.f;

% Check if there is tnseq data
if nargin < 4
    tnPresent = 0;
elseif isempty(tnseq)
    tnPresent = 0;
else
    tnPresent = 1;
end

% Set default coreGenes
if nargin < 5
    coreGenes = {};
elseif isempty(coreGenes)
    coreGenes = {};
end

% Set the default essentiality threshold
if nargin < 6
    essThresh = 3.5;
elseif isempty(essThresh)
    essThresh = 3.5;
end

% Set default binning thresholds
if nargin < 7
    binThresh = {'3.5';'2.5';'1.5'};
elseif isempty(binThresh)
    binThresh = {'3.5';'2.5';'1.5'};
end

% Is there RNA-seq data
if nargin < 8
    rnaPresent = 0;
elseif isempty(rnaseq)
    rnaPresent = 0;
else
    rnaPresent = 1;
end

% Set default RNA-seq exprssion threshold
if rnaPresent == 1
    if nargin < 9
        expressThresh = sum(cell2mat(rnaseq(:,1))) / 5000;
    elseif isempty(expressThresh)
        expressThresh = sum(cell2mat(rnaseq(:,1))) / 5000;
    end
end

% Set default method
if nargin < 10
    method = 1;
elseif isempty(method)
    method = 1;
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
            [deadEndReactions] = ...
                findRxnsFromMets(model,deadEndMetabolites3);
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

% Identify core genes present in the model
testCore = cell2mat(coreGenes);
if testCore ~= 1
    protectedGenes = intersect(coreGenes,model.genes);
    protectedGenesExpanded = protectedGenes;
else
    protectedGenesExpanded = {};
end

%% Parse the tnseq input

if tnPresent == 1
    
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

end

%% Set essential genes as core genes

if testCore == 1
    
    essentialTnseq = {};
    
    for n = 1:length(modelTnseq)
        if modelTnseq{n,1} < medianValue - (essThresh * stdev)
            essentialTnseq = vertcat(essentialTnseq,modelTnseq{n,2});
        end
    end
    
    protectedGenesExpanded = union(protectedGenesExpanded,essentialTnseq);
 
end

%% Parse the RNA-seq input

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
    expressThreshArray = cell(length(modelRNAseq),1);
    isEmpty = cellfun('isempty',expressThreshArray);
    expressThreshArray(isEmpty) = {expressThresh};
    
    isExpressed = cell2mat(modelRNAseq(:,1)) > cell2mat(expressThreshArray);
    genesExpressed = modelRNAseq(isExpressed,2);
    protectedGenesExpanded = union(protectedGenesExpanded,genesExpressed);
    
end

%% Identify essential model genes

% Find essential genes
[~,grRateKO,~,~,~] = singleGeneDeletion(model);

% Expand the protected gene set with the essential genes
lethal = grRateKO < growthThresh;
essGenes = model.genes(lethal);
protectedGenesExpanded = union(protectedGenesExpanded,essGenes);
nonProtectedGenes = setdiff(model.genes,protectedGenesExpanded);

%% Produce randomized core models

randomizedCoreGenes = {};
randomizedCoreGenes2 = {};
coreGeneCategories = {};
randomizedCoreGrowth = {};
randomizedCoreReactions = {};
genePresence = {};
rxnPresence = {};

if tnPresent == 0
    modelTnseq = [];
    medianValue = [];
    stdev = [];
end

global CBT_LP_SOLVER
solver = CBT_LP_SOLVER;
if isempty(solver)
    global CBTLPSOLVER
    solver = CBTLPSOLVER;
end

parfor n = 1:iters
    [randomCore,genesTnseq,geneGrouping,solution,reactions,rxnPresenceList,...
        genePresenceList] = tncore_randomize(model,growthThresh,...
        nonProtectedGenes,modelTnseq,medianValue,stdev,binThresh,solver,method);
    randomizedCoreGenes = horzcat(randomizedCoreGenes,randomCore.genes);
    randomizedCoreGenes2 = horzcat(randomizedCoreGenes2,genesTnseq);
    coreGeneCategories = horzcat(coreGeneCategories,geneGrouping);
    randomizedCoreGrowth = horzcat(randomizedCoreGrowth,num2cell(solution.f));
    randomizedCoreReactions = horzcat(randomizedCoreReactions,reactions);
    genePresence = horzcat(genePresence,genePresenceList);
    rxnPresence = horzcat(rxnPresence,rxnPresenceList);
end

coreGeneVar = randomizedCoreGenes;
coreRxnVar = randomizedCoreReactions;
coreGrowthVar = randomizedCoreGrowth;
coreGeneCat = coreGeneCategories;

%% Find the model most consistent with the tnseq data

if tnPresent == 1

    % Find models with most number of most essential genes
    maxGroup4 = max(cell2mat(coreGeneCategories(4,:)));
    maxGroup4IDs = {};

    for n = 1:iters
        if strmatch(maxGroup4,coreGeneCategories{4,n},'exact')
            maxGroup4IDs = vertcat(maxGroup4IDs,n);
        end
    end

    % Of those models, find models with the most number of second essential
    temp = {};

    for n = 1:length(maxGroup4IDs);
        temp = vertcat(temp,coreGeneCategories{3,maxGroup4IDs{n}});
    end

    maxGroup3 = max(cell2mat(temp));
    maxGroup3IDs = {};

    for n = 1:length(maxGroup4IDs)
        if strmatch(maxGroup3,coreGeneCategories{3,maxGroup4IDs{n}},'exact')
            maxGroup3IDs = vertcat(maxGroup3IDs,maxGroup4IDs{n});
        end
    end

    % Of those models, find models with the most number of third essential
    temp = {};

    for n = 1:length(maxGroup3IDs);
        temp = vertcat(temp,coreGeneCategories{2,maxGroup3IDs{n}});
    end

    maxGroup2 = max(cell2mat(temp));
    maxGroup2IDs = {};

    for n = 1:length(maxGroup3IDs)
        if strmatch(maxGroup2,coreGeneCategories{2,maxGroup3IDs{n}},'exact')
            maxGroup2IDs = vertcat(maxGroup2IDs,maxGroup3IDs{n});
        end
    end

    % If multiple models remain, chose final based on growth rate
    temp = {};

    if length(maxGroup2IDs) > 1
        for n = 1:length(maxGroup2IDs);
            temp = vertcat(temp,randomizedCoreGrowth{1,maxGroup2IDs{n}});
        end
        maxCoreGrowth = max(cell2mat(temp));
        maxCoreGrowthIDs = {};
        for n = 1:length(maxGroup2IDs)
            if strmatch(maxCoreGrowth,...
                    randomizedCoreGrowth{1,maxGroup2IDs{n}},'exact')
                maxCoreGrowthIDs = ...
                    vertcat(maxCoreGrowthIDs,maxGroup2IDs{n});
            end
        end
        coreModelID = maxCoreGrowthIDs{1,1};
    else
        coreModelID = maxGroup2IDs{1,1};
    end

end

%% Reproduce the core model

if tnPresent == 1
    
    % Identify the genes to delete
    geneIDs = cellfun('isempty',randomizedCoreGenes2(:,coreModelID));
    genesToDelete = model.genes(geneIDs);

    % Delete genes
    [coreModel,~,constrRxnNames,~] = ...
        deleteModelGenes(model,genesToDelete);

    % Identify non-essential unknown or non-enzymatic reactions
    y = 0;
    reactionsToRemove = {};
    smallModelSolution = optimizeCbModel(coreModel,'max');

    for n = 1:length(coreModel.grRules)
        if strmatch(coreModel.grRules{n,1},'Unknown','exact')
            x = 1;
        elseif strmatch(coreModel.grRules{n,1},'Spontaneous')
            x = 1;
        elseif isempty(coreModel.grRules{n,1})
            x = 1;
        else
            x = 0;
        end
        if x == 1
            testModel = coreModel;
            testModel = tncore_remove_reactions(testModel,coreModel.rxns{n,1});
            testModelSolution = optimizeCbModel(testModel,'max');
            if testModelSolution.f / smallModelSolution.f > 0.99
                y = y + 1;
                reactionsToRemove{y,1} = coreModel.rxns{n,1};
            end
        end
    end

    % Remove non-essential unknown or non-enzymatic reactions
    if ~isempty(reactionsToRemove)
        for n = 1:length(reactionsToRemove)
            testModel = tncore_remove_reactions(coreModel,reactionsToRemove{n,1});
            testModelSolution = optimizeCbModel(testModel,'max');
            if testModelSolution.f / smallModelSolution.f > 0.99
                coreModel = tncore_remove_reactions(coreModel,reactionsToRemove{n,1});
            end
            smallModelSolution = optimizeCbModel(coreModel,'max');
        end
    end

    % Delete reactions and fix model accordingly
    coreModel = tncore_remove_reactions(coreModel,constrRxnNames);
    coreModel = tncore_delete(coreModel);

end

%% Find the fastest growing core model

if tnPresent == 1
    
    % Find models with fastest growth rate
    maxGrowthRate = max(cell2mat(randomizedCoreGrowth));
    maxGrowthRateIDs = {};

    for n = 1:iters
        if strmatch(maxGrowthRate,randomizedCoreGrowth{1,n},'exact')
            maxGrowthRateIDs = vertcat(maxGrowthRateIDs,n);
        end
    end

    % If there are multiple fastest models, decide based on essentiality
    if length(maxGrowthRateIDs) > 1
        
        % Find models with most number of most essential genes
        temp = {};
        for n = 1:length(maxGrowthRateIDs);
            temp = vertcat(temp,coreGeneCategories{4,maxGrowthRateIDs{n}});
        end
        maxGroup4 = max(cell2mat(temp));
        maxGroup4IDs = {};
        for n = 1:length(maxGrowthRateIDs)
            if strmatch(maxGroup4,coreGeneCategories{4,maxGrowthRateIDs{n}},...
                   'exact')
                maxGroup4IDs = vertcat(maxGroup4IDs,maxGrowthRateIDs{n});
            end
        end
        
        % Find models with most number of most second essential genes
        temp = {};
        for n = 1:length(maxGroup4IDs);
            temp = vertcat(temp,coreGeneCategories{3,maxGroup4IDs{n}});
        end
        maxGroup3 = max(cell2mat(temp));
        maxGroup3IDs = {};
        for n = 1:length(maxGroup4IDs)
            if strmatch(maxGroup3,coreGeneCategories{3,maxGroup4IDs{n}},...
                    'exact')
                maxGroup3IDs = vertcat(maxGroup3IDs,maxGroup4IDs{n});
            end
        end
        
        % Find models with most number of most third essential genes
        temp = {};
        for n = 1:length(maxGroup3IDs);
            temp = vertcat(temp,coreGeneCategories{2,maxGroup3IDs{n}});
        end
        maxGroup2 = max(cell2mat(temp));
        maxGroup2IDs = {};
        for n = 1:length(maxGroup3IDs)
            if strmatch(maxGroup2,coreGeneCategories{2,maxGroup3IDs{n}},...
                    'exact')
                maxGroup2IDs = vertcat(maxGroup2IDs,maxGroup3IDs{n});
            end
        end
        
        coreModelFastID = maxGroup2IDs{1,1};
    else
        coreModelFastID = maxGrowthRateIDs{1,1};
    end

end

%% Reproduce the fast growing core model

if tnPresent == 1

    % Identify the genes to delete
    geneIDs = cellfun('isempty',randomizedCoreGenes2(:,coreModelFastID));
    genesToDelete = model.genes(geneIDs);

    % Delete genes
    [coreModelFast,~,constrRxnNames,~] = ...
        deleteModelGenes(model,genesToDelete);

    % Identify non-essential unknown or non-enzymatic reactions
    y = 0;
    reactionsToRemove = {};
    smallModelSolution = optimizeCbModel(coreModelFast,'max');

    for n = 1:length(coreModelFast.grRules)
        if strmatch(coreModelFast.grRules{n,1},'Unknown','exact')
            x = 1;
        elseif strmatch(coreModelFast.grRules{n,1},'Spontaneous')
            x = 1;
        elseif isempty(coreModelFast.grRules{n,1})
            x = 1;
        else
            x = 0;
        end
        if x == 1
            testModel = coreModelFast;
            testModel = tncore_remove_reactions(testModel,coreModelFast.rxns{n,1});
            testModelSolution = optimizeCbModel(testModel,'max');
            if testModelSolution.f / smallModelSolution.f > 0.99
                y = y + 1;
                reactionsToRemove{y,1} = coreModelFast.rxns{n,1};
            end
        end
    end

    % Remove non-essential unknown or non-enzymatic reactions
    if ~isempty(reactionsToRemove)
        for n = 1:length(reactionsToRemove)
            testModel = tncore_remove_reactions(coreModelFast,reactionsToRemove{n,1});
            testModelSolution = optimizeCbModel(testModel,'max');
            if testModelSolution.f / smallModelSolution.f > 0.99
                coreModelFast = tncore_remove_reactions(coreModelFast,reactionsToRemove{n,1});
            end
            smallModelSolution = optimizeCbModel(coreModelFast,'max');
        end
    end  

    % Delete reactions and fix model accordingly
    coreModelFast = tncore_remove_reactions(coreModelFast,constrRxnNames);
    coreModelFast = tncore_delete(coreModelFast);

end

%% Place holders for if there is no tnseq data

if tnPresent == 0
    coreModel = {};
    coreModelFast = {};
end

