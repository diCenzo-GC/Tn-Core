function [coreModel, genesTnseq, geneGrouping, solution, reactions, ...
            rxnPresence, genePresence] = tncore_randomize(model, ...
            growthThresh, nonProtected, modelTnseq, medianValue, stdev, ...
            binThresh, solver, method)

%
% Takes an input model, and produces a randomly prepared core model that
% produced an objective flux at least as large as the input threshold.
%
% USAGE
%   [coreModel, genesTnseq, geneGrouping, solution, reactions, ...
%        rxnPresence, genePresence] = tncore_randomize(model, ...
%        growthThresh, nonProtected, modelTnseq, medianValue, stdev, ...
%        solver, method)
%
% INPUTS
%   model           The model from which the core models are to be derived.
%
% OPTIONAL INPUTS
%   growthThresh    The minimum allowable objective flux in the generated
%                   core models (Default = 10% the objective flux of the 
%                   input model)
%   nonProtected    The list of genes that can be deleted during
%                   constructin of the core model (Default = genes whose
%                   deletion does not result in an objective flux below the
%                   growthThresh)
%   modelTnseq      Tn-seq data for all genes in the model. First 
%                   column contains the Tn-seq scores, the second column 
%                   contains the gene names (Default = {})
%   medianValue     Median of the log transformed Tn-seq data following
%                   removal of outliers (Default = {})
%   stdev           Standard deviation of the log transformed Tn-seq data 
%                   following removal of outliers (Default = {})
%   binThresh       An array containing three values (number of standard
%                   deviations away from the mean of the log of the Tn-seq 
%                   data) to use for setting the limits of each bin. 
%                   (Default = {'3.5'; '2.5'; '1.5'}) 
%   solver          LP solver to be used (Default = the currently set LP
%                   solver)
%   method:         Indicates if FBA or MOMA should be used during core
%                   model generation. Use 1 for FBA and 2 for MOMA (Default
%                   = 1)
%
% OUTPUT
%   coreModel       The core model produced from the input model.
%   genesTnseq      A list of the Tn-seq values for the genes included in
%                   the model.
%   geneGrouping    The number of core model genes grouped into each
%                   essentiality category.
%   solution        The objective flux rate of the core model.
%   reactions       A list indicating the reactions that are present
%                   (reaction name given) or deleted (reaction name
%                   followed by _deleted) in each core model.
%   rxnPresence:    A binary matrix of reaction presence/absence in the
%                   core model.
%   genePresence:   A binary matrix of gene presence/absence in the
%                   core model. 
%
% AUTHORS
%   George diCenzo and Marco Fondi - 01/11/2017
%   George diCenzo and Marco Fondi - updated - 06/04/2018
%   George diCenzo and Marco Fondi - updated - 10/12/2018
%

%% Check inputs

% Ensure there are enough inputs
assert(nargin >= 1,'This function requires atleast one input');

% Set default growthThresh
if nargin < 2
    solutionOriginal = optimizeCbModel(model,'max');
    growthThresh = 0.1 * solutionOriginal.f;
elseif isempty(growthThresh)
    solutionOriginal = optimizeCbModel(model,'max');
    growthThresh = 0.1 * solutionOriginal.f;
end

% Set default nonProtected list
if nargin < 3
    [grRatio] = singleGeneDeletion(model);
    grRatio = round(grRatio,2);
    isNotLethal = grRatio >= growthThresh;
    nonProtected = model.genes(isNotLethal);
elseif isempty(nonProtected)
    [grRatio] = singleGeneDeletion(model);
    grRatio = round(grRatio,2);
    isNotLethal = grRatio >= growthThresh;
    nonProtected = model.genes(isNotLethal);
end

% Determine if Tnseq data is present
if nargin < 4
    tnPresent = 0;
elseif isempty(modelTnseq)
    tnPresent = 0;
else
    tnPresent = 1;
end

% Check if median and standard deviation values are provided
if tnPresent == 1
    assert(nargin >= 6,...
        'Please provide median and standard deviation values');
    assert(~isempty(medianValue),...
        'Please provide median and standard deviation values');
    assert(~isempty(stdev),...
        'Please provide median and standard deviation values');
end

% Set default binning thresholds
if nargin < 7
    binThresh = {'3.5';'2.5';'1.5'};
elseif isempty(binThresh)
    binThresh = {'3.5';'2.5';'1.5'};
end

% Set solver
if nargin < 8
    global CBT_LP_SOLVER
    solver = CBT_LP_SOLVER;
    if isempty(solver)
        global CBTLPSOLVER
        solver = CBTLPSOLVER;
    end
end

% Set default method
if nargin < 9
    method = 1;
elseif isempty(method)
    method = 1;
end

%% Set solver

changeCobraSolver(solver);

if method == 2
    changeCobraSolver(solver, 'QP');
end

%% Randomize the non-protected gene list

nonProtected = nonProtected(randperm(length(nonProtected)));

%% Delete non-essential non-protected genes

% Determine the non-protected genes to delete
tempGenes = {};
deletableGenes = {};

if method == 1
    
    for n = 1:length(nonProtected)
        tempGenes = vertcat(tempGenes,nonProtected{n,1});
        [testModel,~,constrRxnNames,~] = deleteModelGenes(model,tempGenes);
        testModel = tncore_remove_reactions(testModel,constrRxnNames,false,false);
        testSolution = optimizeCbModel(testModel,'max');
        if testSolution.f > growthThresh
            deletableGenes = tempGenes;
        else
            tempGenes = deletableGenes;
        end
    end
    
elseif method == 2
    
    for n = 1:length(nonProtected)
        tempGenes = vertcat(tempGenes,nonProtected{n,1});
        [testModel] = deleteModelGenes(model,tempGenes);
        testSolution = MOMA(model,testModel,'max',false,true);
        if testSolution.f > growthThresh
            deletableGenes = tempGenes;
        else
            tempGenes = deletableGenes;
        end
    end

end

% Delete the deletable genes
[coreModel,~,constrRxnNames,~] = deleteModelGenes(model,deletableGenes);
coreModel = tncore_remove_reactions(coreModel,constrRxnNames);

%% Record growth rate

solution = optimizeCbModel(coreModel,'max');

%% Determine reactions present

reactions = cell(length(model.rxns),1);
rxnPresence = cell(length(model.rxns),1);

for n = 1:length(model.rxns)
    if strmatch(model.rxns{n},coreModel.rxns)
        reactions{n,1} = model.rxns{n,1};
        rxnPresence{n,1} = '1';
    else
        reactions{n,1} = [model.rxns{n,1} '_deleted'];
        rxnPresence{n,1} = '0';
    end
end

%% Replace gene names with Tnseq data

% Loop for if there is tnseq data
if tnPresent == 1

    genesTnseq = cell(length(model.genes),1);
    genePresence = cell(length(model.genes),1);

    for m = 1:length(coreModel.genes)
        IDs = num2cell(strmatch(model.genes{m},coreModel.genes(:,1),...
            'exact'));
        for o = 1:length(IDs)
            genesTnseq{IDs{o,1},1} = modelTnseq{m,1};
            genePresence{IDs{o,1},1} = '1';
        end
    end

end

% Loop for if there is no tnseq data
if tnPresent == 0

    genePresence = cell(length(model.genes),1);

    for m = 1:length(coreModel.genes)
        IDs = num2cell(strmatch(model.genes{m},coreModel.genes(:,1),...
            'exact'));
        for o = 1:length(IDs)
            genePresence{IDs{o,1},1} = '1';
        end
    end

end

% Replace empty with 0
isEmpty = cellfun('isempty',genePresence);
genePresence(isEmpty) = {'0'};

%% Count the genes in the core models in each essentiality category

if tnPresent == 1

    geneGrouping = cell(4,1);
    isEmpty = cellfun('isempty',geneGrouping);
    geneGrouping(isEmpty) = {0};

    for m = 1:length(coreModel.genes)
        if genesTnseq{m,1} < medianValue - (str2num(binThresh{1}) * stdev)
            geneGrouping{4,1} = geneGrouping{4,1} + 1;
        elseif genesTnseq{m,1} < medianValue - (str2num(binThresh{2}) * stdev)
            geneGrouping{3,1} = geneGrouping{3,1}+ 1;
        elseif genesTnseq{m,1} < medianValue - (str2num(binThresh{3}) * stdev)
            geneGrouping{2,1} = geneGrouping{2,1} + 1;
        elseif genesTnseq{m,1} >= medianValue - (str2num(binThresh{3}) * stdev)
        geneGrouping{1,1} = geneGrouping{1,1} + 1;
        end
    end

end

%% Things to fill in if tnseq is absent

if tnPresent == 0
    genesTnseq = {};
    geneGrouping = {};
end

