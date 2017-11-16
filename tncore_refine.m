function [refinedModel] = tncore_refine(model, tnseq, iters, ...
            growthThresh, coreGenes, deadends)

% Refine a genome-scale metabolic model based on Tn-seq data
%
% USAGE
%   [refinedModel] = tncore_refine(model, tnseq, iters, growthThresh, ...
%       coreGenes, deadends)
%
% INPUTS
%   model           The starting model.
%   tnseq           Tn-seq data, with the data in the first column and the
%                   gene names in the second column.
%
% OPTIONAL INPUTS
%   iters           Number of random core models to generate 
%                   (Default = 1000)
%   growthThresh:   The minimum allowable objective flux in the generated
%                   core models (Default = 10% the objective flux of the 
%                   input model)
%   coreGenes:      A list of all essential genes to be protected during 
%                   random model generation. If not provided, core genes
%                   are calculated from the Tn-seq data as those with a
%                   log-transformed value < the median - 3.5 * the standard
%                   deviation.
%   deadends        Should reactions producing dead-ends be removed.
%                   1 if yes, 0 if no (Default = 1)
%
% OUTPUTS
%   refinedModel    The original model but with the GPRs associated with
%                   core metabolic reactions updated based on the Tn-seq
%                   data, and optionally with dead-ends removed.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 01/11/2017

%% Check input variables

% Check there are enough variables
assert(nargin >= 2,'This function requires at least two input');

% Set default iterations
if nargin < 3
    iters = 1000;
elseif isempty(iters)
    iters = 1000;
end

% Set default growth threshold
if nargin < 4
    solutionOriginal = optimizeCbModel(model,'max');
    growthThresh = 0.1 * solutionOriginal.f;
elseif isempty(growthThresh)
    solutionOriginal = optimizeCbModel(model,'max');
    growthThresh = 0.1 * solutionOriginal.f;
end

% Set default coreGenes
if nargin < 5
    coreGenes = {};
elseif isempty(coreGenes)
    coreGenes = {};
end

% Set the default for deadends
if nargin < 6
    deadends = 0;
end

%% Set core genes if core genes were not provided by the user

% Keep non logtransformed data
tnseqNotLog = tnseq;

if isempty(coreGenes)
        
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
    temp = tnseq(isOutlier);
    nonOutliers = setdiff(cell2mat(tnseq(:,1)),cell2mat(temp));
    medianValue = median(nonOutliers);
    stdev = std(nonOutliers);

    % Set threshold for core genes
    coreThresh = medianValue - (3.5 * stdev);
    
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
    
    modelTnseq = sortrows(modelTnseq,1);
    
    % Get the core genes
    for n = 1:length(modelTnseq)
        if modelTnseq{n,1} < coreThresh
            coreGenes = vertcat(coreGenes,modelTnseq{n,2});
        end
    end

end

%% Produce core model

[~,~,~,coreModel] = tncore_main(model,iters,growthThresh,tnseqNotLog,...
    coreGenes);

%% Change OR to AND where appropriate

% Check which core genes are still in the model and do not have an effect
% when deleted
coreGenesRemaining = intersect(coreGenes,coreModel.genes);
[~,~,~,hasEffect] = singleGeneDeletion(coreModel,'FBA', coreGenesRemaining);
genesToTest = coreGenesRemaining(~hasEffect);

% For genes with no effect, if two are found only in the same reactions,
% and there is no AND statement in those reactions, replace OR with AND
if ~isempty(genesToTest)
    for n = 1:(length(genesToTest)-1)
        for m = (n+1):length(genesToTest)
            A = strmatch(coreModel.grRules,genesToTest{n});
            B = strmatch(coreModel.grRules,genesToTest{m});
            C = cell2mat(A) == cell2mat(B);
            D = sum(C) / length(C);
            if D == 1
                for o = 1:length(C)
                    if ~strfind(coreModel.grRules{C{o},1},'and')
                        coreModel.grRules{C{o},1} = ...
                            strrep(coreModel.grRules{C{o},1},' or ',' and ');
                    end
                end
            end
        end
    end
end
                    
%% Update GPRs in the full model

% Find GPRs to replace
[~,reactionsToUpdate] = findRxnsFromGenes(coreModel,coreGenesRemaining,0,1);
rxnIDcore = num2cell(findRxnIDs(coreModel,reactionsToUpdate));
rxnIDfull = num2cell(findRxnIDs(model,reactionsToUpdate));

% Replace the GPRs in the full model
for n = 1:length(rxnIDcore)
    model.grRules{rxnIDfull{n,1},1} = coreModel.grRules{rxnIDcore{n,1},1};
end

% Output refined model
refinedModel = model;
refinedModel = removeUnusedGenes(refinedModel);
refinedModel.rxnGeneMat = sparse(refinedModel.rxnGeneMat);

%% Trim the model to remove dead ends

if deadends == 1
   
    % Find max upper and lower bounds
    upper = max(refinedModel.ub);
    lower = min(refinedModel.lb);
    upperBounds = num2cell(refinedModel.ub);
    lowerBounds = num2cell(refinedModel.lb);
    reversibility = num2cell(refinedModel.rev);
    reactions = refinedModel.rxns;
    refinedModel.ub = num2cell(refinedModel.ub);
    refinedModel.lb = num2cell(refinedModel.lb);
    refinedModel.rev = num2cell(refinedModel.rev);
    
    % Open all reactions while maintaining directionality
    for n = 1:length(refinedModel.rev)
        refinedModel.ub{n} = upper;
        if refinedModel.rev{n} == 1
            refinedModel.lb{n} = lower;
        end
    end
    
    refinedModel.ub = cell2mat(refinedModel.ub);
    refinedModel.lb = cell2mat(refinedModel.lb);
    refinedModel.rev = cell2mat(refinedModel.rev);
    
    % Find and remove all deadends
    deadEndMetabolites = detectDeadEnds(refinedModel);

    if ~isempty(deadEndMetabolites)
        test = 1;
        while test > 0
            deadEndMetabolites = detectDeadEnds(refinedModel);
            if isempty(deadEndMetabolites)
                test = 0;
            else
                deadEndMetabolites3 = cell(length(deadEndMetabolites),1);
                for n = 1:length(deadEndMetabolites);
                    deadEndMetabolites2 = num2cell(deadEndMetabolites);
                    deadEndMetabolites3{n,1} = ...
                        refinedModel.mets{deadEndMetabolites2{n,1},1};
                end
                [deadEndReactions] = findRxnsFromMets(refinedModel,...
                    deadEndMetabolites3);
                refinedModel = removeRxns(refinedModel,deadEndReactions);
                clear deadEndMetabolites
                clear deadEndMetabolites2
            end
        end
        refinedModel = tncore_remove(refinedModel);
    end
    
    % Restore original reaction bounds
    refinedModel.ub = num2cell(refinedModel.ub);
    refinedModel.lb = num2cell(refinedModel.lb);
    refinedModel.rev = num2cell(refinedModel.rev);

    for n = 1:length(reactions)
        rxnID = findRxnIDs(refinedModel,reactions{n});
        if rxnID ~= 0
            refinedModel.ub{rxnID} = upperBounds{n};
            refinedModel.lb{rxnID} = lowerBounds{n};
            refinedModel.rev{rxnID} = reversibility{n};
        end
    end

    refinedModel.ub = cell2mat(refinedModel.ub);
    refinedModel.lb = cell2mat(refinedModel.lb);
    refinedModel.rev = cell2mat(refinedModel.rev);
    
end
