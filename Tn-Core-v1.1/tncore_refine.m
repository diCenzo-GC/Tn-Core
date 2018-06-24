function [refinedModel] = tncore_refine(model, tnseq, epsilon, ...
    binThresh, essThresh, growthFrac, deadends)

%
% Refine a genome-scale metabolic model based on Tn-seq data
%
% USAGE
%   [refinedModel] = tncore_refine(model, tnseq, epsilon, binThresh, ...
%       essThresh, growthFrac, rnaseq, expressThresh, deadends)
%
% INPUTS
%   model           The starting model.
%   tnseq           Tn-seq data, with the data in the first column and the
%                   gene names in the second column.
%
% OPTIONAL INPUTS
%   epsilon         The epsilon value (the minimum flux rate that is
%                   considered non-zero) to be used during the running of
%                   FASTCC and FASTCORE. (Default = 1.1e-6)
%   binThresh       An array containing three values (number of standard
%                   deviations away from the mean of the log of the Tn-seq 
%                   data) to use for setting the limits of each bin. 
%                   (Default = {'3.5'; '2.5'; '1.5'}) 
%   essThresh       The threshold (in number of standard deviations from
%                   the mean of the log of the Tn-seq data) for a gene to 
%                   be considered essential. (Default = 3.5)
%   growthFrac      The minimum allowable objective flux in the generated
%                   core model, as a fraction of growth of the input model
%                   (Default = 0.5) 
%   deadends        Should reactions producing dead-ends be removed from 
%                   the refined model; 1 if yes, 0 if no (Default = 0)
%
% OUTPUTS
%   refinedModel    The original model but with the GPRs associated with
%                   core metabolic reactions updated based on the Tn-seq
%                   data, and optionally with dead-ends removed.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 01/11/2017
%   George diCenzo and Marco Fondi - updated - 24/06/2018
%

%% Check input variables

% Ensure there are enough inputs
assert(nargin >= 2,'This function requires atleast two inputs');

% Set default for epsilon
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

% Set the default essentiality threshold
if nargin < 5
    essThresh = 3.5;
elseif isempty(essThresh)
    essThresh = 3.5;
end

% Set growth threshold
if nargin < 6
    growthFrac = 0.5;
elseif isempty(growthFrac)
    growthFrac = 0.5;
end
solutionOriginal = optimizeCbModel(model,'max');
growthThresh = growthFrac * solutionOriginal.f;

% Set default for deadends
if nargin < 7
    deadends = 0;
elseif isempty(deadends)
    deadends = 0;
end

%% Produce core model

[coreModel] = tncore_core(model, tnseq, epsilon, binThresh, essThresh, ...
    growthFrac, [], [], [], 0);

%% Find the essential gene list

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

% Get the essential genes
essential = cell2mat(modelTnseq(:,1)) < essThresh;
coreGenes = modelTnseq(essential, 2);

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

% Update the rules
model.rules = model.grRules;

for n = 1:length(model.rules)
    if ~isempty(model.rules{n})
        model.rules{n} = strrep(model.rules{n}, ' or ', ' | ');
        model.rules{n} = strrep(model.rules{n}, ' and ', ' & ');
        for m = 1:length(model.genes)
            if strfind(model.rules{n}, model.genes{m})
                geneRule = ['x(' num2str(m) ')'];
                model.rules{n} = strrep(model.rules{n}, model.genes{m}, ...
                    geneRule);
            end
        end
    end
end

% Update the rxnGeneMat
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

isEmpty = cellfun('isempty',rxnGeneMatNewFull);
rxnGeneMatNewFull(isEmpty) = {0};
rxnGeneMatNewFull = cell2mat(rxnGeneMatNewFull);
rxnGeneMatNew = sparse(double(rxnGeneMatNewFull));
model.rxnGeneMat = rxnGeneMatNew;

% Fix the model
refinedModel = tncore_remove(model);

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
