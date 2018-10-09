function [refinedModel, report] = tncore_refine(model, tnseq, epsilon, ...
    binThresh, essThresh, growthFrac, deadends)

%
% Refine a genome-scale metabolic model based on Tn-seq data
%
% USAGE
%   [refinedModel] = tncore_refine(model, tnseq, epsilon, binThresh, ...
%       essThresh, growthFrac, deadends)
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
%   report          A table displaying the changes made in preparing the
%                   refined model. Also includes a list of additional
%                   potential changes that are not included in the refined
%                   model.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 01/11/2017
%   George diCenzo and Marco Fondi - updated - 24/06/2018
%   George diCenzo and Marco Fondi - updated - 09/10/2018

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

%% Store the original model

origModel = model;

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
essentialityThreshold = medianValue - (stdev * essThresh);
essential = cell2mat(modelTnseq(:,1)) < essentialityThreshold;
coreGenes = modelTnseq(essential, 2);

%% Fix incorrectly modified GPRs

% Check which core genes are still in the model 
coreGenesRemaining = intersect(coreGenes, coreModel.genes);

% Find GPRs with essential core gene
testEss = singleGeneDeletion(coreModel, 'FBA', coreGenesRemaining);
testEss = round(testEss, 2);
coreGenesRemainingTemp = {};
for n = 1:length(testEss)
    if testEss(n) == 0
        coreGenesRemainingTemp = vertcat(coreGenesRemainingTemp, coreGenesRemaining{n});
    end
end
coreGenesRemaining = coreGenesRemainingTemp;

% Find genes with rxns with essential core gene that are inappropriately
% essential
[~, essGeneRxns] = findRxnsFromGenes(coreModel, coreGenesRemaining, [], 1);
essGeneRxns = unique(essGeneRxns(:,1));
essGeneRxnGenes = unique(findGenesFromRxns(coreModel, essGeneRxns));
otherGenes = setdiff(essGeneRxnGenes, coreGenesRemaining);
otherGenesFitness = singleGeneDeletion(coreModel, 'FBA', otherGenes);
otherGenesFitness = round(otherGenesFitness, 2);
otherGenesTemp = {};
for n = 1:length(otherGenesFitness)
    if otherGenesFitness(n) == 0
        otherGenesTemp = vertcat(otherGenesTemp, otherGenes{n});
    end
end
otherGenes = otherGenesTemp;
otherGenesFitnessOrig = singleGeneDeletion(model, 'FBA', otherGenes);
otherGenesFitnessOrig = round(otherGenesFitnessOrig, 2);
otherGenesTemp = {};
for n = 1:length(otherGenesFitnessOrig)
    if otherGenesFitnessOrig(n) > 0
        otherGenesTemp = vertcat(otherGenesTemp, otherGenes{n});
    end
end
otherGenes = otherGenesTemp;

% Find reactions associated with otherGenes
[~, otherGenesRxns] = findRxnsFromGenes(coreModel, otherGenes, [], 1);
otherGenesRxns = unique(otherGenesRxns(:,1));
otherGenesRxnsInv = setdiff(coreModel.rxns, otherGenesRxns);

% Make expanded temporary core model
reactionsToDelete = setdiff(model.rxns, coreModel.rxns);
tempCoreModel = removeRxns(model, reactionsToDelete);
tempCoreModel = tncore_remove(tempCoreModel);

% Change GPRs in temporary core model
otherGenesRxnIDsCore = findRxnIDs(coreModel, otherGenesRxnsInv);
otherGenesRxnIDsFull = findRxnIDs(tempCoreModel, otherGenesRxnsInv);

% Replace the GPRs in the core model
for n = 1:length(otherGenesRxnIDsCore)
    tempCoreModel.grRules{otherGenesRxnIDsCore(n),1} = ...
        coreModel.grRules{otherGenesRxnIDsFull(n),1};
end

% Update the rules
tempCoreModel.rules = tempCoreModel.grRules;

for n = 1:length(tempCoreModel.rules)
    if ~isempty(tempCoreModel.rules{n})
        tempCoreModel.rules{n} = strrep(tempCoreModel.rules{n}, ' or ', ' | ');
        tempCoreModel.rules{n} = strrep(tempCoreModel.rules{n}, ' and ', ' & ');
        for m = 1:length(tempCoreModel.genes)
            if strfind(tempCoreModel.rules{n}, tempCoreModel.genes{m})
                geneRule = ['x(' num2str(m) ')'];
                tempCoreModel.rules{n} = strrep(tempCoreModel.rules{n}, tempCoreModel.genes{m}, ...
                    geneRule);
            end
        end
    end
end

% Update the rxnGeneMat
rxnGeneMatNewFull = cell(length(tempCoreModel.rxns),length(tempCoreModel.genes));

for n = 1:length(tempCoreModel.rxns)
    if ~isempty(tempCoreModel.rules{n})
        rulesTemp = strrep(tempCoreModel.rules{n}, '&', '|');
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
tempCoreModel.rxnGeneMat = rxnGeneMatNew;

% Fix the model
coreModelTemp2 = tncore_remove(tempCoreModel);
coreModel = coreModelTemp2;

%% Update GPRs in the full model

% Check which core genes are still in the model and do not have an effect
% when deleted
coreGenesRemaining = intersect(coreGenes,coreModel.genes);

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

%% Prepare the first part of the report

% Prepare an array to hold the report
report = cell(1,5);
report{1,1} = 'Reaction';
report{1,2} = 'Type_of_change';
report{1,3} = 'Original_GPR';
report{1,4} = 'New_GPR';
report{1,5} = 'Is_change_made_in_refined_model';

% Identify and record differences in GPRs
for n = 1:length(refinedModel.rxns)
    if strmatch(refinedModel.grRules{n}, origModel.grRules{n}, 'exact')
    else
        report{end+1,1} = refinedModel.rxns{n};
        report{end,2} = 'GPR_modification';
        report{end,3} = origModel.grRules{n};
        report{end,4} = refinedModel.grRules{n};
        report{end,5} = 'Yes';
    end
end

%% Change OR to AND where appropriate
% Don't make the change in the actual model, but note the change in the
% report

% Set a temporary model
tempCoreModel = coreModel;

% Check which core genes are still in the model and do not have an effect
% when deleted
coreGenesRemaining = intersect(coreGenes,tempCoreModel.genes);
[~,~,~,hasEffect] = singleGeneDeletion(tempCoreModel,'FBA', coreGenesRemaining);
genesToTest = coreGenesRemaining(~hasEffect);

% For genes with no effect, if two are found only in the same reactions,
% and there is no AND statement in those reactions, replace OR with AND
if ~isempty(genesToTest)
    for n = 1:(length(genesToTest)-1)
        for m = (n+1):length(genesToTest)
            [~, A] = findRxnsFromGenes(tempCoreModel, genesToTest{n}, [], 1);
            [~, B] = findRxnsFromGenes(tempCoreModel, genesToTest{m}, [], 1);
            A = A(:,1);
            B = B(:,1);
            if length(A) == length(B)
                C = cell2mat(A) == cell2mat(B);
                C = C(:,1);
                D = sum(C) / length(C);
                E = findRxnIDs(coreModel, A);
                if D == 1
                    for o = 1:length(E)
                        if strfind(tempCoreModel.grRules{E(o),1},'and')
                        else
                            tempCoreModel.grRules{E(o),1} = ...
                                strrep(tempCoreModel.grRules{E(o),1},' or ',' and ');
                        end
                    end
                end
            end
        end
    end
end

%% Prepare the second part of the report

% Identify and record potential 'or' to 'and' transitions
for n = 1:length(coreModel.rxns)
    if strmatch(coreModel.grRules{n}, tempCoreModel.grRules{n}, 'exact')
    else
        report{end+1,1} = coreModel.rxns{n};
        report{end,2} = 'OR_to_AND_transition';
        report{end,3} = coreModel.grRules{n};
        report{end,4} = tempCoreModel.grRules{n};
        report{end,5} = 'No';
    end
end

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

%% Prepare the third part of the report

% Find and record reactions removed during dead-end removal
for n = 1:length(origModel)
    if strmatch(origModel.rxns{n}, refinedModel.rxns(:), 'exact')
    else
        report{end+1,1} = origModel.rxns{n};
        report{end,2} = 'Deadend_reaction_deletion';
        report{end,3} = '-';
        report{end,4} = '-';
        report{end,5} = 'Yes';
    end
end

% Convert report to a table
report = cell2table(report(2:end,:), 'VariableNames', report(1,:));

