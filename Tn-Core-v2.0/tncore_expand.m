function [expandedModel] = tncore_expand(coreModel, fullModel, deadends)

%
% Takes two input models and combines them. One model is a small model and
% the second is a large model. Combining them produces an expanded model -
% the small model with the additional reactions of the large model. GPRs
% are transferred with the reactions
%
% USAGE
%   [expandedModel] = tncore_expand(coreModel, fullModel, deadends)
%
% INPUTS
%   coreModel       The small model to serve as the core into which the 
%                   missing reactions from the large model will be added.
%   fullModel       The large model from which the additional reactions 
%                   will be taken.
%
% OPTIONAL INPUTS
%   deadends        Should reactions producing dead-ends be removed.
%                   1 if yes, 0 if no (Default = 1)
%
% OUTPUTS
%   expandedModel   The output model, consisting of the core model with the
%                   additional reactions of the full model, and optionally
%                   containing no dead-end metabolites.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 12/11/2018
%

%% Check inputs

% Ensure there are enough inputs
assert(nargin >= 2,'This function requires at least two inputs - the two models');

% Set default deadends
if nargin < 3
    deadends = 1;
end

%% Compare models by reaction abbreviations

% Set variables
reactionsCompared = {};

% Add headers to this array
reactionsCompared{1,1} = 'Reaction';
reactionsCompared{1,2} = 'SameName';
reactionsCompared{1,3} = 'SameFormula';
reactionsCompared{1,4} = 'SameReaction';

% Compare reaction names
for n = 1:length(fullModel.rxns);
    reactionsCompared{n+1,1} = fullModel.rxns{n};
    if strmatch(fullModel.rxns{n},coreModel.rxns,'exact');
        reactionsCompared{n+1,2} = '1';
    else
        reactionsCompared{n+1,2} = '0';
    end
end

%% Compare reaction equations

% Get reaction equation information for core model
coreEqs = printRxnFormula(coreModel);
metMatrix = {};
stoichMatrix = {};

for n = 1:length(coreEqs)
    [metList, stoichList] = parseRxnFormula(coreEqs{n});
    stoichList = num2cell(stoichList);
    metList = transpose(metList);
    stoichList = transpose(stoichList);
    info = horzcat(metList, stoichList);
    for m = 1:size(info, 1)
        metMatrix{m,n} = info{m,1};
        stoichMatrix{m,n} = info{m,2};
    end
end

% Compare each reaction of full model to the matrixes produced above
fullEqs = printRxnFormula(fullModel);

for n = 1:length(fullEqs)
    a = 0;
    b = 0;
    output = {};
    reactionsCompared{n+1,3} = '0';
    [metList, stoichList] = parseRxnFormula(fullEqs{n});
    stoichList = num2cell(stoichList);
    metList = transpose(metList);
    stoichList = transpose(stoichList);
    info = horzcat(metList, stoichList);
    % For each met in the rxn, find all associated rxns in core model
    for m = 1:size(info, 1)
        for o = 1:size(metMatrix,2)
            temp = metMatrix(:,o);
            empty = cellfun('isempty', temp);
            temp = temp(~empty);
            if strmatch(info{m,1}, temp, 'exact')
                b = a + 1;
                output{end+1,b} = o;
            end
        end
        a = b;
    end
    % If all mets are found in a rxn in the core model, find rxns with all
    if size(output,2) == size(info,1)
        if size(output,2) > 2
            output1 = output(:,1);
            isOutputEmpty = cellfun('isempty', output1);
            output1 = output1(~isOutputEmpty);
            output1 = coreModel.rxns(cell2mat(output1));
            output2 = output(:,2);
            isOutputEmpty = cellfun('isempty', output2);
            output2 = output2(~isOutputEmpty);
            output2 = coreModel.rxns(cell2mat(output2));
            reactions = intersect(output1, output2);
            for y = 3:size(output,2)
                output3 = output(:,y);
                isOutputEmpty = cellfun('isempty', output3);
                output3 = output3(~isOutputEmpty);
                output3 = coreModel.rxns(cell2mat(output3));
                reactions = intersect(reactions, output3);
            end
            reactions = findRxnIDs(coreModel, reactions);
        elseif size(output,2) == 2
            output1 = output(:,1);
            isOutputEmpty = cellfun('isempty', output1);
            output1 = output1(~isOutputEmpty);
            output1 = coreModel.rxns(cell2mat(output1));
            output2 = output(:,2);
            isOutputEmpty = cellfun('isempty', output2);
            output2 = output2(~isOutputEmpty);
            output2 = coreModel.rxns(cell2mat(output2));
            reactions = intersect(output1, output2);
            reactions = findRxnIDs(coreModel, reactions);
        else
            reactions = output(:,1);
            isOutputEmpty = cellfun('isempty', reactions);
            reactions = cell2mat(reactions(~isOutputEmpty));
        end
        % If at least one rxn has all mets, check stoichiometry the same
        if ~isempty(reactions)
            test = {};
            for y = 1:length(reactions)
                temp = metMatrix(:,reactions(y));
                temp = temp(~cellfun('isempty', temp));
                temp2 = stoichMatrix(:,reactions(y));
                temp2 = temp2(~cellfun('isempty', temp2));
                if size(info, 1) == length(temp)
                    for z = 1:size(info, 1)
                        pos = strmatch(info{z,1}, temp, 'exact');
                        stoichCoreRxn = temp2(pos);
                        if stoichCoreRxn{1,1} == info{z,2}
                            test{end+1,1} = 0;
                        elseif stoichCoreRxn{1,1} == -1 * info{z,2}
                            test{end+1,1} = 2;
                        else
                            test{end+1,1} = 1;
                        end
                    end
                    if sum(cell2mat(test)) == 0
                        reactionsCompared{n+1,3} = '1';
                    elseif sum(cell2mat(test)) / 2 == size(test,1)
                        reactionsCompared{n+1,3} = '1';
                   end
                end
            end
        end
    end
end

%% Compare the models

% Determine reactions unique to fullModel
for n = 1:(length(reactionsCompared)-1);
   if reactionsCompared{n+1,3} == '1';
      reactionsCompared{n+1,4} = '1';
   else
       reactionsCompared{n+1,4} = '0';
   end
end

% If a full model rxn has same name as core model rxn but a different
% formula, rename the rxn
for n = 1:(length(reactionsCompared)-1)
    x = 1;
    if reactionsCompared{n+1,3} == '0' && reactionsCompared{n+1,2} == '1'
        x = x + 1;
        number = ['_' num2str(x)];
        reactionsCompared{n+1,1} = strcat(reactionsCompared{n+1,1}, number);
    end
end

%% Add new metabolites

originalCoreMets = coreModel.mets;

% Get equations of reactions to add
newEquations = {};
fullModel.equation = printRxnFormula(fullModel);
coreModel.equation = printRxnFormula(coreModel);

for n = 1:(length(reactionsCompared)-1)
    if reactionsCompared{n+1,4} == '0'
        newEquations = vertcat(newEquations,fullModel.equation{n,1});
    end
end

% Find all metabolites in those reactions
allMets = {};

for n = 1:length(fullModel.mets)
    isPresent = cellfun(@isempty, strfind(newEquations,fullModel.mets{n}));
    if sum(isPresent) < length(isPresent)
        allMets = vertcat(allMets,fullModel.mets{n});
    end
end

% Identify only the missing metabolites
missingMets = {};

for n = 1:length(allMets)
    isPresent = cellfun(@isempty, strfind(coreModel.mets,allMets{n}));
    if sum(isPresent) == length(isPresent)
        missingMets = vertcat(missingMets,allMets{n});
    end
end

% Add missing metabolites to the model
coreModel.mets = vertcat(coreModel.mets,missingMets);

% Enlarge the S field
coreModel.S = full(coreModel.S);
coreModel.S = num2cell(coreModel.S);
Sfull = cell(length(coreModel.mets),length(coreModel.rxns));
Sfull(1:length(originalCoreMets),1:length(coreModel.rxns)) = coreModel.S;
isEmpty = cellfun('isempty',Sfull);
Sfull(isEmpty) = {0};
Sfull = cell2mat(Sfull);
Sfull = sparse(double(Sfull));
coreModel.S = Sfull;

% Add the info for the other met fields for the new metabolites
temp_fullModel = fullModel;
temp_fullModel.b = num2cell(temp_fullModel.b);

if (isfield(temp_fullModel,'metCharge') && ...
        ~isempty(temp_fullModel.metCharge))
    temp_fullModel.metCharge = num2cell(temp_fullModel.metCharge);
end

for n = (length(originalCoreMets)+1):length(coreModel.mets)
    metID = findMetIDs(fullModel,coreModel.mets{n});
    if (isfield(temp_fullModel,'metNames') && ...
            ~isempty(temp_fullModel.metNames)&& ...
            isfield(coreModel,'metNames') && ...
            ~isempty(coreModel.metNames))
        coreModel.metNames{n} = fullModel.metNames{metID};
    end
    if (isfield(temp_fullModel,'metChEBIID') && ...
            ~isempty(temp_fullModel.metChEBIID)&& ...
            isfield(coreModel,'metChEBIID') && ...
            ~isempty(coreModel.metChEBIID))
        coreModel.metChEBIID{n} = fullModel.metChEBIID{metID};
    end
    if (isfield(temp_fullModel,'metCharge') && ...
            ~isempty(temp_fullModel.metCharge)&& ...
            isfield(coreModel,'metCharge') && ...
            ~isempty(coreModel.metCharge))
        coreModel.metCharge = ...
            vertcat(coreModel.metCharge,temp_fullModel.metCharge{metID});
    end
    if (isfield(temp_fullModel,'metHMDB') && ...
            ~isempty(temp_fullModel.metHMDB) && ...
            isfield(coreModel,'metHMDB') && ...
            ~isempty(coreModel.metHMDB))
        coreModel.metHMDB{n} = fullModel.metHMDB{metID};
    end
    if (isfield(temp_fullModel,'metKEGGID') && ...
            ~isempty(temp_fullModel.metKEGGID) && ...
            isfield(coreModel,'metKEGGID') && ...
            ~isempty(coreModel.metKEGGID))
        coreModel.metKEGGID{n} = fullModel.metKEGGID{metID};
    end
    if (isfield(temp_fullModel,'metPubChemID') && ...
            ~isempty(temp_fullModel.metPubChemID) && ...
            isfield(coreModel,'metPubChemID') && ...
            ~isempty(coreModel.metPubChemID))
        coreModel.metPubChemID{n} = fullModel.metPubChemID{metID};
    end
    if (isfield(temp_fullModel,'metInChIString') && ...
            ~isempty(temp_fullModel.metInChIString) && ...
            isfield(coreModel,'metInChIString') && ...
            ~isempty(coreModel.metInChIString))
        coreModel.metInChIString{n} = fullModel.metInChIString{metID};
    end
    if (isfield(temp_fullModel,'b') && ~isempty(temp_fullModel.b) && ...
            isfield(coreModel,'b') && ~isempty(coreModel.b))
        coreModel.b = vertcat(coreModel.b,temp_fullModel.b{metID});
    end
end

%% Add all missing genes

missingGenes = setdiff(fullModel.genes,coreModel.genes);
coreModel.genes = vertcat(coreModel.genes,missingGenes);

%% Add the missing reactions

% Prepare a new file starting from core model
expModelDraft = coreModel;

% Change full model sections to array
temp_fullModel.rev = num2cell(temp_fullModel.rev);
temp_fullModel.lb = num2cell(temp_fullModel.lb);
temp_fullModel.ub = num2cell(temp_fullModel.ub);
temp_fullModel.c = num2cell(temp_fullModel.c);

% Make S matrix a full matrix
expModelDraft.S = full(expModelDraft.S);
expModelDraft.S = num2cell(expModelDraft.S);

% Change expModelDraft fields to cells
expModelDraft.rev = num2cell(expModelDraft.rev);
expModelDraft.lb = num2cell(expModelDraft.lb);
expModelDraft.ub = num2cell(expModelDraft.ub);
expModelDraft.c = num2cell(expModelDraft.c);

% Add the additional reactions
for n = 1:(length(reactionsCompared)-1);
    if reactionsCompared{n+1,4} == '0';        
        % Update reaction fields
        expModelDraft.rxns{end+1} = reactionsCompared{n+1,1};
        expModelDraft.rev{end+1} = temp_fullModel.rev{n,1};
        expModelDraft.c{end+1} = temp_fullModel.c{n,1};
        expModelDraft.lb{end+1} = temp_fullModel.lb{n,1};
        expModelDraft.ub{end+1} = temp_fullModel.ub{n,1};        
        expModelDraft.equation{end+1} = temp_fullModel.equation{n,1};
        if (isfield(temp_fullModel,'subSystems') && ...
                ~isempty(temp_fullModel.subSystems) && ...
                isfield(coreModel,'subSystems') && ...
                ~isempty(coreModel.subSystems))
            expModelDraft.subSystems{end+1} = fullModel.subSystems{n,1};
        end
        if (isfield(temp_fullModel,'confidenceScores') && ...
                ~isempty(temp_fullModel.confidenceScores) && ...
                isfield(coreModel,'confidenceScores') && ...
                ~isempty(coreModel.confidenceScores))
            expModelDraft.confidenceScores{end+1} = ...
                fullModel.confidenceScores{n,1};
        end
        if (isfield(temp_fullModel,'rxnReferences') && ...
                ~isempty(temp_fullModel.rxnReferences) && ...
                isfield(coreModel,'rxnReferences') && ...
                ~isempty(coreModel.rxnReferences))
            expModelDraft.rxnReferences{end+1} = fullModel.rxnReferences{n,1};
        end
        if (isfield(temp_fullModel,'rxnECNumbers') && ...
                ~isempty(temp_fullModel.rxnECNumbers) && ...
                isfield(coreModel,'rxnECNumbers') && ...
                ~isempty(coreModel.rxnECNumbers))
            expModelDraft.rxnECNumbers{end+1} = fullModel.rxnECNumbers{n,1};
        end
        if (isfield(temp_fullModel,'rxnNotes') && ...
                ~isempty(temp_fullModel.rxnNotes) && ...
                isfield(coreModel,'rxnNotes') && ...
                ~isempty(coreModel.rxnNotes))
            expModelDraft.rxnNotes{end+1} = fullModel.rxnNotes{n,1};
        end
        if (isfield(temp_fullModel,'rxnNames') && ...
                ~isempty(temp_fullModel.rxnNames) && ...
                isfield(coreModel,'rxnNames') && ...
                ~isempty(coreModel.rxnNames))
            expModelDraft.rxnNames{end+1} = fullModel.rxnNames{n,1};
        end
        
        % Update S matrix
        [metaboliteList,stoichCoeffList] = parseRxnFormula(fullModel.equation{n});
        metIDs = findMetIDs(expModelDraft,metaboliteList);
        metIDs = num2cell(metIDs);
        stoichCoeffList = num2cell(stoichCoeffList);
        lengthOfS = size(expModelDraft.S,2);
        for m = 1:length(metIDs)
            expModelDraft.S{metIDs{m},lengthOfS+1} = stoichCoeffList{m};
        end
        
        % Update gene fields
        expModelDraft.grRules{end+1} = temp_fullModel.grRules{n,1};
        expModelDraft.rules{end+1} = temp_fullModel.rules{n,1};
   end
end

% Convert S field to sparse
isEmpty = cellfun('isempty',expModelDraft.S);
expModelDraft.S(isEmpty) = {0};
expModelDraft.S = cell2mat(expModelDraft.S);
expModelDraft.S = sparse(double(expModelDraft.S));

% Restore format of different model sections
expModelDraft.rev = cell2mat(expModelDraft.rev);
expModelDraft.c = cell2mat(expModelDraft.c);
expModelDraft.lb = cell2mat(expModelDraft.lb);
expModelDraft.ub = cell2mat(expModelDraft.ub);

%% Fix gene fields

% Fixes rules section
rulesNew = expModelDraft.rules;

for n = 1:length(fullModel.genes);
    for m = size(expModelDraft.rxnGeneMat,1)+1:length(expModelDraft.rules);
        if strfind(expModelDraft.grRules{m},fullModel.genes{n});
            newNumber = strmatch(fullModel.genes{n},expModelDraft.genes,...
                'exact');
            old = ['x(' num2str(n) ')'];
            new = ['x(' num2str(newNumber) ')'];
            rulesNew{m,1} = strrep(rulesNew{m,1},old,new);
        end
    end
end

expModelDraft.rules = rulesNew;

% Update the rxnGeneMat field
rxnGeneMatFull = cell(length(expModelDraft.rxns),...
    length(expModelDraft.genes));

for n = 1:length(expModelDraft.rxns)
    if ~isempty(expModelDraft.rules{n})
        rulesTemp = strrep(expModelDraft.rules{n}, '&', '|');
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

isEmpty = cellfun('isempty',rxnGeneMatFull);
rxnGeneMatFull(isEmpty) = {0};
rxnGeneMatFull = cell2mat(rxnGeneMatFull);
rxnGeneMatFull = sparse(double(rxnGeneMatFull));
expModelDraft.rxnGeneMat = rxnGeneMatFull;
[expModelDraft] = tncore_remove(expModelDraft);

%% Find and remove dead ends
if deadends == 1
    
    % Find max upper and lower bounds
    upper = max(expModelDraft.ub);
    lower = min(expModelDraft.lb);
    upperBounds = num2cell(expModelDraft.ub);
    lowerBounds = num2cell(expModelDraft.lb);
    reversibility = num2cell(expModelDraft.rev);
    reactions = expModelDraft.rxns;
    expModelDraft.ub = num2cell(expModelDraft.ub);
    expModelDraft.lb = num2cell(expModelDraft.lb);
    expModelDraft.rev = num2cell(expModelDraft.rev);

    % Open all reactions while maintaining directionality
    for n = 1:length(expModelDraft.rev)
        expModelDraft.ub{n} = upper;
        if expModelDraft.rev{n} == 1
            expModelDraft.lb{n} = lower;
        end
    end
    
    expModelDraft.ub = cell2mat(expModelDraft.ub);
    expModelDraft.lb = cell2mat(expModelDraft.lb);
    expModelDraft.rev = cell2mat(expModelDraft.rev);
    
    % Find and remove all deadends
    test = 1;
    while test > 0
        deadEndMetabolites = detectDeadEnds(expModelDraft, true);
        deadEndMetabolites = num2cell(deadEndMetabolites);
        if isempty(deadEndMetabolites)
            test = 0;
        end
        if test == 1
            deadEndMetabolites2 = cell(length(deadEndMetabolites),1);
            for n = 1:length(deadEndMetabolites);
                deadEndMetabolites2{n,1} = ...
                    expModelDraft.mets{deadEndMetabolites{n,1},1};
            end
            [rxnList,~] = ...
                findRxnsFromMets(expModelDraft,deadEndMetabolites2);
            expModelDraft = tncore_remove_reactions(expModelDraft,rxnList);
            clear deadEndMetabolites
            clear deadEndMetabolites2
        end
    end

    % Restore original reaction bounds
    expModelDraft.ub = num2cell(expModelDraft.ub);
    expModelDraft.lb = num2cell(expModelDraft.lb);
    expModelDraft.rev = num2cell(expModelDraft.rev);

    for n = 1:length(reactions)
        rxnID = findRxnIDs(expModelDraft,reactions{n});
        if rxnID ~= 0
            expModelDraft.ub{rxnID} = upperBounds{n};
            expModelDraft.lb{rxnID} = lowerBounds{n};
            expModelDraft.rev{rxnID} = reversibility{n};
        end
    end
    
    expModelDraft.ub = cell2mat(expModelDraft.ub);
    expModelDraft.lb = cell2mat(expModelDraft.lb);
    expModelDraft.rev = cell2mat(expModelDraft.rev);

    % Remove genes no longer present in the model from the genes, rules, 
    % and rxnGeneMat fields
    expandedModel = tncore_remove(expModelDraft);

else

    expandedModel = expModelDraft;

end
