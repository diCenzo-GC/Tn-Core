function [modelNew] = tncore_remove(model)

%
% Remove genes in the gene list that are no longer present a model, due to 
% deletion of reactions. Genes are removed from the genes, rules, and 
% rxnGeneMat fields.
%
% USAGE
% [modelNew] = tncore_remove(model)
%
% Inputs
%   model     The model from which to remove the genes.
%
% Outputs
%   modelNew  The model lacking the removed genes.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 01/11/2017
%   George diCenzo and Marco Fondi - updated - 10/12/2018
%   George diCenzo and Marco Fondi - updated - 25/03/2019
%   George diCenzo - updated - 3/11/2022
%

%% Remove no longer existant genes

% Remove no longer existant genes from gene field
genesNew = cell(1,1);
x = 0;
allGrRules = {};
for n = 1:length(model.grRules)
    rxnGrRules = strrep(model.grRules{n}, '(', '');
    rxnGrRules = strrep(rxnGrRules, ')', '');
    rxnGrRules = strrep(rxnGrRules, ' and ', ' or ');
    rxnGrRules = transpose(strsplit(rxnGrRules, ' or '));
    rxnGrRules = strrep(rxnGrRules, ' ', '');
    allGrRules = vertcat(allGrRules, rxnGrRules);
end
allGrRules = unique(allGrRules);
genesNew = intersect(model.genes, allGrRules);

% Update the rules field based on updated gene field
rulesNew = model.rules;
for n = 1:length(model.genes);
    for m = 1:length(model.rules);
        if strfind(model.grRules{m},model.genes{n});
            x = x + 1;
            newNumber = strmatch(model.genes{n},genesNew,'exact');
            old = ['x(' num2str(n) ')'];
            new = ['x(' cell2mat(genesNew(newNumber)) ')'];
            rulesNew{m,1} = strrep(rulesNew{m,1},old,new);
        end
    end
end
for n = 1:length(model.genes);
    for m = 1:length(model.rules);
        if strfind(model.grRules{m},model.genes{n});
            x = x + 1;
            newNumber = strmatch(model.genes{n},genesNew,'exact');
            old = ['x(' cell2mat(genesNew(newNumber)) ')'];
            new = ['x(' num2str(newNumber) ')'];
            rulesNew{m,1} = strrep(rulesNew{m,1},old,new);
        end
    end
end

% Update the rxnGeneMat field based on updated gene field
rxnGeneMatNewFull = cell(length(model.rxns),length(genesNew));

for n = 1:length(model.rxns)
    if ~isempty(rulesNew{n})
        rulesTemp = strrep(rulesNew{n}, '&', '|');
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

%% Prepare the final model

modelNew = model;
modelNew.genes = genesNew;
modelNew.rules = rulesNew;
modelNew.rxnGeneMat = rxnGeneMatNew;

