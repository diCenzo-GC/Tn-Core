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
%

%% Remove no longer existant genes

% Remove no longer existant genes from gene field
genesNew = cell(1,1);
x = 0;
for n = 1:length(model.genes);
    for m = 1:length(model.grRules);
        if strfind(model.grRules{m},model.genes{n});
            x = x + 1;
            genesNew{x,1} = model.genes{n,1};
            break
        end
    end
end

% Update the rules field based on updated gene field
rulesNew = model.rules;
for n = 1:length(model.genes);
    for m = 1:length(model.rules);
        if strfind(model.grRules{m},model.genes{n});
            x = x + 1;
            newNumber = strmatch(model.genes{n},genesNew,'exact');
            old = ['x(' num2str(n) ')'];
            new = ['x(' num2str(newNumber) ')'];
            rulesNew{m,1} = strrep(rulesNew{m,1},old,new);
        end
    end
end

% Update the rxnGeneMat field based on updated gene field
rxnGeneMatNewFull = cell(length(model.rxns),length(genesNew));
for n = 1:length(genesNew)
    [results ListResults] = findRxnsFromGenes(model, genesNew{n},0,1);
    rxnIDs = num2cell(findRxnIDs(model,ListResults(:,1)));
    for m = 1:length(rxnIDs)
        rxnGeneMatNewFull{rxnIDs{m},n} = 1;
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
