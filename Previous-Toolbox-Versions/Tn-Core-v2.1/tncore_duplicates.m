function [modelNew] = tncore_duplicates(model)

%
% Removes duplicate genes (i.e., two genes with the same name) from a
% model. Updates the genes, rules, and rxnGeneMat fields. Also removes any
% genes not associated with a reaction.
%
% USAGE
% [modelNew] = tncore_duplicates(model)
%
% Inputs
%   model     The model from which to remove the duplicate genes.
%
% Outputs
%   modelNew  The model lacking the duplicated genes.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 10/12/2018
%

%% Check a model is provided

assert(nargin == 1, 'Pleaes provide an input model');

%% Remove no longer existant genes

% Remove duplicate genes from gene field
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
genesNew = unique(genesNew);

% Update rules
model.rules = model.grRules;
for n = 1:length(model.rxns)
    if ~isempty(model.rules{n})
        model.rules{n} = strrep(model.rules{n}, ' or ', ' | ');
        model.rules{n} = strrep(model.rules{n}, ' and ', ' & ');
        rulesTemp = model.rules{n};
        rulesTemp = strrep(rulesTemp, ' & ', ' | ');
        rulesSplit = strsplit(rulesTemp, ' | ');
        for m = 1:length(rulesSplit)
            rulesSplit{m} = strrep(rulesSplit{m}, '(', '');
            rulesSplit{m} = strrep(rulesSplit{m}, ')', '');
            rulesSplit{m} = strrep(rulesSplit{m}, ' ', '');
        end
        for m = 1:length(rulesSplit)
            pos = strmatch(rulesSplit{m}, genesNew, 'exact');
            newRule = ['x(' num2str(pos) ')'];
            model.rules{n} = strrep(model.rules{n}, genesNew{pos}, newRule);
        end
    end
end

% Update the rxnGeneMat field based on updated gene field
rxnGeneMatNewFull = cell(length(model.rxns),length(genesNew));
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

%% Prepare the final model

modelNew = model;
modelNew.genes = genesNew;
modelNew.rxnGeneMat = rxnGeneMatNew;

