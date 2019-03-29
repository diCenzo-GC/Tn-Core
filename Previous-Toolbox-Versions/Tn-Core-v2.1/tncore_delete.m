function [modelNew] = tncore_delete(model)

% 
% Any gene marked with '_deleted' is removed from the gene list and the 
% grRules, and the rules and rxnGeneMat fields are updated accordingly.
% Designed for use to remove deleted genes after removing reactions 
% constrained with deleteModelGenes.
%
% USAGE
%   [modelNew] = tncore_delete(model)
%
% INPUT
%   model       The model from which the genes will be deleted.
%
% OUTPUT
%   modelNew    The model with the genes deleted.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 01/11/2017
%   George diCenzo and Marco Fondi - updated - 12/04/2018
%   George diCenzo and Marco Fondi - updated - 10/12/2018
%

%% Check input

assert(nargin == 1,'This function requires one input');

% Check if there is at least one gene to delete
x = 0;
for n = 1:length(model.genes);
    if strfind(model.genes{n,1},'_deleted')
        x = x + 1;
    end
end
assert(x >= 1,'No genes to be deleted were found');


%% Make list of genes to remove and remove _deleted from gene names

genesToDelete = cell(1,1);
x = 0;

for n = 1:length(model.genes);
    if strfind(model.genes{n,1},'_deleted')
        x = x + 1;
        model.genes{n,1} = strrep(model.genes{n,1},'_deleted','');
        genesToDelete{x,1} = model.genes{n,1};
        genesToDelete{x,2} = n;
    end
end

%% Ensure the grRules do not contain _deleted

for n = 1:length(model.grRules)
    model.grRules{n} = strrep(model.grRules{n}, '_deleted', '');
end

%% Remove deleted genes from grRules

% Save the original model
origModel = model;

% Remove genes from grRules
for n = 1:size(genesToDelete,1);
    for m = 1:length(model.grRules);
        if strfind(model.grRules{m}, genesToDelete{n,1})
            if strfind(model.grRules{m}, ' and ')
                if strfind(model.grRules{m}, ' or ')
                    grRulesTemp = strrep(origModel.grRules{m}, ' and ', ' or ');
                    genesOfRxn = transpose(strsplit(grRulesTemp, ' or '));
                    for i = 1:size(genesOfRxn, 1)
                        genesOfRxn{i,1} = strrep(genesOfRxn{i,1}, '(', '');
                        genesOfRxn{i,1} = strrep(genesOfRxn{i,1}, ')', '');
                        genesOfRxn{i,1} = strrep(genesOfRxn{i,1}, ' ', '');
                    end
                    genesOfRxn = genesOfRxn(~cellfun('isempty', genesOfRxn));
                    genesTemp = {};
                    for i = 1:size(genesOfRxn, 1)
                        [~,~,constrRxns] = deleteModelGenes(model, genesOfRxn{i,1});
                        if isempty(constrRxns)
                            genesTemp = vertcat(genesTemp, genesOfRxn{i,1});
                        else
                            matchedConstrRxns = strmatch(model.rxns{n}, constrRxns, 'exact');
                            if isempty(matchedConstrRxns)
                                genesTemp = vertcat(genesTemp, genesOfRxn{i,1});
                            end
                        end
                    end
                    genesEss = setdiff(genesOfRxn, genesTemp);
                    genesTemp = setdiff(genesTemp, genesToDelete{n});
                    genesTempPerm = {};
                    for i = 1:(5*length(genesTemp))
                        genesTempPerm(i,:) = transpose(genesTemp(randperm(length(genesTemp))));
                    end
                    essentiality = {};
                    for i = 1:size(genesTempPerm, 1)
                        genesTempDelete = {};
                        genesTempDelete{1,1} = genesToDelete{n};
                        for j = 1:size(genesTempPerm, 2)
                            genesTempDelete = vertcat(genesTempDelete, genesTempPerm{i,j});
                            [~,~,constrRxns] = deleteModelGenes(model, genesTempDelete);
                            if strmatch(model.rxns{m}, constrRxns, 'exact')
                                essentiality{i,j} = 1;
                                genesTempDelete = genesTempDelete(1:end-1);
                            else
                                essentiality{i,j} = 0;
                            end
                        end
                    end
                    genesTempPerm = transpose(genesTempPerm);
                    essentiality = transpose(essentiality);
                    for i = 1:size(genesTemp, 1)
                        total = {};
                        for j = 1:size(genesTempPerm, 2)
                            pos = strmatch(genesTemp{i}, genesTempPerm(:,j), 'exact');
                            total = vertcat(total, essentiality(pos, j));
                        end
                        if sum(cell2mat(total)) >= 1
                            genesEss = vertcat(genesEss, genesTemp{i});
                        end
                    end
                    genesNonEss = setdiff(genesOfRxn, genesEss);
                    for i = 1:length(genesNonEss)
                        model.grRules{m} = strrep(model.grRules{m},genesNonEss{i,1},'');
                    end
                else
                    model.grRules{m} = strrep(model.grRules{m},genesToDelete{n,1},'');
                end
            else
                model.grRules{m} = strrep(model.grRules{m},genesToDelete{n,1},'');
            end
        end
    end
end

% Remove unwanted statements from grRules
for m = 1:50
    for n = 1:length(model.grRules);
        model.grRules{n} = strrep(model.grRules{n},'()','');
        model.grRules{n} = strrep(model.grRules{n},'  ',' ');
        model.grRules{n} = strrep(model.grRules{n},'( ','(');
        model.grRules{n} = strrep(model.grRules{n},' )',')');
        model.grRules{n} = strrep(model.grRules{n},'or  or',' or ');
        model.grRules{n} = strrep(model.grRules{n},' or or ',' or ');
        model.grRules{n} = strrep(model.grRules{n},'( or )','');
        model.grRules{n} = strrep(model.grRules{n},'(or)','');
        model.grRules{n} = strrep(model.grRules{n},'( or ','(');
        model.grRules{n} = strrep(model.grRules{n},' or )',')');
        model.grRules{n} = strrep(model.grRules{n},'(or ','(');
        model.grRules{n} = strrep(model.grRules{n},' or)',')');
        model.grRules{n} = strrep(model.grRules{n},'and  and',' and ');
        model.grRules{n} = strrep(model.grRules{n},' and and ',' and ');
        model.grRules{n} = strrep(model.grRules{n},'( and )','');
        model.grRules{n} = strrep(model.grRules{n},'(and)','');
        model.grRules{n} = strrep(model.grRules{n},'( and ','(');
        model.grRules{n} = strrep(model.grRules{n},' and )',')');
        model.grRules{n} = strrep(model.grRules{n},'(and ','(');
        model.grRules{n} = strrep(model.grRules{n},' and)',')');
        if strmatch(' or ', model.grRules{n})
            model.grRules{n} = extractAfter(model.grRules{n}, 4);
        end
        if strmatch(' and ', model.grRules{n})
            model.grRules{n} = extractAfter(model.grRules{n}, 5);
        end
        endStr = model.grRules{n};
        if length(endStr) > 4
            endStr2 = endStr(end-3:end);
            if strcmp(endStr2, ' or ')
                model.grRules{n} = endStr(1:end-4);
            end
            endStr2 = endStr(end-4:end);
            if strcmp(endStr2, ' and ')
                model.grRules{n} = endStr(1:end-5);
            end
        end
        if isstring(model.grRules{n})
            model.grRules{n} = str2mat(model.grRules{n});
        end
    end
end

%% Remove deleted genes from rules

% Update the rules based on the grRules
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
            rulesSplit{m} = strrep(rulesSplit{m}, 'x', '');
            rulesSplit{m} = strrep(rulesSplit{m}, ' ', '');
        end
        for m = 1:length(rulesSplit)
            pos = findGeneIDs(model, rulesSplit{m});
            newRule = ['x(' num2str(pos) ')'];
            model.rules{n} = strrep(model.rules{n}, model.genes{pos}, newRule);
        end
    end
end

%% Remove genes no longer present in the model from the genes, rules, 
%  and rxnGeneMat fields

% Remove genes and produce final model
modelNew = tncore_remove(model);
