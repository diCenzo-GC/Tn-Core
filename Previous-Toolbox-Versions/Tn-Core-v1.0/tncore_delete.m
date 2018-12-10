function [modelNew] = tncore_delete(model)
% 
% Any gene marked with '_deleted' is removed from the gene list and the 
% grRules, and the rules and rxnGeneMat fields are updated accordingly.
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

%% Check input

assert(nargin == 1,'This function requires one input');

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

%% Remove deleted genes from grRules

% Remove genes from grRules
for n = 1:size(genesToDelete,1);
    for m = 1:length(model.grRules);
        model.grRules{m} = strrep(model.grRules{m},genesToDelete{n,1},'');
    end
end

% Remove unwanted statements from grRules
for m = 1:50
    for n = 1:length(model.grRules);
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
        model.grRules{n} = strrep(model.grRules{n},'()','');
        model.grRules{n} = strrep(model.grRules{n},'  ',' ');
        model.grRules{n} = strrep(model.grRules{n},'( ','(');
        model.grRules{n} = strrep(model.grRules{n},' )',')');
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
    end
end

%% Remove deleted genes from rules

% Remove genes from rules
for n = 1:length(genesToDelete);
    geneID = findGeneIDs(model,genesToDelete{n});
    geneToRemove = ['x(' num2str(geneID) ')'];
    for m = 1:length(model.grRules);
        model.rules{m} = strrep(model.rules{m},geneToRemove,'');
    end
end

% Remove unwanted statements from grRules
for m = 1:50
    for n = 1:length(model.rules);
        model.rules{n} = strrep(model.rules{n},'|  |',' | ');
        model.rules{n} = strrep(model.rules{n},' | | ',' | ');
        model.rules{n} = strrep(model.rules{n},'( | )','');
        model.rules{n} = strrep(model.rules{n},'(|)','');
        model.rules{n} = strrep(model.rules{n},'( | ','(');
        model.rules{n} = strrep(model.rules{n},' | )',')');
        model.rules{n} = strrep(model.rules{n},'(| ','(');
        model.rules{n} = strrep(model.rules{n},' |)',')');
        model.rules{n} = strrep(model.rules{n},'&  &',' & ');
        model.rules{n} = strrep(model.rules{n},' & & ',' & ');
        model.rules{n} = strrep(model.rules{n},'( & )','');
        model.rules{n} = strrep(model.rules{n},'(&)','');
        model.rules{n} = strrep(model.rules{n},'( & ','(');
        model.rules{n} = strrep(model.rules{n},' & )',')');
        model.rules{n} = strrep(model.rules{n},'(& ','(');
        model.rules{n} = strrep(model.rules{n},' &)',')');
        model.rules{n} = strrep(model.rules{n},'()','');
        model.rules{n} = strrep(model.rules{n},'  ',' ');
        model.rules{n} = strrep(model.rules{n},'( ','(');
        model.rules{n} = strrep(model.rules{n},' )',')');
        if strmatch(' | ', model.rules{n})
            model.rules{n} = extractAfter(model.rules{n}, 3);
        end
        if strmatch(' & ', model.rules{n})
            model.rules{n} = extractAfter(model.rules{n}, 3);
        end
        endStr = model.rules{n};
        if length(endStr) > 3
            endStr2 = endStr(end-2:end);
            if (strcmp(endStr2, ' | ') || strcmp(endStr2, ' & '))
                model.rules{n} = endStr(1:end-3);
            end
        end
    end
end

%% Remove genes no longer present in the model from the genes, rules, 
%  and rxnGeneMat fields

% Remove genes and produce final model
modelNew = tncore_remove(model);
