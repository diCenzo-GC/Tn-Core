function [modelNew] = tncore_fix(model)

%
% For any of the following fields that are absent in the input model, they
% will be added to the export model: rxnGeneMat, grRules, rules. At least
% one of the grRules or rules field must be present.
%
% USAGE
% [modelNew] = tncore_fix(model)
%
% Inputs
%   model     The model to be updated.
%
% Outputs
%   modelNew  The updated model.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 12/11/2018
%

%% Check a model is provided

assert(nargin == 1, 'Pleaes provide an input model');

%% Check which fields are present

% Look for rxnGeneMat
if isfield(model, 'rxnGeneMat')
    rxnGeneMatPres = true;
else
    rxnGeneMatPres = false;
end

% Look for rules
if isfield(model, 'rules')
    rulesPres = true;
else
    rulesPres = false;
end

% Look for grRules
if isfield(model, 'grRules')
    grRulesPres = true;
else
    grRulesPres = false;
end

% Error if none exist
if rulesPres == false
    assert(grRulesPres == 1, 'The input model must have at least one of rules, grRules, or rxnGeneMat');
end

% Make sure there is a genes field
assert(isfield(model, 'genes'), 'The input model must have a genes field');

%% If the grRules exists

if grRulesPres == true

    % Build rules if necessary
    if rulesPres == false
        
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
                    pos = strmatch(rulesSplit{m}, model.genes, 'exact');
                    newRule = ['x(' num2str(pos) ')'];
                    model.rules{n} = strrep(model.rules{n}, model.genes{pos}, newRule);
                end
            end
        end
        
    end
    
    % Build the rxnGeneMat if necessary
    if rxnGeneMatPres == false
        
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
    end

end

%% If the rules field exists

if rulesPres == true
    
    % Build rxnGeneMat if necessary
    if rxnGeneMatPres == false
        
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
    end
    
    % Build grRules if necessary
    if grRulesPres == false
    
        model.grRules = model.rules;
        for n = 1:length(model.grRules)
            model.grRules{n} = strrep(model.grRules{n}, '|', 'or');
            model.grRules{n} = strrep(model.grRules{n}, '&', 'and');
        end
        for n = 1:length(model.genes)
            for m = 1:length(model.grRules)
                rule = ['x(' num2str(n) ')'];
                model.grRules{m} = strrep(model.grRules{m}, rule, model.genes{n});
            end
        end
    
    end
    
end

%% Export the model

modelNew = model;
