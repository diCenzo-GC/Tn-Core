function modelNew = tncore_remove_reactions(model, reactions)

%
% A function to remove reactions and ensure that reaction removal is done
% correctly when the number of reactions and metabolites are equal.
%
% USAGE
%   [modelNew] = tncore_delete(model)
%
% INPUT
%   model       The model from which the rections will be deleted.
%   reactions   An array of reactions to delete
%
% OUTPUT
%   modelNew    The model with the reactions deleted, and unused 
%               metabolites and genes also deleted.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 10/12/2018
%

%% Inputs

% Check number of inputs
assert(nargin == 2, 'This function requires two inputs')

% Rename the reaction input
rxnsToDelete = reactions;

%% If reactions and metabolites are different lengths

if length(model.rxns) ~= length(model.mets)
    modelNew = removeRxns(model, rxnsToDelete);
end

%% If reactions and metabolites are the same lengths

if length(model.rxns) == length(model.mets)
    
    % Count reactions
    rxnCount = length(model.rxns);
    
    % Find all fields
    fields = fieldnames(model);
    
    % Get rxn fields
    rxnFields = {};
    for n = 1:length(fields)
        if strmatch('rxn', fields{n})
            if length(model.(fields{n})) == rxnCount
                if strmatch('rxnGeneMat', fields{n})
                else
                    rxnFields = vertcat(rxnFields, fields{n});
                end
            end
        elseif strmatch('rev', fields{n}, 'exact')
            rxnFields = vertcat(rxnFields, fields{n});
        elseif strmatch('c', fields{n}, 'exact')
            rxnFields = vertcat(rxnFields, fields{n});
        elseif strmatch('lb', fields{n}, 'exact')
            rxnFields = vertcat(rxnFields, fields{n});
        elseif strmatch('ub', fields{n}, 'exact')
            rxnFields = vertcat(rxnFields, fields{n});
        elseif strmatch('rules', fields{n}, 'exact')
            rxnFields = vertcat(rxnFields, fields{n});
        elseif strmatch('grRules', fields{n}, 'exact')
            rxnFields = vertcat(rxnFields, fields{n});
        elseif strmatch('subSystems', fields{n}, 'exact')
            rxnFields = vertcat(rxnFields, fields{n});
        elseif strmatch('confidenceScores', fields{n}, 'exact')
            rxnFields = vertcat(rxnFields, fields{n});
        end
    end
    
    % Get met fields
    metFields = {};
    for n = 1:length(fields)
        if strmatch('met', fields{n})
            if length(model.(fields{n})) == length(model.mets)
                if strmatch('rxnGeneMat', fields{n})
                else
                    metFields = vertcat(metFields, fields{n});
                end
            end
        elseif strmatch('b', fields{n}, 'exact')
            metFields = vertcat(metFields, fields{n});
        end
    end
    
    % Get deleted reaction indexes
    rxnIDs = findRxnIDs(model, rxnsToDelete);
        
    % Delete rxns from rxn fields
    for n = 1:length(rxnFields)
        model.(rxnFields{n})(rxnIDs) = [];
    end
    
    % Delete rxns from rxnGeneMat and S fields
    model.rxnGeneMat(rxnIDs,:) = [];
    model.S(:,rxnIDs) = [];
        
    % Find removed mets
    metsToDelete = {};
    for n = 1:length(model.mets)
        if sum(abs(model.S(n,:))) == 0
            metsToDelete = vertcat(metsToDelete, model.mets{n});
        end
    end
    metIDs = findMetIDs(model, metsToDelete);
    
    % Delete mets from met fields
    for n = 1:length(metFields)
        model.(metFields{n})(metIDs) = [];
    end
    
    % Delete mets from the S field
    model.S(metIDs,:) = [];
    
    % Update genes
    modelNew = tncore_remove(model);
    
end
