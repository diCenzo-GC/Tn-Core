function [coreModel] = tncore_reconstruct(model, genePresence, trim)

% Reconstruct a core model based on a gene presence/absence list.
% Compatible with producing a core model based on one column of the
% coreGeneVar array exported by the tncore_main function.
%
% USAGE
%   [coreModel] = tncore_reconstruct(model, genePresence, trim)
%
% INPUTS
%   model           The model from which the core model is derived. If
%                   reconstructing a core model generated by the  
%                   tncore_main function, use the reducedModel exported by 
%                   the tncore_main function.
%   genePresence    A binary array indicating the presence (1) or absence
%                   (0) of the gene in the core model.
%
% OPTIONAL INPUTS
%   trim            Indicates if non-essential reactions with no GPR or  
%                   with only an 'Unknwon' GPR should be removed from the 
%                   core model. Use 1 for trim, and 0 for do not trim 
%                   (Default = 1)
%
% OUTPUTS
%   coreModel       The core model
%
% AUTHORS
%   George diCenzo and Marco Fondi - 01/11/2017

%% Check variables

% Check there are enough inputs
assert(nargin >= 2,'This function requires at least two input');

% Set the default for trimming non-essential non-gene reactions
if nargin < 3
    trim = 1;
end

%% Produce the core model

% Identify genes to remove
if isnan(str2double(genePresence))
    isPresent = logical(cell2mat(genePresence));
else
    isPresent = logical(str2double(genePresence));
end
genes = model.genes(~isPresent);

% Delete the reactions constrained when the genes are deleted
[coreModel,~,constrRxnNames,~] = ...
    deleteModelGenes(model,genes);
coreModel = removeRxns(coreModel,constrRxnNames);

if trim == 1

    % Identify non-essential unknown or non-enzymatic reactions
    y = 0;
    reactionsToRemove = {};
    smallModelSolution = optimizeCbModel(coreModel,'max');

    for n = 1:length(coreModel.grRules)
        if strmatch(coreModel.grRules{n,1},'Unknown','exact')
            x = 1;
        elseif strmatch(coreModel.grRules{n,1},'Spontaneous')
            x = 1;
        elseif isempty(coreModel.grRules{n,1})
            x = 1;
        else
            x = 0;
        end
        if x == 1
            testModel = coreModel;
            testModel = removeRxns(testModel,coreModel.rxns{n,1});
            testModelSolution = optimizeCbModel(testModel,'max');
            if testModelSolution.f / smallModelSolution.f > 0.99
                y = y+1;
                reactionsToRemove{y,1} = coreModel.rxns{n,1};
            end
        end
    end

    % Remove non-essential unknown or non-enzymatic reactions
    if ~isempty(reactionsToRemove)
        for n = 1:length(reactionsToRemove)
            testModel = removeRxns(coreModel,reactionsToRemove{n,1});
            testModelSolution = optimizeCbModel(testModel,'max');
            if testModelSolution.f / smallModelSolution.f > 0.99
                coreModel = removeRxns(coreModel,reactionsToRemove{n,1});
            end
            smallModelSolution = optimizeCbModel(coreModel,'max');
        end
    end
        
end

% Delete genes from remaining grRules and fix model accordingly
coreModel = tncore_delete(coreModel);
