function [reducedModel] = tncore_deadends(model)

%
% Iteratively remove all reactions producing deadend metabolites to ...
% produce a model free of deadend metabolites.
%
% USAGE
%   [reducedModel] = tncore_deadends(model)
%
% INPUTS
%   model           The model from which the deadends will be removed.
%
% OUTPUTS
%   reducedModel    The input model with deadends removed.
%
% AUTHORS
%   George diCenzo and Marco Fondi - 10/12/2018
%

%% Check a model is provided

assert(nargin == 1, 'Pleaes provide an input model');

%% Refine the model to remove dead ends

deadEndMetabolites = detectDeadEnds(model, true);

if ~isempty(deadEndMetabolites)
    test = 1;
    while test > 0
        deadEndMetabolites = detectDeadEnds(model, true);
        if isempty(deadEndMetabolites)
            test = 0;
        else
            deadEndMetabolites3 = cell(length(deadEndMetabolites),1);
            for n = 1:length(deadEndMetabolites);
                deadEndMetabolites2 = num2cell(deadEndMetabolites);
                deadEndMetabolites3{n,1} = ...
                    model.mets{deadEndMetabolites2{n,1},1};
            end
            [deadEndReactions] = findRxnsFromMets(model,deadEndMetabolites3);
            model = tncore_remove_reactions(model,deadEndReactions);
            if length(model.mets) < size(model.S, 1)
                model.S(deadEndMetabolites,:) = [];
                model.b(deadEndMetabolites,:) = [];
            end
            clear deadEndMetabolites
            clear deadEndMetabolites2
        end
    end
    reducedModel = tncore_remove(model);
else
    reducedModel = tncore_remove(model);
    fprintf('The input model has no deadends\n');
end
