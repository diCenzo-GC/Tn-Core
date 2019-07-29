function [model, tnseq, rnaseq] = tncore_import()

%
% An optional function to import the data to be used with Tn-Core, and to 
% update the model to be in the correct format. Data to be imported must be
% present in the current working directory. Files to import should be named
% as indicated below. Only the files for model import are required; the
% others are optional.
%
% USAGE
%   [model, tnseq, rnaseq] = tncore_import()
% 
% INPUT FILES FOR MODEL IMPORT
%   inputModel.xml      A SBML formatted metabolic model
%   ExRxns.txt          A tab-delimited file indicating the exchange
%                       reactions to be on (first column), and the lower 
%                       boundary (second column). Exchange reactions not
%                       listed will have a lower bound of 0, unless it is
%                       essential
%   ObjectiveRxn.txt    A tab-delimited file indicating the reaction to be
%                       set as the objective function (first column) and
%                       the objective value (second column)
%
% INPUT FILES FOR IMPORTING TN-SEQ DATA
%   TnSeqData.txt       A tab-delimited file with the gene essentiality
%                       index in the first column, and the gene name in the
%                       second column
%
% INPUT FILES FOR IMPORTING TNA-SEQ DATA
%   RnaSeqData.txt      A tab-delimited file with the gene expression data
%                       in the first column, and the gene name in the
%                       second column
%
% OUTPUTS
%   model               A COBRA formatted metabolic model with the
%                       objective function set, the exchange bounds set,
%                       and with the necessary fields to work with Tn-Core
%   tnseq               The Tn-seq data formatted for use with Tn-Core
%   rnaseq              The RNA-seq data formatted for use with Tn-Core
%
% AUTHORS
%   George diCenzo and Marco Fondi - 06/04/2018
%   George diCenzo and Marco Fondi - updated - 25/03/2019
%

%% Check inputs

% Is there a model file
if exist('inputModel.xml') == 2
    inModelPres = true;
else
    inModelPres = false;
end
assert(inModelPres == true, 'Please prepare an inputModel.xml file with the exchange reactions to set, and place in the current directory');

% Is there a ExRxns file
if exist('ExRxns.txt') == 2
    ExRxnsPres = true;
else
    ExRxnsPres = false;
end
assert(ExRxnsPres == true, 'Please prepare a ExRxns.txt file with the exchange reactions to set, and place in the current directory');

% Is there a ObjectiveRxn file
if exist('ObjectiveRxn.txt') == 2
    ObjRxnPres = true;
else
    ObjRxnPres = false;
end
assert(ObjRxnPres == true, 'Please prepare a ObjectiveRxn.txt file with the reactions to set at the objective, and place in the current directory');

%% Inport and prepare the metabolic model

% Load the model
model = readCbModel('inputModel.xml');

% Prepare a grRules field if it is missing
if ~isfield(model, 'grRules')
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

% Prepare a rxnGeneMat field if it is missing
if ~isfield(model, 'rxnGeneMat')
    rxnGeneMatNewFull = cell(length(model.rxns),length(model.genes));
    for n = 1:length(model.rxns)
        for m = 1:length(model.genes)
            if strfind(model.grRules{n}, model.genes{m})
                rxnGeneMatNewFull{n,m} = 1;
            end
        end
    end
    isEmpty = cellfun('isempty',rxnGeneMatNewFull);
    rxnGeneMatNewFull(isEmpty) = {0};
    rxnGeneMatNewFull = cell2mat(rxnGeneMatNewFull);
    rxnGeneMatNew = sparse(double(rxnGeneMatNewFull));
    model.rxnGeneMat = rxnGeneMatNew;
end

% Get the exchange reactions
exRxns = table2cell(readtable('ExRxns.txt', 'ReadVariableNames', false, ...
    'ReadRowNames', false, 'Delimiter', '\t'));
allExRxns = model.rxns(findExcRxns(model));
if ~isempty(exRxns)
    allExRxns = setdiff(allExRxns, exRxns(:,1));
end

% Set the objective function
objfun = table2cell(readtable('ObjectiveRxn.txt', 'ReadVariableNames', false, ...
    'ReadRowNames', false, 'Delimiter', '\t'));
model = changeObjective(model, objfun{1,1}, objfun{1,2});

% Set the exchange reactions
model = changeRxnBounds(model, allExRxns, 0, 'l');
if ~isempty(exRxns)
    model = changeRxnBounds(model, exRxns(:,1), cell2mat(exRxns(:,2)), 'l');
end

% Check if model grows, and if not, find missing Ex rxns
sol = optimizeCbModel(model);
if round(sol.f, 6) == 0
    model = changeRxnBounds(model, allExRxns, -10, 'l');
    if ~isempty(exRxns)
        model = changeRxnBounds(model, exRxns(:,1), cell2mat(exRxns(:,2)), 'l');
    end
    for n = 1:length(allExRxns)
        model = changeRxnBounds(model, allExRxns{n}, 0, 'l');
        testSol = optimizeCbModel(model);
        if round(testSol.f, 6) == 0
            model = changeRxnBounds(model, allExRxns{n}, -10, 'l');
        end
    end
end

%% Get the Tn-seq and RNA-seq data

% Get the Tn-seq data
if exist('TnSeqData.txt') == 2
    tnseq = table2cell(readtable('TnSeqData.txt', 'ReadVariableNames', false, ...
    'ReadRowNames', false, 'Delimiter', '\t'));
else
    tnseq = {};
end

% Get the RNA-seq data
if exist('RnaSeqData.txt') == 2
    rnaseq = table2cell(readtable('RnaSeqData.txt', 'ReadVariableNames', false, ...
    'ReadRowNames', false, 'Delimiter', '\t'));
else
    rnaseq = {};
end
