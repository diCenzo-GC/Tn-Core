function [coocurMatrix, coocurMatrixAdj, growthMatrix, ...
            varPresenceNew, uniqueModels] = tncore_matrix(model, ...
            variable, growth, varPresence, type, minPresence)

%
% Produces matrixes (cooccurrence matrix, growth matrix) summarizing the
% variability in the core models produced by the tncore_main function.
%
% USAGE
%   [coocurMatrix, coocurMatrixAdj, growthMatrix, varPresenceLabel, ...
%       uniqueModels] = tncore_matrix(model, variable, growth, ...
%       varPresence, type, minPresence)
%
% INPUTS
%   model               The reducedModel exported from tncore_main.
%   variable            The coreGeneVar or coreRxnVar arrays exported from
%                       tncore_redundancy.
%   growth              The coreGrowthVar array exported from tncore_main.
%   varPresence         The genePresence or rxnPresence arrays exported
%                       from tncore_redundancy.
%
% OPTIONAL INPUTS
%   type                Type of matrixes to produce.  Use 1 for analysis of
%                       the genes and use 2 for an analysis of the
%                       reactions (Default = 1)
%   minPresence         The minimum number of models a gene/reaction must
%                       be present in to be included output matrixes
%                       (Default = 1)
%
% OUTPUTS
%   coocurMatrix        An array indicating the total number of times each
%                       pair of genes or reactions occured in the same
%                       model.
%   coocurMatrixAdj     An array indicating the ratio of models that have
%                       gene/rxn X that also have gene/rxn Y
%   growthMatrix        For each gene, indicates the average growth rate of
%                       core models containing the gene, and the total
%                       number of core models containing the gene
%   varPresenceNew      The varPresence input variable containing row and
%                       column labels, and also with genes/rxns not 
%                       occuring in any model removed.
%   unique              The number of unique genes returned, based on
%                       comparison of the genes if input is type 1, or on a
%                       comparison of the reactions if input is type 2
%
% AUTHORS
%   George diCenzo and Marco Fondi - 01/11/2017
%

%% Deal with input variables

% Ensure there are enough inputs
assert(nargin >= 4,'This function requires atleast four inputs');

% Determine default variable type
if nargin < 5
    type = 1;
elseif isempty(type)
    type = 1;
end

% Set default minPresence
if nargin < 6
    minPresence = 1;
elseif isempty(minPresence)
    minPresence = 1;
end

%% Pull out only the genes variably present in the core models

if type == 1

    % Variable to hold genes present in at least one model
    genes = {};
    
    % Identify genes variably present
    for n = 1:length(model.genes)
        if isnan(sum(str2double(varPresence(n,1:end))))
            if sum(cell2mat(varPresence(n,1:end))) >= minPresence
                if sum(cell2mat(varPresence(n,1:end))) < size(varPresence,2)
                    genes = vertcat(genes,{model.genes{n}});
                end
            end
        else
            if sum(str2double(varPresence(n,1:end))) >= minPresence
                if sum(str2double(varPresence(n,1:end))) < size(varPresence,2)
                    genes = vertcat(genes,{model.genes{n}});
                end
            end
        end
    end
    
end

%% Pull out only the reactions variably present in the core models

if type == 2

    % Variable to hold reactions present in at least one model
    reactions = {};
    
    % Identify reactions variably present
    for n = 1:length(model.rxns)
        if isnan(sum(str2double(varPresence(n,1:end))))
            if sum(cell2mat(varPresence(n,1:end))) >= minPresence
                if sum(cell2mat(varPresence(n,1:end))) < size(varPresence,2)
                    reactions = vertcat(reactions,{model.rxns{n}});
                end
            end
        else
            if sum(str2double(varPresence(n,1:end))) >= minPresence
                if sum(str2double(varPresence(n,1:end))) < size(varPresence,2)
                    reactions = vertcat(reactions,{model.rxns{n}});
                end
            end
        end
    end    
    
end

%% Prepare a gene co-ocurrance matrix

if type == 1
   
    % Set output variable
    coocurMatrix = cell(length(genes),length(genes));
    
    % Transpose the input variable table
    variableTransposed = transpose(variable);
    
    % Find position of each gene in list
    geneIDs = num2cell(findGeneIDs(model,genes));
    
    % Count the cooccurrence of the genes
    for n = 1:length(geneIDs)
        nModels = strmatch(model.genes{geneIDs{n}},...
            variableTransposed(:,geneIDs{n}),'exact');
        for m = n:length(geneIDs)
            mModels = strmatch(model.genes{geneIDs{m}},...
                variableTransposed(:,geneIDs{m}),'exact');
            intersectModels = intersect(nModels,mModels);
            if isempty(intersectModels)
                count = 0;
            else
                count = length(intersectModels);
            end
            coocurMatrix{n,m} = count;
        end
    end
    
    % Make the matrix a square
    transposedCoocurMatrix = transpose(coocurMatrix);
    matrixIsEmpty = cellfun('isempty',coocurMatrix);
    coocurMatrix(matrixIsEmpty) = ...
        transposedCoocurMatrix(matrixIsEmpty);
    
end

%% Prepare a reaction cooccurrence matrix

if type == 2

    % Set output variable
    coocurMatrix = cell(length(reactions),length(reactions));
    
    % Transpose the input variable table
    variableTransposed = transpose(variable);
    
    % Find position of each gene in list
    rxnIDs = num2cell(findRxnIDs(model,reactions));
    
    % Count the cooccurrence of the genes
    for n = 1:length(rxnIDs)
        nModels = strmatch(model.rxns{rxnIDs{n}},...
            variableTransposed(:,rxnIDs{n}),'exact');
        for m = n:length(rxnIDs)
            mModels = strmatch(model.rxns{rxnIDs{m}},...
                variableTransposed(:,rxnIDs{m}),'exact');
            intersectModels = intersect(nModels,mModels);
            if isempty(intersectModels)
                count = 0;
            else
                count = length(intersectModels);
            end
            coocurMatrix{n,m} = count;
        end
    end
    
    % Make the matrix a square
    transposedCoocurMatrix = transpose(coocurMatrix);
    matrixIsEmpty = cellfun('isempty',coocurMatrix);
    coocurMatrix(matrixIsEmpty) = ...
        transposedCoocurMatrix(matrixIsEmpty);
    
end

%% Prepare a gene-growth matrix

if type == 1
    
    % Set storage variable
    growthMatrix = cell(3,size(variable,2));
    
    % Record average growth and number of models with each gene
    growthMean = mean(cell2mat(growth));
    for n = 1:length(model.genes)
        IDs = num2cell(strmatch(model.genes{n},variableTransposed(:,n),...
            'exact'));
        temp = cell(length(IDs),1);
        if ~isempty(IDs)
            for m = 1:length(IDs)
                temp{m,1} = growth{1,IDs{m}};
            end
            growthMatrix{1,n} = mean(cell2mat(temp));
            growthMatrix{3,n} = length(IDs);
            growthMatrix{2,n} = ((growthMean * size(variable,2)) - ...
                (growthMatrix{1,n} * growthMatrix{3,n})) / ...
                (size(variable,2) - growthMatrix{3,n});
        else
            growthMatrix{1,n} = [];
            growthMatrix{3,n} = [];
            growthMatrix{2,n} = [];
        end
    end
    
    % Remove the empty columns
    isNotEmpty = any(~cellfun('isempty',growthMatrix),1);
    growthMatrix = growthMatrix(:,isNotEmpty);
    
    % Remove data for genes in all models
    growthMatrixTemp = {};
    
    for n = 1:size(growthMatrix,2)
        if growthMatrix{3,n} ~= size(varPresence,2)
            growthMatrixTemp = horzcat(growthMatrixTemp,growthMatrix(:,n));
        end
    end
    
    growthMatrix = growthMatrixTemp;
    
end

%% Prepare a reaction-growth matrix

if type == 2
    
    % transpose the input variable table
    variableTransposed = transpose(variable);
    
    % Set storage variable
    growthMatrix = cell(3,size(variable,2));

    % Record average growth and number of models with each gene
    growthMean = mean(cell2mat(growth));
    for n = 1:length(model.rxns)
        IDs = num2cell(strmatch(model.rxns{n},variableTransposed(:,n),...
            'exact'));
        temp = cell(length(IDs),1);
        if ~isempty(IDs)
            for m = length(IDs)
                temp{m,1} = growth{1,IDs{m}};
            end
            growthMatrix{1,n} = mean(cell2mat(temp));
            growthMatrix{3,n} = length(IDs);
            growthMatrix{2,n} = ((growthMean * size(variable,2)) - ...
                (growthMatrix{1,n} * growthMatrix{3,n})) / ...
                (size(variable,2) - growthMatrix{3,n});
        else
            growthMatrix{1,n} = [];
            growthMatrix{3,n} = [];
            growthMatrix{2,n} = [];
        end
    end
    
    % Remove the empty columns
    isNotEmpty = any(~cellfun('isempty',growthMatrix),1);
    growthMatrix = growthMatrix(:,isNotEmpty);

    % Remove data for reactions in all models
    growthMatrixTemp = {};
    
    for n = 1:size(growthMatrix,2)
        if growthMatrix{3,n} ~= size(varPresence,2)
            growthMatrixTemp = horzcat(growthMatrixTemp,growthMatrix(:,n));
        end
    end
    
    growthMatrix = growthMatrixTemp;

end

%% Prepare adjusted gene or reaction matrix

% Variable to hold the data
coocurMatrixAdj = ...
    cell(length(coocurMatrix),length(coocurMatrix));
iters = size(variable,2);

% Calculate a chi statistic
for n = 1:length(coocurMatrix)
    for m = 1:length(coocurMatrix)
        observed = coocurMatrix{n,m};
        expected = (iters * ((growthMatrix{3,n} / iters) * (growthMatrix{3,m} / iters)));
        coocurMatrixAdj{n,m} = (observed - expected)^2 / expected;
        if observed < expected
            coocurMatrixAdj{n,m} = coocurMatrixAdj{n,m} * -1;
        end 
    end
end

%% And gene or reaction names to the matrixes

% Collect the correct labels
if type == 1
    labels = genes;
    labels2 = model.genes;
else
    labels = reactions;
    labels2 = model.rxns;
end

% Add labels to the rows for the co-ocurrance matrixes
coocurMatrix = horzcat(labels,coocurMatrix);
coocurMatrixAdj = horzcat(labels,coocurMatrixAdj);

% Add labels to the columns
labelsTrans = transpose(labels);
labelsTrans = horzcat({[]},labelsTrans);
coocurMatrix = vertcat(labelsTrans,coocurMatrix);
coocurMatrixAdj = vertcat(labelsTrans,coocurMatrixAdj);

% Add labels to the rows of the growth matrix
growthRows = {'GrowthWithGene';'GrowthWithOutGene';'Incidence'};
growthMatrix = horzcat(growthRows,growthMatrix);

% Add labels to the columns of the growth matrix
growthMatrix = vertcat(labelsTrans,growthMatrix);

%% Modify the varPresence matrix

% Add labels to the rows of the varPresence binary input
varPresenceTemp = horzcat(labels2,varPresence);

% Add labels to the columns of the varPresence input
headers = cell(1,size(varPresenceTemp,2));

for n = 1:size(varPresence,2)
    temp = ['Model_' num2str(n)];
    headers{1,n+1} = temp;
end

varPresenceTemp = vertcat(headers,varPresenceTemp);

% Remove genes/rxns in no core model
varPresenceNew = {};

for n = 2:size(varPresenceTemp,1)
    if isnan(sum(str2double(varPresenceTemp(n,2:end))))
        if sum(cell2mat(varPresenceTemp(n,2:end))) ~= 0
            varPresenceNew = vertcat(varPresenceNew,varPresenceTemp(n,:));
        end
    else  
        if sum(str2double(varPresenceTemp(n,2:end))) ~= 0
            varPresenceNew = vertcat(varPresenceNew,varPresenceTemp(n,:));
        end
    end
end

varPresenceNew = vertcat(headers,varPresenceNew);

%% Number of unique models

if isnan(sum(str2double(varPresence(1,1:end))))
    uniqueModels = size(unique(cell2mat(transpose(varPresence)),'rows'),1);
else
    uniqueModels = size(unique(str2double(transpose(varPresence)),'rows'),1);
end

