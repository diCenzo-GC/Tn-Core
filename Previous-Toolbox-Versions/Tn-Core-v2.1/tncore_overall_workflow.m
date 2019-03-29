function tncore_overall_workflow(varargin)

%
% An overall pipeline for running Tn-Core, from initializing the workspace
% and importing the data, through to exporting of the output.
%
% USAGE
%   tncore_overall_workflow(varargin)
%
% BASIC SET-UPS
%   For preparation of a core model without analysis of redundancy:
%   tncore_overall_workflow()
%
%   For preparation of a core model with analysis of redundancy:
%   tncore_overall_workflow('testRedun', true)
%
% OPTIONAL INPUTS
%   solver          Solver to be used with the COBRA and TIGER Toolboxes. 
%                   (Default = 'gurobi')
%   prepCore        A true/false variable indicating if a context-specific 
%                   core metebolic model should be prepared. 
%                   (Default = true)
%   testRedun       A true/false variable indicating if an analysis of
%                   redundancy in core metabolic pathways should be
%                   performed. (Default = false)
%   tnseqData       A true/false variable indicating if a Tn-seq data is 
%                   available in the current working directory for core 
%                   model generation. (Default = true)
%   rnaseqData      A true/false variable indicating if a RNA-seq data is 
%                   available in the current working directory for core 
%                   model generation. (Default = false)
%   epsilon         The epsilon value (the minimum flux rate that is
%                   considered non-zero) to be used during the running of
%                   FASTCC and FASTCORE when preparing a context-specific 
%                   core model. (Default = 1.1e-6)
%   binThresh       An array containing three values (number of standard
%                   deviations away from the mean of the log of the Tn-seq 
%                   data) to use for setting the limits of each bin. 
%                   (Default = {'3.5'; '2.5'; '1.5'}) 
%   essThresh       The threshold (in number of standard deviations from
%                   the mean of the log of the Tn-seq data) for a gene to 
%                   be considered essential when running the redundancy 
%                   analysis. (Default = 3.5)
%   growthFrac      The minimum allowable objective flux in the generated
%                   core model, as a fraction of growth of the input model
%                   (Default = 0.5) 
%   expressThresh   The threshold for a gene to be considered expressed
%                   (Default = 0.02% of the sum of all expression values)
%   keepHigh        A true/false variable indicating if highly expressed 
%                   genes should be added back into the context-specific 
%                   model, after the context-specific model has otherwise 
%                   been prepared. (Default = false)
%   deadends        A true/false variable indicating if reactions producing
%                   dead-end metabolites should be removed from the core 
%                   model. (Default = true)
%   iters           Number of random core models to generate, if choosing
%                   to examine redundancy. (Default = 1000)
%   redunMethod     Indicates if FBA or MOMA should be used during core
%                   model generation, if choosing to examine redundancy. 
%                   Use 1 for FBA and 2 for MOMA. (Default = 1)
%   tnseqData2      A true/false variable indicating if Tn-seq data 
%                   should be used when examining redundancy. If true, then
%                   the variable tnseqData must also be set to true.
%                   (Default = true)
%   rnaseqData2     A true/false variable indicating if RNA-seq data 
%                   should be used when examining redundancy. If true, then
%                   the variable rnaseqData must also be set to true.
%                   (Default = false)
% AUTHORS
%   George diCenzo and Marco Fondi - 06/04/2018
%   George diCenzo and Marco Fondi - 12/11/2018
%

%% Parse input variables

% Set up input parser
p = inputParser();

% Set up default values
defaultSolver = 'gurobi';
defaultTnseqData = true;
defaultRnaseqData = false;
defaultPrepCore = true;
defaultEpsilon = 1.1e-6;
defaultBinThresh = {'3.5'; '2.5'; '1.5'};
defaultEssThresh = 3.5;
defaultGrowthFrac = 0.5;
defaultDeadends = true;
defaultTestRedun = false;
defaultIters = 1000;
defaultRedunMethod = 1;
defaultExpressThresh = 0;
defaultKeepHigh = false;
defaultTnseqData2 = true;
defaultRnaseqData2 = false;

% Prepare parser
p.addParameter('solver', defaultSolver);
p.addParameter('tnseqData', defaultTnseqData);
p.addParameter('rnaseqData', defaultRnaseqData);
p.addParameter('prepCore', defaultPrepCore);
p.addParameter('epsilon', defaultEpsilon);
p.addParameter('binThresh', defaultBinThresh);
p.addParameter('essThresh', defaultEssThresh);
p.addParameter('growthFrac', defaultGrowthFrac);
p.addParameter('deadends', defaultDeadends);
p.addParameter('testRedun', defaultTestRedun);
p.addParameter('iters', defaultIters);
p.addParameter('redunMethod', defaultRedunMethod);
p.addParameter('expressThresh', defaultExpressThresh);
p.addParameter('keepHigh', defaultKeepHigh);
p.addParameter('tnseqData2', defaultTnseqData2);
p.addParameter('rnaseqData2', defaultRnaseqData2);

% Parse input variables
p.parse(varargin{:});

% Set up variables
solver = p.Results.solver;
tnseqData = p.Results.tnseqData;
rnaseqData = p.Results.rnaseqData;
prepCore = p.Results.prepCore;
epsilon = p.Results.epsilon;
binThresh = p.Results.binThresh;
essThresh = p.Results.essThresh;
growthFrac = p.Results.growthFrac;
deadends = p.Results.deadends;
testRedun = p.Results.testRedun;
iters = p.Results.iters;
redunMethod = p.Results.redunMethod;
expressThresh = p.Results.expressThresh;
keepHigh = p.Results.keepHigh;
tnseqData2 = p.Results.tnseqData2;
rnaseqData2 = p.Results.rnaseqData2;

% If preparing core model, check that Tn-seq data is provided
if prepCore == true
    assert(tnseqData == true, 'Core model can only be prepared if Tn-seq data is provided');
end

% Check if Tn-seq and RNA-seq properly provided
if testRedun == true
    if tnseqData2 == true
        assert(tnseqData == true, 'Please also set tnseqData to true');
    end
    if rnaseqData2 == true
        assert(rnaseqData == true, 'Please also set rnaseqData to true');
    end
end

% Change format of deadends
if deadends == true
    deadends = 1;
elseif deadends == false
    deadends = 0;
end

%% Initialize the system

% Turn things on and set solvers
tncore_init(solver);

%% Import the data

% Import the data
if tnseqData == true && rnaseqData == true
   [model, tnseq, rnaseq] = tncore_import();
elseif tnseqData == true && rnaseqData == false
   [model, tnseq] = tncore_import();
elseif tnseqData == false && rnaseqData == false
   [model] = tncore_import();
elseif tnseqData == false && rnaseqData == true
   [model, ~, rnaseq] = tncore_import();
end

% Set default RNA-seq exprssion threshold
if rnaseqData == 1 && expressThresh == 0
    expressThresh = sum(cell2mat(rnaseq(:,1))) / 5000;
end

% Set growth threshold
solutionOriginal = optimizeCbModel(model,'max');
growthThresh = growthFrac * solutionOriginal.f;

%% Prepare the core reconstruction

% Prepare the core model
if prepCore == true && rnaseqData == true
    [contextModel, reducedModel] = tncore_core(model, tnseq, epsilon, binThresh, ...
        growthFrac, rnaseq, expressThresh, keepHigh, deadends);
elseif prepCore == true && rnaseqData == false
    [contextModel, reducedModel] = tncore_core(model, tnseq, epsilon, binThresh, ...
        growthFrac, [], [], keepHigh, deadends);
end

%% Analysis of core model redundancy

% Run the redunancy analysis script
if testRedun == true 
    if tnseqData2 == true && rnaseqData2 == true
        [coreGeneVar, coreRxnVar, coreGrowthVar, ~, ~, reducedModel, coreGeneCat, ...
            genePresence, rxnPresence] = tncore_redundancy(model, ...
            iters, growthFrac, tnseq, {1}, essThresh, binThresh, rnaseq, expressThresh, redunMethod);
    elseif tnseqData2 == true && rnaseqData2 == false
        [coreGeneVar, coreRxnVar, coreGrowthVar, ~, ~, reducedModel, coreGeneCat, ...
            genePresence, rxnPresence] = tncore_redundancy(model, ...
            iters, growthFrac, tnseq, {1}, essThresh, binThresh, [], [], redunMethod);
    elseif tnseqData2 == false && rnaseqData2 == true
        [coreGeneVar, coreRxnVar, coreGrowthVar, ~, ~, reducedModel, coreGeneCat, ...
            genePresence, rxnPresence] = tncore_redundancy(model, ...
            iters, growthFrac, [], [], binThresh, rnaseq, expressThresh, redunMethod);
    elseif tnseqData2 == false && rnaseqData2 == false
        [coreGeneVar, coreRxnVar, coreGrowthVar, ~, ~, reducedModel, coreGeneCat, ...
            genePresence, rxnPresence] = tncore_redundancy(model, ...
            iters, growthFrac, [], [], binThresh, [], [], redunMethod);
    end
end

% Prepare the output matrixes
if testRedun == true
    [coocurMatrixGenes, coocurMatrixAdjGenes, growthMatrixGenes, ...
        genePresenceNew, uniqueModelsGenes] = tncore_matrix(reducedModel, ...
        coreGeneVar, coreGrowthVar, genePresence, 1, 1);
    [coocurMatrixRxns, coocurMatrixAdjRxns, growthMatrixRxns, ...
        rxnPresenceNew, uniqueModelsRxns] = tncore_matrix(reducedModel, ...
        coreRxnVar, coreGrowthVar, rxnPresence, 2, 1);
end

%% Export the data

% Export everything
if prepCore == true
    if tnseqData == true && rnaseqData == true
        if testRedun == true
            tncore_export(contextModel, reducedModel, tnseq, rnaseq, coocurMatrixGenes, ...
                coocurMatrixAdjGenes, growthMatrixGenes, iters, coocurMatrixRxns, ...
                coocurMatrixAdjRxns, growthMatrixRxns, genePresenceNew, rxnPresenceNew);
        elseif testRedun == false
            tncore_export(contextModel, reducedModel, tnseq, rnaseq);
        end
    elseif tnseqData == true && rnaseqData == false
        if testRedun == true
            tncore_export(contextModel, reducedModel, tnseq, [], coocurMatrixGenes, ...
                coocurMatrixAdjGenes, growthMatrixGenes, iters, coocurMatrixRxns, ...
                coocurMatrixAdjRxns, growthMatrixRxns, genePresenceNew, rxnPresenceNew);
        elseif testRedun == false
            tncore_export(contextModel, reducedModel, tnseq, []);
        end
    end
elseif prepCore == false
    if testRedun == true
        tncore_export([], [], [], [], coocurMatrixGenes, ...
            coocurMatrixAdjGenes, growthMatrixGenes, iters, coocurMatrixRxns, ...
            coocurMatrixAdjRxns, growthMatrixRxns, genePresenceNew, rxnPresenceNew);
    elseif testRedun == false
    end 
end

