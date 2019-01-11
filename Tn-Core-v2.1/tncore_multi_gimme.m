function [gene_states, genes, sol, tiger, weightStruct] = tncore_multiGIMME(tiger, ...
    fields, expressStruct, genesStruct, threshStruct, frac)

%
% An adaptation of tncore_gimme (which is a slight modification of the 
% TIGER Toolbox [Jensen PA, Lutz KA, Papin JA. 2011. BMC Syst Biol. 5:147] 
% implementation of GIMME [Becker SA, Palsson BØ. 2008. PLOS Comput Biol. 
% 4:e1000082]) to use with a multi-species model. Takes separate gene
% lists, expression data, and expression thresholds for each species in the
% model. Requires the TIGER Toolbox to work. Preparation of the input
% structures is described in the Tn-Core manual.
%
% USAGE
%   [gene_states, genes, sol, tiger, weightStruct] = tncore_multi_gumme(tiger, ...
%       fields, expressStruct, genesStruct, threshStruct, frac)
%
% INPUTS
%   tiger           The TIGER model from which the context-specific model 
%                   is to be extracted
%   fields          A cell array containing the names of each of the fields
%                   present in the expressStruct, geneStruct, and
%                   threshStruct fields
%   expressStruct   A structure with one matrix per species in the model,
%                   and containing the expression data for all model genes
%                   of the corresponding species
%   genesStruct     A structure with one cell array per species in the
%                   model, and containing the model gene names of the
%                   corresponding species in teh same order as in
%                   expressStruct
%
% OPTIONAL INPUTS
%   threshStruct    A structure with one matrix per species in the model,
%                   and containing the expression threshold for a gene in
%                   that species to be considered expressed.
%                   (Default = the mean expression of genes in the species)
%   obj_frac        The minimum allowable objective flux in the generated
%                   context specific model, as a fraction of growth of the 
%                   input model. (Default = 0.5)
%
% OUTPUTS
%   gene_states     Binary expression states calculated by GIMME
%   genes           Cell array of the gene names corresponding to
%                   gene_states
%   sol             The solution structure with details from the MILP
%                   solver
%   tiger           TIGER model with gene_states applied
%   weightStruct    A structure with one matrix per species in the model,
%                   and the weighting applied to the expression values of
%                   each species
%
% AUTHORS
%   George diCenzo and Marco Fondi - 10/12/2018
%

%% Check input

% Check tiger
tiger = assert_tiger(tiger);

% Ensure there are enough inputs
assert(nargin >= 4,'tncore_multiGIMME requires at least four inputs');
assert(numel(fieldnames(expressStruct)) == numel(fieldnames(genesStruct)), ...
   'number of fields in express must match fields in genes');
if ~isempty(threshStruct)
    assert(numel(fieldnames(expressStruct)) == numel(fieldnames(threshStruct)), ...
       'number of fields in express must match fields in thresh');
end
assert(numel(fieldnames(expressStruct)) == length(fields), ...
   'number of fields in express must match length of fields array');

% Check length of the gene and express structures
for n = 1:length(fields)
    assert(length(expressStruct.(fields{n})) == length(genesStruct.(fields{n})), ...
        'the length of each field in the expressStruct should be equal to the corresponding field in the genesStruct');
end
   
%% Set up default and initial parameters

% Set default thresholds
if nargin < 5
    threshStruct = struct();
    for n = 1:length(fields)
        threshStruct.(fields{n}) = mean(expressStruct.(fields{n}));
    end
elseif isempty(threshStruct)
    threshStruct = struct();
    for n = 1:length(fields)
        threshStruct.(fields{n}) = mean(expressStruct.(fields{n}));
    end
end

% Set default objective fraction
if nargin < 6
    frac = 0.5;
end

% Combine gene names
genes = {};
for n = 1:length(fields)
    genes = vertcat(genes, genesStruct.(fields{n}));
end

% Combine expression
express = {};
express = cell2mat(express);
for n = 1:length(fields)
    express = vertcat(express, expressStruct.(fields{n}));
end

%% Turn off genes

% Find off genes for each data set
offLocsStruct = struct();
offGenesStruct = struct();
for n = 1:length(fields)
    offLocsStruct.(fields{n}) = expressStruct.(fields{n}) < threshStruct.(fields{n});
    offGenesStruct.(fields{n}) = genesStruct.(fields{n})(offLocsStruct.(fields{n}));
end

% Combine off gene lists
off_locs = {};
off_genes = {};
off_locs = cell2mat(off_locs);
for n = 1:length(fields)
    off_locs = vertcat(off_locs, offLocsStruct.(fields{n}));
    off_genes = vertcat(off_genes, offGenesStruct.(fields{n}));
end
off_locs = logical(off_locs);

% Get weighting of each data set versus the first data set
weightStruct = struct();
weightStruct.(fields{1}) = 1;
for n = 2:length(fields)
    weightStruct.(fields{n}) = threshStruct.(fields{n}) / threshStruct.(fields{1});
end

% Calculate 'w' values for each dataset
wStruct = struct();
for n = 1:length(fields)
    wStruct.(fields{n}) = threshStruct.(fields{n}) - expressStruct.(fields{n})(offLocsStruct.(fields{n}));
end

% Adjust the 'w' values for each dataset
for n = 1:length(fields)
    wStruct.(fields{n}) = wStruct.(fields{n}) / weightStruct.(fields{n});
end

% Combine the 'w' values
w = {};
w = cell2mat(w);
for n = 1:length(fields)
    w = vertcat(w, wStruct.(fields{n})(:,1));
end

%% From original GIMME

[~,off_idxs] = convert_ids(tiger.varnames,off_genes);

if frac > 0
    model = add_growth_constraint(tiger,frac);
else
    model = tiger;
end
model.obj(:) = 0;
model.obj(off_idxs) = w;
sol = cmpi.solve_mip(model);

if cmpi.is_acceptable_exit(sol)
    gene_states = ones(size(express));
    x = 0;
    for n = 1:length(gene_states)
        if off_locs(n) == 1
            x = x + 1;
            if sol.x(off_idxs(x)) > 0
                gene_states(n,1) = 1;
            else
                gene_states(n,1) = 0;
            end
        end
    end
    tiger = set_var(tiger,genes,gene_states);
else
    gene_states = [];
end
