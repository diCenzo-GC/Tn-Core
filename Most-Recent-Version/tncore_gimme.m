function [gene_states, genes, sol, tiger] = tncore_gimme(tiger, express, ...
    thresh, gene_names, obj_frac)

%
% A slight modification of the TIGER Toolbox (Jensen PA, Lutz KA, Papin JA. 
% 2011. BMC Syst Biol. 5:147) implementation of GIMME (Becker SA, Palsson 
% BØ. 2008. PLOS Comput Biol. 4:e1000082). After solving the MIP problem, 
% instead of rounding the gene states, all gene states above zero are set 
% to one, and the rest are set to zero. Requires the TIGER Toolbox to work.
%
% USAGE
%   [gene_states, genes, sol, tiger] = tncore_gimme(tiger, express, ...
%       thresh, gene_names, obj_frac)
%
% INPUTS
%   tiger           The TIGER model from which the context-specific model 
%                   is to be extracted
%   express         Expression data for all model genes (or properly
%                   modified essentiality data)
%
% OPTIONAL INPUTS
%   thresh          The threshold for a gene to be considered expressed
%                   (Default = 0.02% of the sum of all expression values).
%                   Note, if using essentiality (Tn-seq) data, this
%                   variable should be manually set in the inputs.
%   gene_names      Cell array of the gene names in the same order as their
%                   expression data in the express variable
%   obj_frac        The minimum allowable objective flux in the generated
%                   core model, as a fraction of growth of the input model
%                   (Default = 0.5)
%
% OUTPUTS
%
%   gene_states     Binary expression states calculated by GIMME
%   genes           Cell array of the gene names corresponding to
%                   gene_states
%   sol             The solution structure with details from the MILP
%                   solver
%   tiger           TIGER model with gene_states applied
%

%% Check input variables

% Check input is a TIGER formatted model
tiger = assert_tiger(tiger);

% Check there are enough inputs
assert(nargin >= 2,'GIMME requires at least two inputs');

% Set default expression threshold
if nargin < 3
    thresh = sum(express) / 5000;
elseif isempty(thresh)
    thresh = sum(express) / 5000;
end

% Set default for gene_names
if nargin < 4
    gene_names = tiger.genes;
elseif isempty(gene_names)
    gene_names = tiger.genes;
end
genes = gene_names;

% Check length of gene_names is correct
assert(length(express) == length(gene_names), 'gene_names must match size of EXPRESS');

% Set default for obj_frac
if nargin < 5
    obj_frac = 0.5;
elseif isempty(gene_names)
    obj_frac = 0.5;
end
frac = obj_frac;

%% Run GIMME

off_locs = express < thresh;
off_genes = genes(off_locs);
w = thresh - express(off_locs);

[~,off_idxs] = convert_ids(tiger.varnames,off_genes);

if frac > 0
    model = add_growth_constraint(tiger,frac);
else
    model = tiger;
end
model.obj(:) = 0;
model.obj(off_idxs) = w;
sol = cmpi.solve_mip(model);

%% Set gene states

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

