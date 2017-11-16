function [gene_states,genes,sol,tiger] = ...
                                      gimme_mine(tiger,express,thresh,varargin)
% GIMME  Gene Inactivity Moderated by Metabolism and Expression
%
%   [GENE_STATES,GENES,SOL,TIGER] = 
%       GIMME(TIGER,EXPRESS,THRESH,...parameters...)
%
%   Integrate expression data with a metabolic model using a modified
%   GIMME algorithm [Becker & Palsson (2008), PLoS Comput Biol].  This
%   implementation uses an integrated model to map expression values to
%   the gene instead of averaging over each reaction.
%
%   Inputs
%   TIGER    TIGER model.  COBRA models will be converted to TIGER models
%            with a warning.
%   EXPRESS  Expression value for each gene.
%   THRESH   Threshold for genes to be turned "on".  If not given, the
%            default is mean(EXPRESS).
%
%   Parameters
%   'gene_names'  Cell array of names for genes in EXPRESS.  If not given,
%                 the default is TIGER.genes.
%   'obj_frac'    Fraction of metabolic objective required in the 
%                 resulting model (v_obj >= frac*v_obj_max). 
%                 Default is 0.3.
%
%   Outputs
%   GENE_STATES  Binary expression states calculated by GIMME.
%   GENES        Cell of gene names corresponding to GENE_STATES.
%   SOL          Solution structure with details from the MILP solver.
%   TIGER        TIGER model with GENE_STATES applied.

tiger = assert_tiger(tiger);

p = inputParser;
p.addParamValue('gene_names',tiger.genes);
p.addParamValue('obj_frac',0.3);
p.parse(varargin{:});

genes = p.Results.gene_names;
frac = p.Results.obj_frac;

assert(nargin >= 2,'GIMME requires at least two inputs');
assert(length(express) == length(genes), ...
       'gene_names must match size of EXPRESS');
   
if nargin < 3
    thresh = mean(express);
end

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

