function tncore_init(solver)

%
% An optional simple function to initialize the workspace
%
% USAGE
%   tncore_init(solver)
% 
% OPTIONAL INPUTS:
%   solver      The solver to be used (Default = 'gurobi')
%
% AUTHORS
%   George diCenzo and Marco Fondi - 06/04/2018
%

%% Initialize everything

if nargin == 0
    solver = 'gurobi';
elseif isempty(solver)
    solver = 'gurobi';
end

initCobraToolbox;
changeCobraSolver(solver, 'all');
start_tiger(solver);
