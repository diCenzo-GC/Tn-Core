function [varPresenceMatrix] = ...
            tncore_compare(varNames, varPresenceArray, headers)

% To create a matrix indicating the percent of models containing a gene or
% reaction between core model generation attempts
%
% USAGE
%   [varPresenceMatrix] = ...
%       tncore_compare(geneNames, varPresenceArray, headers, type);
%
% INPUTS
%   geneNames           A cell array of the variable names in the same 
%                       order as the data for the variable presence/absence
%                       in the varPresence arrays.
%   varPresenceArray    A cell array containing the names of the cell
%                       arrays holding the variable presence/absence data
%
% OPTIONAL INPUTS
%   headers             A cell array containing names to use as the headers
%                       in the output file (Default = varPresenceArray)
%
% OUTPUTS
%   varPresenceMatrix   A matrix indicating the percent of models that
%                       contain the gene/reaction for each set of models
%
% AUTHORS
%   George diCenzo and Marco Fondi - 28/09/2017

%% Check inputs

% Check there are enough inputs
assert(nargin >= 2, 'This function requires at least two inputs');

% Set default headers
if nargin < 3
    headers = varPresenceArray;
elseif isempty(headers)
    headers = varPresenceArray;
end

%% Count gene occurence and prepare matrix

varPresenceMatrix = cell(size(evalin('base',varPresenceArray{1,1}),1),...
    length(varPresenceArray));

for n = 1:length(varPresenceArray)
    genePresenceTemp = evalin('base',varPresenceArray{n});
    if ~isnan(sum(str2double(genePresenceTemp(1,:))))
        for m = 1:size(genePresenceTemp,1)
            varPresenceMatrix{m,n} = sum(str2double(genePresenceTemp(m,:)));
        end
    else
        for m = 1:size(genePresenceTemp,1)
            varPresenceMatrix{m,n} = sum(cell2mat(genePresenceTemp(m,:)));
        end
    end
end

%% Add headers and gene names to matrix

% Add gene names
varPresenceMatrix = horzcat(varNames,varPresenceMatrix);

% Add headers
headers = horzcat({[]},headers);
varPresenceMatrix = vertcat(headers,varPresenceMatrix);
