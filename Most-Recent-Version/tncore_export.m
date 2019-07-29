function tncore_export(model, reducedModel, tnseq, rnaseq, coMatGenes, ...
    coMatGenesAdj, geneGrowMat, iters, coMatRxns, coMatRxnsAdj, rxnGrowMat, ...
    genePres, rxnPres)

%
% A function to export the outputs of Tn-Core. Exports one or more of the 
% following items:
% 1) The core metabolic model as an Excel file
% 2) Several files related to the occurrence and co-occurrence of each 
%    gene, as text files
% 3) Several files related to the occurrence and co-occurrence of each 
%    reaction, as text files
% 4) Files related to gene presence/absece, as text files
% 5) Files related to reaction presence/absece, as text files
%
% USAGE
%   tncore_export(model, reducedModel, tnseq, rnaseq, coMatGenes, ...
%       coMatGenesAdj, geneGrowMat, iters, coMatRxns, coMatRxnsAdj, ...
%       rxnGrowMat, genePres, rxnPres)
%
% INPUTS FOR MODEL EXPORT
%   model           The core metabolic model returned by the tncore_main
%                   function
%   reducedModel    The reduced genome-scale metabolic model returned by 
%                   the tncore_main function
%
% OPTIONAL INPUTS FOR MODEL EXPORT
%   tnseq           The tnseq data used during core model preparation
%   rnaseq          The RNA-seq data used during core model preparation
%
% INPUTS FOR GENE OR REACTION (CO)OCURENCE AND PRESENCE TABLE AND MATRIX EXPORT
%   iters           The number of iterations used in the tncore_redundancy
%                   function
%
% INPUTS FOR GENE (CO)OCURENCE TABLE AND MATRIX EXPORT
%   coMatGenes      The genes coocurence matrix produced by the
%                   tncore_matrix function
%   coMatGenesAdj   The adjusted gene coocurence matrix (with the
%                   chi-squared statistics) produced by the tncore_matrix
%                   function
%   geneGrowMat     The genes growth matrix produced by the tncore_matrix
%                   function
%
% INPUTS FOR REACTION (CO)OCURENCE TABLE AND MATRIX EXPORT
%   coMatRxns       The reaction coocurence matrix produced by the
%                   tncore_matrix function
%   coMatRxnsAdj    The adjusted reaction coocurence matrix (with the
%                   chi-squared statistics) produced by the tncore_matrix
%                   function
%   rxnGrowMat      The reaction growth matrix produced by the tncore_matrix
%                   function
%
% INPUTS FOR GENE PRESENCE TABLE EXPORT
%   genePres        The gene presence variable returned by the
%                   tncore_matrix function
%
% INPUTS FOR GENE PRESENCE TABLE EXPORT
%   rxnPres         The reaction presence variable returned by the
%                   tncore_matrix function
%
% AUTHORS
%   George diCenzo and Marco Fondi - 06/04/2018
%

%% Check the input variables

% Check that there is at least one input
assert(nargin >= 1,'Please provide at least one input');

% Check if there is a model
if isempty(model)
    isModel = 0;
else
    isModel = 1;
end

% Check the original model is provided
if isModel == 1
    assert(nargin >= 2,'Please provide the reduced model as well');
    assert(~isempty(reducedModel),'Please provide the reduced model as well');
end

% Check if there is tnseq data
if nargin < 3
    isTnseq = 0;
elseif isempty(tnseq)
    isTnseq = 0;
else
    isTnseq = 1;
end

% Check if there is rnaseq data
if nargin < 4
    isRnaseq = 0;
elseif isempty(rnaseq)
    isRnaseq = 0;
else
    isRnaseq = 1;
end

% Check if gene coocurence matrix
if nargin < 5
    isGeneCo = 0;
elseif isempty(coMatGenes)
    isGeneCo = 0;
else
    isGeneCo = 1;
end

% Check if gene adjusted coocurence matrix
if nargin < 6
    isGeneCoAdj = 0;
elseif isempty(coMatGenesAdj)
    isGeneCoAdj = 0;
else
    isGeneCoAdj = 1;
end

% Check if both gene coocurence matrixes are provided
if isGeneCo == 1
    assert(nargin >= 6,'Please provide the gene adjusted coocurence matrix as well');
    assert(isGeneCoAdj == 1,'Please provide the gene adjusted coocurence matrix as well');
end
if isGeneCoAdj == 1
    assert(isGeneCo == 1,'Please provide the unadjusted gene coocurence matrix as well');
end

% Check if growth  matrix is provided
if isGeneCo == 1
    assert(nargin >= 7,'Please provide the gene growth matrix as well');
    assert(~isempty(geneGrowMat),'Please provide the gene growth matrix as well');
end

% Check if rxn coocurence matrix
if nargin < 9
    isRxnCo = 0;
elseif isempty(coMatRxns)
    isRxnCo = 0;
else
    isRxnCo = 1;
end

% Check if rxn adjusted coocurence matrix
if nargin < 10
    isRxnCoAdj = 0;
elseif isempty(coMatRxnsAdj)
    isRxnCoAdj = 0;
else
    isRxnCoAdj = 1;
end

% Check if both rxn coocurence matrixes are provided
if isRxnCo == 1
    assert(nargin >= 10,'Please provide the gene adjusted coocurence matrix as well');
    assert(isRxnCoAdj == 1,'Please provide the gene adjusted coocurence matrix as well');
end
if isRxnCoAdj == 1
    assert(isRxnCo == 1,'Please provide the unadjusted gene coocurence matrix as well');
end

% Check if reaction growth matrix is provided
if isGeneCo == 1
    assert(nargin >= 11,'Please provide the reaction growth matrix as well');
    assert(~isempty(rxnGrowMat),'Please provide the reaction growth matrix as well');
end

% Check if iters is provided
if isGeneCo == 1
    assert(nargin >= 8,'Please provide iterations of tncore_redundancy that were run');
    assert(~isempty(iters),'Please provide iterations of tncore_redundancy that were run')
elseif isRxnCo == 1
    assert(~isempty(iters),'Please provide iterations of tncore_redundancy that were run')
end

% Check if gene presence matrix
if nargin < 12
    isGenePres = 0;
elseif isempty(genePres)
    isGenePres = 0;
else
    isGenePres = 1;
end

% Check if reaction presence matrix
if nargin < 13
    isRxnPres = 0;
elseif isempty(rxnPres)
    isRxnPres = 0;
else
    isRxnPres = 1;
end

%% Export the model reactions

if isModel == 1
    
    % Get gene-protein pairings
    if isfield(reducedModel, 'rxnNotes')
        
        rxnNotesNew = reducedModel.rxnNotes;
        proteinsAll = {};
        for n = 1:length(rxnNotesNew)
            if strfind(rxnNotesNew{n}, 'PROTEIN_ASSOCIATION')
                rxnNotesTemp = strsplit(rxnNotesNew{n}, ';');
                for m = 1:length(rxnNotesTemp)
                    if strmatch('PROTEIN_ASSOCIATION', rxnNotesTemp{m})
                        rxnNotesTemp{m} = strrep(rxnNotesTemp{m}, 'PROTEIN_ASSOCIATION:', '');
                        proteinsTemp = strsplit(rxnNotesTemp{m}, ' ');
                        for x = 1:length(proteinsTemp)
                            proteinsTemp{x} = strrep(proteinsTemp{x}, '(', '');
                            proteinsTemp{x} = strrep(proteinsTemp{x}, ')', '');
                            if strmatch('or', proteinsTemp{x}, 'exact')
                            elseif strmatch('and', proteinsTemp{x}, 'exact')
                            else
                                proteinsAll = vertcat(proteinsAll, proteinsTemp{x});
                            end
                        end
                    end
                end
            end
        end
        
        grRulesNew = reducedModel.grRules;
        genesAll = {};
        for n = 1:length(grRulesNew)
            genesTemp = strsplit(grRulesNew{n}, ' ');
            for x = 1:length(genesTemp)
                genesTemp{x} = strrep(genesTemp{x}, '(', '');
                genesTemp{x} = strrep(genesTemp{x}, ')', '');
                if strmatch('or', genesTemp{x}, 'exact')
                elseif strmatch('and', genesTemp{x}, 'exact')
                else
                    genesAll = vertcat(genesAll, genesTemp{x});
                end
            end
        end
        
        if length(genesAll) == length(proteinsAll)
            combinedAll = horzcat(genesAll, proteinsAll);
            [~, idx] = unique(combinedAll(:,1));
            geneProteinPairs = combinedAll(idx,:);
            geneProteinExists = true;
        else
            geneProteinExists = false;
        end
        
    end
    
    % Get reaction names
    RxnID = model.rxns;
    RxnName = model.rxnNames;
    
    % Get reaction formulas
    Reaction_A = printRxnFormula(model);
    Reaction_B = Reaction_A;
    for n = 1:length(Reaction_B)
        for m = 1:length(model.mets)
            Reaction_B{n} = strrep(Reaction_B{n}, model.mets{m}, model.metNames{m});
        end
    end
    
    % Get the reversibility
    Reversible = cell(length(model.rxns), 1);
    for n = 1:length(model.rxns)
        if model.rev(n) == 1
            Reversible{n} = 'true';
        elseif model.rev(n) == 0
            Reversible{n} = 'false';
        end
    end
    
    % Modify the rxnNotes field
    if isfield(model, 'rxnNotes')
        rxnNotes = model.rxnNotes;
        for n = 1:length(rxnNotes)
            rxnNotes{n} = strrep(rxnNotes{n}, '  ', ' ');
            rxnNotes{n} = strrep(rxnNotes{n}, '  ', ' ');
            rxnNotes{n} = strrep(rxnNotes{n}, '  ', ' ');
            rxnNotes{n} = strrep(rxnNotes{n}, ' or ', 'ZZZZZ');
            rxnNotes{n} = strrep(rxnNotes{n}, ' and ', 'XXXXX');
            rxnNotes{n} = strrep(rxnNotes{n}, ') )', '))');
            rxnNotes{n} = strrep(rxnNotes{n}, '( (', '))');
            rxnNotes{n} = strrep(rxnNotes{n}, ' ', ';');
            rxnNotes{n} = strrep(rxnNotes{n}, 'ZZZZZ', ' or ');
            rxnNotes{n} = strrep(rxnNotes{n}, 'XXXXX', ' and ');
        end
    end
    
    % Get KEGG data
    if isfield(model, 'rxnNotes')
        KEGG_RID = cell(length(rxnNotes), 1);
        for n = 1:length(rxnNotes)
            if strfind(rxnNotes{n}, 'KEGG_RID')
                keggTemp = strsplit(rxnNotes{n}, ';');
                for m = 1:length(keggTemp)
                    if strmatch('KEGG_RID', keggTemp{m})
                        KEGG_RID{n} = strrep(keggTemp{m}, 'KEGG_RID:', '');
                    end
                end
            end
        end
    end
    
    % Get protein class data
    if isfield(model, 'rxnNotes')
        EnzymeClass = cell(length(rxnNotes), 1);
        for n = 1:length(rxnNotes)
            if strfind(rxnNotes{n}, 'PROTEIN_CLASS')
                protTemp = strsplit(rxnNotes{n}, ';');
                for m = 1:length(protTemp)
                    if strmatch('PROTEIN_CLASS', protTemp{m})
                        EnzymeClass{n} = strrep(protTemp{m}, 'PROTEIN_CLASS:', '');
                    end
                end
            end
        end
    end
    
    % Get the gene associations
    Genes = model.grRules;
    if geneProteinExists == true
        Proteins = Genes;
        for n = 1:length(Proteins)
            for m = 1:length(geneProteinPairs)
                Proteins{n} = strrep(Proteins{n}, geneProteinPairs{m,1}, geneProteinPairs{m,2});
            end
        end
    end
    
    % Combine the reaction information
    if isfield(model, 'rxnNotes')
        if geneProteinExists == true
            reactionTable = table(RxnID, RxnName, Reaction_A, Reaction_B, Reversible, ...
                KEGG_RID, EnzymeClass, Genes, Proteins);
        else
            reactionTable = table(RxnID, RxnName, Reaction_A, Reaction_B, Reversible, ...
                KEGG_RID, EnzymeClass, Genes);
        end
    else
        reactionTable = table(RxnID, RxnName, Reaction_A, Reaction_B, Reversible, ...
            Genes);
    end
    
    % Export the reaction table
    writetable(reactionTable, 'exportedModel.xlsx', 'Sheet', 1);
    
end

%% Export the model genes

if isModel == 1
    
    % Get the genes
    Gene = model.genes;
    
    % Get the proteins
    if geneProteinExists == true
        Protein = Gene;
        for n = 1:length(Protein)
            for m = 1:length(geneProteinPairs)
                Protein{n} = strrep(Protein{n}, geneProteinPairs{m,1}, geneProteinPairs{m,2});
            end
        end
    end
    
    % Get the Tn-seq data
    if isTnseq == 1
        
        tnseqOut = cell(length(Gene),1);
        for n = 1:length(Gene)
            idx = strmatch(Gene{n}, tnseq(:,2), 'exact');
            if ~isempty(idx)
                tnseqOut{n} = tnseq{idx,1};
            else
                tnseqOut{n} = 'None';
            end
        end
        Tn_Seq_Data = tnseqOut;
        
    end
    
    % Get the RNA-seq data
    if isRnaseq == 1
        
        rnaseqOut = cell(length(Gene),1);
        for n = 1:length(Gene)
            idx = strmatch(Gene{n}, rnaseq(:,2), 'exact');
            if ~isempty(idx)
                rnaseqOut{n} = rnaseq{idx,1};
            else
                rnaseqOut{n} = 'None';
            end
        end
        RNA_Seq_Data = rnaseqOut;
        
    end
    
    % Combine the gene information
    if geneProteinExists == true
        if isTnseq == 1 && isRnaseq == 1
            geneTable = table(Gene, Protein, Tn_Seq_Data, RNA_Seq_Data);
        elseif isTnseq == 1 && isRnaseq == 0
            geneTable = table(Gene, Protein, Tn_Seq_Data);
        elseif isTnseq == 0 && isRnaseq == 1
            geneTable = table(Gene, Protein, RNA_Seq_Data);
        elseif isTnseq == 0 && isRnaseq == 0
            geneTable = table(Gene, Protein);
        end
    else
        if isTnseq == 1 && isRnaseq == 1
            geneTable = table(Gene, Tn_Seq_Data, RNA_Seq_Data);
        elseif isTnseq == 1 && isRnaseq == 0
            geneTable = table(Gene, Tn_Seq_Data);
        elseif isTnseq == 0 && isRnaseq == 1
            geneTable = table(Gene, RNA_Seq_Data);
        elseif isTnseq == 0 && isRnaseq == 0
            geneTable = table(Gene);
        end
    end
    
    % Export the gene table
    writetable(geneTable, 'exportedModel.xlsx', 'Sheet', 2);
    
end

%% Export the model metabolites

if isModel == 1
    
    % Get the metabolites and their names
    MetIDs = model.mets;
    MetNames = model.metNames;
    
    % Get metabolite charge, if present
    if isfield(model, 'metCharge')
        if ~isempty(model.metCharge)
            Charge = model.metCharge;
        else
            Charge = [];
        end
            Charge = [];
    end
    
    % Combine the metabolite information
    if ~isempty(Charge)
        metTable = table(MetIDs, MetNames, Charge);
    else
        metTable = table(MetIDs, MetNames);
    end

    % Export the metabolite table
    writetable(metTable, 'exportedModel.xlsx', 'Sheet', 3);

end

%% Export gene co-occurance table

if isGeneCo == 1

    % Get occurance of each gene in the reduced model
    for n = 1:size(geneGrowMat, 2)-1
        geneCount{n,1} = geneGrowMat{1, n+1};
        geneCount{n,2} = geneGrowMat{4, n+1};
    end

    % Get list of variable genes
    GeneA = coMatGenes(2:end, 1);
    GeneB = transpose(coMatGenes(1, 2:end));

    % Make the coocurence array
    Coocurence_Table = cell(nchoosek(length(GeneA), 2), 8);
    x = 0;
    for n = 1:(length(GeneA)-1)
        for m = n+1:length(GeneB)
            x = x + 1;
            posA = strmatch(GeneA{n}, geneCount(:,1), 'exact');
            posB = strmatch(GeneB{m}, geneCount(:,1), 'exact');
            Coocurence_Table{x,1} = GeneA{n};
            Coocurence_Table{x,2} = GeneB{m};
            Coocurence_Table{x,3} = geneCount{posA,2};
            Coocurence_Table{x,4} = geneCount{posB,2};
            Coocurence_Table{x,5} = coMatGenes{n+1,m+1};
            Coocurence_Table{x,6} = ...
                ((geneCount{posA,2} / iters) * (geneCount{posB,2} / iters)) * iters;
            Coocurence_Table{x,7} = Coocurence_Table{x,5} / Coocurence_Table{x,6};
            Coocurence_Table{x,8} = coMatGenesAdj{n+1,m+1};
        end
    end
    Coocurence_Table = sortrows(Coocurence_Table, 8);

    % Prepare the table to export
    Coocurence_Table_export = cell2table(Coocurence_Table);
    Coocurence_Table_export.Properties.VariableNames = ...
        {'Gene_1', 'Gene_2', 'Gene_1_Occurrence', 'Gene_2_Occurrence', ...
        'Cooccurrence_Count', 'Expected_Cooccurrence', 'Expected_Ratio', ...
        'Chi_Squared_Metric'};
    
    % Modify gene presence table
    for n = 1:length(geneCount)
        geneCount{n,3} = 100 * geneCount{n,2} / iters;
    end
    geneCount = vertcat({'Gene', 'Occurence_Count', 'Occurence_Percent'}, geneCount);
    
    % Prepare other tables for export
    coMatGenes_export = cell2table(coMatGenes);
    coMatGenesAdj_export = cell2table(coMatGenesAdj);
    geneCount_export = cell2table(geneCount);

    % Export tables
    writetable(Coocurence_Table_export, 'geneCoocurrenceTable.txt', 'Delimiter', '\t');
    writetable(coMatGenes_export, 'geneCoocurrenceMatrix.txt', 'Delimiter', '\t', ...
        'WriteRowNames', false, 'WriteVariableNames', false);
    writetable(coMatGenesAdj_export, 'geneCoocurrenceAdjustedMatrix.txt', 'Delimiter', '\t', ...
        'WriteRowNames', false, 'WriteVariableNames', false);

end

%% Export reaction co-occurance table

if isRxnCo == 1

    % Get occurance of each reaction in the reduced model
    for n = 1:size(rxnGrowMat, 2)-1
        rxnCount{n,1} = rxnGrowMat{1, n+1};
        rxnCount{n,2} = rxnGrowMat{4, n+1};
    end

    % Get list of variable rxns
    RxnA = coMatRxns(2:end, 1);
    RxnB = transpose(coMatRxns(1, 2:end));

    % Make the coocurence array
    Coocurence_Table = cell(nchoosek(length(RxnA), 2), 8);
    x = 0;
    for n = 1:(length(RxnA)-1)
        for m = n+1:length(RxnB)
            x = x + 1;
            posA = strmatch(RxnA{n}, rxnCount(:,1), 'exact');
            posB = strmatch(RxnB{m}, rxnCount(:,1), 'exact');
            Coocurence_Table{x,1} = RxnA{n};
            Coocurence_Table{x,2} = RxnB{m};
            Coocurence_Table{x,3} = rxnCount{posA,2};
            Coocurence_Table{x,4} = rxnCount{posB,2};
            Coocurence_Table{x,5} = coMatRxns{n+1,m+1};
            Coocurence_Table{x,6} = ...
                ((rxnCount{posA,2} / iters) * (rxnCount{posB,2} / iters)) * iters;
            Coocurence_Table{x,7} = Coocurence_Table{x,5} / Coocurence_Table{x,6};
            Coocurence_Table{x,8} = coMatRxnsAdj{n+1,m+1};
        end
    end
    Coocurence_Table = sortrows(Coocurence_Table, 8);

    % Prepare the table to export
    Coocurence_Table_export = cell2table(Coocurence_Table);
    Coocurence_Table_export.Properties.VariableNames = ...
        {'Reaction_1', 'Reaction_2', 'Reaction_1_Occurrence', ...
        'Reaction_2_Occurrence', 'Cooccurrence_Count', 'Expected_Cooccurrence', ...
        'Expected_Ratio', 'Chi_Squared_Metric'};
    
    % Modify rxn presence table
    for n = 1:length(rxnCount)
        rxnCount{n,3} = 100 * rxnCount{n,2} / iters;
    end
    rxnCount = vertcat({'Reaction', 'Occurence_Count', 'Occurance_Percent'}, rxnCount);
    
    % Prepare other tables for export
    coMatRxns_export = cell2table(coMatRxns);
    coMatRxnsAdj_export = cell2table(coMatRxnsAdj);
    rxnCount_export = cell2table(rxnCount);

    % Export tables
    writetable(Coocurence_Table_export, 'rxnCoocurrenceTable.txt', 'Delimiter', '\t');
    writetable(coMatRxns_export, 'rxnCoocurrenceMatrix.txt', 'Delimiter', '\t', ...
        'WriteRowNames', false, 'WriteVariableNames', false);
    writetable(coMatRxnsAdj_export, 'rxnCoocurrenceAdjustedMatrix.txt', 'Delimiter', '\t', ...
        'WriteRowNames', false,'WriteVariableNames', false);

end

%% Export gene presence matrix

if isGenePres == 1

    geneCount = {};
    if isnan(sum(str2double(genePres(2,2:end))))
        for n = 2:size(genePres,1)
            presence = sum(cell2mat(genePres(n, 2:end)));
            geneCount{n-1,1} = genePres{n,1};
            geneCount{n-1,2} = presence;
            geneCount{n-1,3} = 100 * presence / iters;
        end
    else
        for n = 2:size(genePres,1)
            presence = sum(str2double(genePres(n, 2:end)));
            geneCount{n-1,1} = genePres{n,1};
            geneCount{n-1,2} = presence;
            geneCount{n-1,3} = 100 * presence / iters;
        end
    end

    % Export tables
    genePres = cell2table(genePres);
    writetable(genePres, 'genePresenceMatrix.txt', 'Delimiter', '\t', ...
        'WriteRowNames', false, 'WriteVariableNames', false);
    geneCount = vertcat({'Gene', 'Occurence_Count', 'Occurence_Percent'}, geneCount);
    geneCount = cell2table(geneCount);
    writetable(geneCount, 'geneCountTable.txt', 'Delimiter', '\t', ...
           'WriteRowNames', false, 'WriteVariableNames', false);

end


%% Export reaction presence matrix

if isRxnPres == 1

    rxnCount = {};
    if isnan(sum(str2double(rxnPres(2,2:end))))
        for n = 2:size(rxnPres,1)
            presence = sum(cell2mat(rxnPres(n, 2:end)));
            rxnCount{n,1} = rxnPres{n,1};
            rxnCount{n,2} = presence;
            rxnCount{n,3} = 100 * presence / iters;
        end
    else
        for n = 2:size(rxnPres,1)
            presence = sum(str2double(rxnPres(n, 2:end)));
            rxnCount{n,1} = rxnPres{n,1};
            rxnCount{n,2} = presence;
            rxnCount{n,3} = 100 * presence / iters;
        end
    end

    % Export table
    rxnPres = cell2table(rxnPres);
    writetable(rxnPres, 'rxnPresenceMatrix.txt', 'Delimiter', '\t', ...
        'WriteRowNames', false, 'WriteVariableNames', false);
    rxnCount = vertcat({'Reaction', 'Occurence_Count', 'Occurence_Percent'}, rxnCount);
    rxnCount = cell2table(rxnCount);
    writetable(rxnCount, 'rxnCountTable.txt', 'Delimiter', '\t', ...
           'WriteRowNames', false, 'WriteVariableNames', false);

end
