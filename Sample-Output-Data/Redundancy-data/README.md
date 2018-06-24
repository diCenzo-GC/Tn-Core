# Redundancy data

The files listed below are included in this directory.

geneCoocurrenceAdjustedMatrix.txt: A matrix containing a Chi-squared statistic for the likelihood that the genes on the X and Y axes are more (positive) or less (negative) likely to occur in the same core model than by chance.

geneCoocurrenceMatrix.txt: The number of core models that contain both the gene on the X axis and the gene one the Y axis. 

geneCoocurrenceTable.txt: The first two columns represent each pair of genes for the genes variably present/absent in the core models. The third and fourth columns indicate the number of core models that contain the genes. The fifth column indicates the number of core models containing both of the genes, while the sixth column indicates the number of core models expected to have both genes by chance based on their individual occurrence. The seventh column is the ratio of the observed cooccurrence (numerator) to the expected concurrence (denominator). The eighth column contains a Chi-squared statistic for the likelihood that the genes are more (positive) or less (negative) likely to occur in the same core model than by chance.

geneCountTable.txt: A table containing all genes found in at least one core model. The first column is the gene name. The second column indicates the number of core models containing the gene. The third column indicates the percentage of core models containing the gene.

genePresenceMatrix.txt: The fist column contains the names of all genes found in at least one core model. The rest of the columns represent all of the generated core models. For each row (i.e., for each gene), the presence (1) or absence (0) of the gene in the model is indicated.

geneCountThresholdTable.txt: Similar to the geneCountTable.txt file, except that it: i) contains all genes in the reduced model (not only those found in at least one core model), and ii) contains gene counts for multiple core model populations, each produced using different growth threshold.

rxnCoocurrenceAdjustedMatrix.txt: A matrix containing a Chi-squared statistic for the likelihood that the reactions on the X and Y axes are more (positive) or less (negative) likely to occur in the same core model than by chance.

rxnCoocurrenceMatrix.txt: The number of core models that contain both the reaction on the X axis and the reaction one the Y axis. 

rxnCoocurrenceTable.txt: The first two columns represent each pair of reactions for the reactions variably present/absent in the core models. The third and fourth columns indicate the number of core models that contain the reactions. The fifth column indicates the number of core models containing both of the reactions, while the sixth column indicates the number of core models expected to have both reactions by chance based on their individual occurrence. The seventh column is the ratio of the observed cooccurrence (numerator) to the expected concurrence (denominator). The eighth column contains a Chi-squared statistic for the likelihood that the reactions are more (positive) or less (negative) likely to occur in the same core model than by chance.

rxnCountTable.txt: A table containing all reactions found in at least one core model. The first column is the reaction name. The second column indicates the number of core models containing the reaction. The third column indicates the percentage of core models containing the reaction.

rxnPresenceMatrix.txt: The fist column contains the names of all reactions found in at least one core model. The rest of the columns represent all of the generated core models. For each row (i.e., for each reaction), the presence (1) or absence (0) of the reaction in the model is indicated.

rxnCountThresholdTable.txt: Similar to the rxnCountTable.txt file, except that it: i) contains all reactions in the reduced model (not only those found in at least one core model), and ii) contains reaction counts for multiple core model populations, each produced using different growth threshold.
