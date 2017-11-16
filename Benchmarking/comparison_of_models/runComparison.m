%% Initialize matlab

clear;

%% Prepare original models

prepareExistingModels;

%% Load models

load('iGD1575.mat');
load('iGD726.mat');
load('coreModelA.mat');
load('coreModelB.mat');
load('coreModelC.mat');
load('coreModelD.mat');
load('coreModel_fastcore.mat');
load('coreModel_gimme.mat');
load('coreModel_minnw.mat');

%% Compare gene essentialities

% Single gene deletion analysis
iGD1575_grRatio = singleGeneDeletion(iGD1575, 'MOMA');
iGD726_grRatio = singleGeneDeletion(iGD726, 'MOMA');
coreModelA_grRatio = singleGeneDeletion(coreModelA, 'MOMA');
coreModelB_grRatio = singleGeneDeletion(coreModelB, 'MOMA');
coreModelC_grRatio = singleGeneDeletion(coreModelC, 'MOMA');
coreModelD_grRatio = singleGeneDeletion(coreModelD, 'MOMA');
coreModel_fastcore_grRatio = singleGeneDeletion(coreModel_fastcore, 'MOMA');
coreModel_gimme_grRatio = singleGeneDeletion(coreModel_gimme, 'MOMA');
coreModel_minnw_grRatio = singleGeneDeletion(coreModel_minnw, 'MOMA');

% List of genes of interest
carbonGenes = {'smc03070'; 'smc03069'; 'smc03153'; 'smc04262'; 'smc00152'; ...
    'smc03978'; 'smc02495'; 'smc03979'; 'smc03981'; 'smc01028'; 'smc04005'; ...
    'smc01030'; 'smc02087'; 'smc03846'; 'smc00480'; 'smc02482'; 'smc02465'; ...
    'smc00149'; 'smc02479'; 'smc03895'; 'smc02500'};

% Make output variable
outputCompared = cell(length(carbonGenes), 9);

% Pull out gene numbers
iGD1575_IDs = findGeneIDs(iGD1575, carbonGenes);
iGD726_IDs = findGeneIDs(iGD726, carbonGenes);
coreModelA_IDs = findGeneIDs(coreModelA, carbonGenes);
coreModelB_IDs = findGeneIDs(coreModelB, carbonGenes);
coreModelC_IDs = findGeneIDs(coreModelC, carbonGenes);
coreModelD_IDs = findGeneIDs(coreModelD, carbonGenes);
coreModel_fastcore_IDs = findGeneIDs(coreModel_fastcore, carbonGenes);
coreModel_gimme_IDs = findGeneIDs(coreModel_gimme, carbonGenes);
coreModel_minnw_IDs = findGeneIDs(coreModel_minnw, carbonGenes);

% Get the data
for n = 1:length(carbonGenes)
    
    if iGD1575_IDs(n) ~= 0
        outputCompared{n,1} = round(iGD1575_grRatio(iGD1575_IDs(n)),3);
    else
        outputCompared{n,1} = -1;
    end
    
    if iGD726_IDs(n) ~= 0
        outputCompared{n,2} = round(iGD726_grRatio(iGD726_IDs(n)),3);
    else
        outputCompared{n,2} = -1;
    end
    
    if coreModelA_IDs(n) ~= 0
        outputCompared{n,3} = round(coreModelA_grRatio(coreModelA_IDs(n)),3);
    else
        outputCompared{n,3} = -1;
    end
    
    if coreModelB_IDs(n) ~= 0
        outputCompared{n,4} = round(coreModelB_grRatio(coreModelB_IDs(n)),3);
    else
        outputCompared{n,4} = -1;
    end     
    
    if coreModelC_IDs(n) ~= 0
        outputCompared{n,5} = round(coreModelC_grRatio(coreModelC_IDs(n)),3);
    else
        outputCompared{n,5} = -1;
    end     
    
    if coreModelD_IDs(n) ~= 0
        outputCompared{n,6} = round(coreModelD_grRatio(coreModelD_IDs(n)),3);
    else
        outputCompared{n,6} = -1;
    end     

    if coreModel_fastcore_IDs(n) ~= 0
        outputCompared{n,7} = round(coreModel_fastcore_grRatio(coreModel_fastcore_IDs(n)),3);
    else
        outputCompared{n,7} = -1;
    end   
    
    if coreModel_gimme_IDs(n) ~= 0
        outputCompared{n,9} = round(coreModel_gimme_grRatio(coreModel_gimme_IDs(n)),3);
    else
        outputCompared{n,9} = -1;
    end     
    
    if coreModel_minnw_IDs(n) ~= 0
        outputCompared{n,8} = round(coreModel_minnw_grRatio(coreModel_minnw_IDs(n)),3);
    else
        outputCompared{n,8} = -1;
    end   
    
end

% Add labels
headers = {'Gene', 'iGD1575', 'iGD726', 'coreModelA', 'coreModelB', ...
    'coreModelC', 'coreModelD', 'fastcore' , 'minNW', 'gimme'};
outputCompared = horzcat(carbonGenes, outputCompared);
outputCompared = vertcat(headers, outputCompared);
outputCompared = transpose(outputCompared);

% Save the data
save('allComparisonData.mat');
save('outputComparison.mat', 'outputCompared');

% Export as Excell table
toExport = cell2table(outputCompared);
writetable(toExport, 'carbonGenesCompared.xlsx', 'WriteVariableNames', false);

% Clear
clear;




