%% Set up workspave

clear all
addpath(genpath('../Software/cobratoolbox/'));
rmpath(genpath('../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../Software/tiger/'));
addpath(genpath('../Software/Tn-Core-v1.2/'));
addpath(genpath('../Software/FastCore/'));
initCobraToolbox;
changeCobraSolver('gurobi', 'all');
rmpath(genpath('../Software/FastCore/'));
addpath(genpath('../Software/FastCore/'));

%% Perform the S. meliloti refinement

% Get the iGD1575 model
load('iGD1575.mat');

% Get the draft meliloti model
load('draftMeliloti.mat');

% Get tnseq data
tnseq_meliloti = table2cell(readtable('TnSeqData.txt', 'ReadVariableNames', false, ...
    'ReadRowNames', false, 'Delimiter', '\t'));

% Refine iGD1575
[refinedModelA, reportA] = tncore_refine(iGD1575, tnseq_meliloti);

% Refine draft meliloti
[refinedModelB, reportB] = tncore_refine(draftMeliloti, tnseq_meliloti);

%% Perform the P. aeruginosa refinement

% Get the iPae1146 model
load('iPAE1146.mat');

% Get tnseq data
load('inputVariables_iPAE1146.mat');
tnseq_aeruginosa = tnseq;
clearvars tnseq

% Refine iPae1146
[refinedModelC, reportC] = tncore_refine(iPAE1146, tnseq_aeruginosa);

%% Perform the R. sphaeroides refinement for aero

% Get the iRsp1140 model
iRsp1140 = readCbModel('iRsp1140.xml');

% Get tnseq data
tnseq_sphaeroides = table2cell(readtable('TnSeqData_iRsp1140.txt', 'ReadVariableNames', false, ...
    'ReadRowNames', false, 'Delimiter', '\t'));

% Set objective
iRsp1140 = changeObjective(iRsp1140, 'RXN1391');

% Set medium
EX_list = iRsp1140.rxns(findExcRxns(iRsp1140));
iRsp1140.rev(findExcRxns(iRsp1140)) = 1;
iRsp_medium = {'RXN0217';'RXN0223';'RXN0213';'RXN0222';'RXN0219';'RXN0191';...
    'RXN0190';'RXN0224';'RXN0196';'RXN0188';'RXN0221';'RXN1395';'RXN1326';...
    'RXN0193';'RXN1167';'RXN1329';'RXN1158';'RXN1354';'RXN0220';'RXN0192'};
iRsp1140 = changeRxnBounds(iRsp1140, EX_list, 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, EX_list, 1000, 'u');
iRsp1140 = changeRxnBounds(iRsp1140, iRsp_medium, -100, 'l');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0224', -5, 'l');

% Set additional model constraints
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0109', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0631', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN1827', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN1843', 0, 'b');

% Refine iRsp1140
[refinedModelD, reportD] = tncore_refine(iRsp1140, tnseq_sphaeroides);

%% Compare iRsp1140 automated refinement with manual refinement for aero

% Get manually refined iRsp1140
iRsp1140_opt = readCbModel('iRsp1140_opt');

% Identify reactions manually refined
manuallyRefinedListAero = {};
for n = 1:length(iRsp1140.rxns)
    pos = strmatch(iRsp1140.rxns{n}, iRsp1140_opt.rxns, 'exact');
    if ~isempty(pos)
        if strmatch(iRsp1140.grRules{n}, iRsp1140_opt.grRules{pos}, 'exact')
        else
            manuallyRefinedListAero = vertcat(manuallyRefinedListAero, iRsp1140.rxns{n});
        end
    end
end

% Identify reactions automatedly refined
automatedRefinedListAero = {};
for n = 1:length(refinedModelD.rxns)
    pos = strmatch(refinedModelD.rxns{n}, iRsp1140.rxns, 'exact');
    if ~isempty(pos)
        if strmatch(refinedModelD.grRules{n}, iRsp1140.grRules{pos}, 'exact')
        else
            automatedRefinedListAero = vertcat(automatedRefinedListAero, refinedModelD.rxns{n});
        end
    end
end

% Overlap between automated and manual refinement
overlapRefinedListsAero = intersect(manuallyRefinedListAero, automatedRefinedListAero);

% Reactions only manually refined
onlyManuallyRefinedAero = setdiff(manuallyRefinedListAero, automatedRefinedListAero);

% Reactions only automatedly refined
onlyAtomatedRefinedAero = setdiff(automatedRefinedListAero, manuallyRefinedListAero);

%% Perform the R. sphaeroides refinement for photo

% Get the iRsp1140 model
iRsp1140 = readCbModel('iRsp1140.xml');

% Get tnseq data
tnseq_sphaeroidesPhoto = table2cell(readtable('TnSeqData_iRsp1140_photo.txt', 'ReadVariableNames', false, ...
    'ReadRowNames', false, 'Delimiter', '\t'));

% Set objective
iRsp1140 = changeObjective(iRsp1140, 'RXN1306');

% Set medium
EX_list = iRsp1140.rxns(findExcRxns(iRsp1140));
iRsp1140.rev(findExcRxns(iRsp1140)) = 1;
iRsp_medium = {'RXN0217';'RXN0223';'RXN0213';'RXN0222';'RXN0222';'RXN0219';'RXN0191';...
    'RXN0190';'RXN0196';'RXN0188';'RXN0221';'RXN1395';'RXN1326';'RXN1205';...
    'RXN0193';'RXN1167';'RXN1329';'RXN1158';'RXN1354';'RXN0220';'RXN0192'; 'RXN1331'};
iRsp1140 = changeRxnBounds(iRsp1140, EX_list, 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, EX_list, 1000, 'u');
iRsp1140 = changeRxnBounds(iRsp1140, iRsp_medium, -100, 'l');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0219', -2.5, 'l');

% Set additional model constraints
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0949', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN1318', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0661', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0630', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0663', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN1063', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN1064', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0515', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN1020', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN1345', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0648', 0, 'b');
iRsp1140 = changeRxnBounds(iRsp1140, 'RXN0648', 1000, 'u');

% Refine iRsp1140
[refinedModelE, reportE] = tncore_refine(iRsp1140, tnseq_sphaeroidesPhoto);

%% Compare iRsp1140 automated refinement with manual refinement for photo

% Get manually refined iRsp1140
iRsp1140_opt = readCbModel('iRsp1140_opt');

% Identify reactions manually refined
manuallyRefinedListPhoto = {};
for n = 1:length(iRsp1140.rxns)
    pos = strmatch(iRsp1140.rxns{n}, iRsp1140_opt.rxns, 'exact');
    if ~isempty(pos)
        if strmatch(iRsp1140.grRules{n}, iRsp1140_opt.grRules{pos}, 'exact')
        else
            manuallyRefinedListPhoto = vertcat(manuallyRefinedListPhoto, iRsp1140.rxns{n});
        end
    end
end

% Identify reactions automatedly refined
automatedRefinedListPhoto = {};
for n = 1:length(refinedModelE.rxns)
    pos = strmatch(refinedModelE.rxns{n}, iRsp1140.rxns, 'exact');
    if ~isempty(pos)
        if strmatch(refinedModelE.grRules{n}, iRsp1140.grRules{pos}, 'exact')
        else
            automatedRefinedListPhoto = vertcat(automatedRefinedListPhoto, refinedModelE.rxns{n});
        end
    end
end

% Overlap between automated and manual refinement
overlapRefinedListsPhoto = intersect(manuallyRefinedListPhoto, automatedRefinedListPhoto);

% Reactions only manually refined
onlyManuallyRefinedPhoto = setdiff(manuallyRefinedListPhoto, automatedRefinedListPhoto);

% Reactions only automatedly refined
onlyAtomatedRefinedPhoto = setdiff(automatedRefinedListPhoto, manuallyRefinedListPhoto);

% Total overlap between automated and manual refinement
totalAuto = vertcat(automatedRefinedListPhoto, automatedRefinedListAero);
overlapAll = intersect(manuallyRefinedListPhoto, totalAuto);

%% Save and clear

save('refinement_validation_output.mat');
clear

