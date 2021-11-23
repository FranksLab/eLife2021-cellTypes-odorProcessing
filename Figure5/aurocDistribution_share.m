clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

VOI = 2:11;
Conc = 3;
TrialSets{1} = 2:11;
Cycle = 3;

%% Heatmap one column per conc, all odors

for k = 1:length(KWIKfiles)
    for lset = 1:2 %LR/NLR cells
        clear Scores
        clear LRcells
        Scores = SCOmaker_Beast_PreInh(KWIKfiles{k},TrialSets);
        LRcells = LRcellPicker_chgPt(KWIKfiles{k},[-.1 .1]);
        LR_idx{1} = LRcells.primLR;
        LR_idx{2} = LRcells.nonLR;
        
        auroc{lset}{k} = Scores.auROC(VOI,Conc,LR_idx{lset},Cycle)*2-1;
    end
end

%% Reshape

for lset = 1:2 %LR/NLR cells
    aur{lset} = reshape(squeeze(cat(3,auroc{lset}{:})),[],1);
end
    
%% Plotting

figure; subplot(2,2,1); hold on
edges = -1:.08:1;
histogram(aur{1},edges,'Normalization','probability','DisplayStyle','stairs')
histogram(aur{2},edges,'Normalization','probability','DisplayStyle','stairs')
ax=gca; box off; axis square; ax.XAxis.Limits = [-1 1];
ax.XLabel.String = 'response index';
ax.YLabel.String = 'fraction responses';

