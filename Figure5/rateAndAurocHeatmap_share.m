clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng.txt';
T = readtable(Catalog, 'Delimiter', ' ');
ROIfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

VOI = 2:11;
Conc = 3;
TOI{1} = 2:11;
Cycle = 3;

%% Heatmap one column per conc, all odors

for R = 1:length(ROIfiles)
    for lset = 1:2 %LR/NLR cells
        clear Scores
        clear LRcells
                
        Scores = SCOmaker_Beast_PreInh(ROIfiles{R},TOI);
        
        LRcells = LRcellPicker_chgPt(ROIfiles{R},[-.1 .1]);
        LR_idx{1} = LRcells.primLR;
        LR_idx{2} = LRcells.nonLR;
        
        rawRate{lset,R} = Scores.RawRate(VOI,Conc,LR_idx{lset},Cycle);
        auroc{lset,R} = Scores.auROC(VOI,Conc,LR_idx{lset},Cycle)*2-1;
    end
end

for lset = 1:2 %LR/NLR cells
    HeatVar{lset} = permute(cat(3,auroc{lset,:}),[3,1,2]);
    RR{lset} = permute(cat(3,rawRate{lset,:}),[3,1,2]);
end

%% Unsorted rate heatmaps

figure; hold on
for lset = 1:2
    subplot(1,2,lset); axis on; hold on
    imagesc(RR{lset});
    HT = hot(32);
    HT = HT(1:end-4,:);
    colormap(HT)
    colorbar('southoutside')
    caxis([0 15])
end

%% Sorted auroc heatmaps

figure; hold on

for lset = 1:2
    clear cts; clear ups; clear downs;
    for r = 1:size(HeatVar{lset},1)
        ups(r) = sum(HeatVar{lset}(r,:)>0);
        downs(r) = sum(HeatVar{lset}(r,:)<0);
        cts = ups - downs;
    end
    [B,I{lset}] = sort(cts);
    Heatsort{lset} = HeatVar{lset}(I{lset},:);
    subplot(1,2,lset); hold on
    imagesc(Heatsort{lset});
    CT = flipud(cbrewer('div','RdBu',64));
    CT = CT(3:end-3,:);
    colormap(CT)
    colorbar('southoutside')
    caxis([-1 1])
end

%% Sorted rate heatmaps

clear cts; clear ups; clear downs; 

figure; hold on
for lset = 1:2
    Ratesort{lset} = RR{lset}(I{lset},:);
    subplot(1,2,lset); hold on
    imagesc(Ratesort{lset});
    HT = hot(32);
    HT = HT(1:end-4,:);
    colormap(HT)
    colorbar('southoutside')
    caxis([0 15])
end