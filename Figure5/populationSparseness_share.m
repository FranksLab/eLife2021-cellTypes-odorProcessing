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
Cycle = 3;
TOI = {2:11};

%% Population sparseness

for R = 1:length(ROIfiles)
    clear SparseVar
    clear LRcells
    [efd{R}] = EFDmaker_Beast(ROIfiles{R},'bhv');
    Scores{R} = SCOmaker_Beast_PreInh(ROIfiles{R},TOI);
    
    LRcells = LRcellPicker_chgPt(ROIfiles{R},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;    
    LR_idx{2} = LRcells.nonLR;
   
    for lset = 1:length(LR_idx)
        spRate{lset} = Scores{R}.RawRate(VOI,Conc,LR_idx{lset},Cycle);
        SparseVar{lset} = squeeze(spRate{lset});
        [~,SP{lset}{R}] = Sparseness(SparseVar{lset});
    end
end

%% Means for each experiment

for lset = 1:length(LR_idx)
    popSp{lset} = cat(1,SP{lset}{:});
    meanPopSp{lset} = mean(popSp{lset},1);
end

% bootstrap CI
for lset = 1:length(LR_idx)
    CI{lset} = bootci(1000,{@mean,meanPopSp{lset}},'Type','per');
    CI{lset}(1) = mean(meanPopSp{lset})-CI{lset}(1);
    CI{lset}(2) = CI{lset}(2)-mean(meanPopSp{lset});
end

%% Plotting

figure; subplot(2,2,1); hold on
for lset = 1:length(LR_idx)
    x = repmat(lset,1,size(meanPopSp{1},2));
    scatter(x,meanPopSp{lset},'MarkerFaceColor','k','MarkerEdgeColor','none')%,'jitter','on','jitterAmount',0.2)
    errorbar(lset,mean(meanPopSp{lset}),CI{lset}(1),CI{lset}(2),'rx');
end

ylim([0.5 1])
xlim([0 3])
set(gca,'XTick',[0 1 2 3])
set(gca,'YTick',[0 0.5 1])
ylabel('population sparseness')
box off; axis square

for v = 1:length(VOI)
    plot([1 2],[meanPopSp{1}(v) meanPopSp{2}(v)])
end
