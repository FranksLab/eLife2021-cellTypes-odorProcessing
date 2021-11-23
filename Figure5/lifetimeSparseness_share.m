clearvars
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

VOI = 2:11;
Conc = 3;
Cycle = 3;
TOI = {2:11};

%% Lifetime sparseness

for k = 1:length(KWIKfiles)
    clear SparseVar
    clear SpikeTimes    
    clear LRcells
    
    [efd] = EFDmaker_Beast(KWIKfiles{k},'bhv');
    Scores = SCOmaker_Beast_PreInh(KWIKfiles{k},TOI);
    spRate = Scores.RawRate(VOI,Conc,:,Cycle);
    SparseVar = squeeze(spRate);
    [SL,~] = Sparseness(spRate);
    
    LRcells = LRcellPicker_chgPt(KWIKfiles{k},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    
    SL_LR{k} = SL(LR_idx{1});
    SL_nLR{k} = SL(LR_idx{2});

end

%% Plotting lifetime sparseness histogram

edges = 0:.03:1;
figure; hold on
colors = {rgb('DarkGoldenRod')};

SL_LR_all = cat(1,SL_LR{:});
SL_nLR_all = cat(1,SL_nLR{:});

subplot(2,2,1); hold on
meanLR = nanmean(SL_LR_all);
meanNLR = nanmean(SL_nLR_all);
plot(meanLR,.1,'bx'); plot(meanNLR,.1,'kx');
histogram(SL_LR_all,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor',colors{1},'LineWidth',1);
histogram(SL_nLR_all,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',1);
ax = gca; axis square; box off;
ax.YAxis.Limits = [0 .1];
ax.XTick = [min(edges) 0.5 max(edges)];
ax.YTick = [0 .1];
ax.XTickLabel = {'0', '0.5', '1'};
ax.XLabel.String = 'lifetime sparseness';
ax.YLabel.String = ['fraction' ' ' ROI];

% bootstrap confidence interval
CI_LR = bootci(1000,{@nanmean,SL_LR_all},'Type','per');
CI_LR(1) = meanLR-CI_LR(1);
CI_LR(2) = CI_LR(2)-meanLR;
CI_NLR = bootci(1000,{@nanmean,SL_nLR_all},'Type','per');
CI_NLR(1) = meanNLR-CI_NLR(1);
CI_NLR(2) = CI_NLR(2)-meanNLR;

%% Plotting means by experiment

for k = 1:length(KWIKfiles)
    mean_LR{k} = nanmean(SL_LR{k});
    mean_nLR{k} = nanmean(SL_nLR{k});
end

CI_LR = bootci(1000,{@mean,cell2mat(mean_LR)},'Type','per');
CI_LR(1) = mean(cell2mat(mean_LR))-CI_LR(1);
CI_LR(2) = CI_LR(2)-mean(cell2mat(mean_LR));
CI_NLR = bootci(1000,{@mean,cell2mat(mean_nLR)},'Type','per');
CI_NLR(1) = mean(cell2mat(mean_nLR))-CI_NLR(1);
CI_NLR(2) = CI_NLR(2)-mean(cell2mat(mean_nLR));

subplot(2,2,3); hold on
x = 1:2;
scatter((repmat(x(1),length(mean_LR),1)),cell2mat(mean_LR),'MarkerEdgeColor','none','MarkerFaceColor',colors{1})%,'jitter','on','jitterAmount',.1)
scatter((repmat(x(2),length(mean_nLR),1)),cell2mat(mean_nLR),'MarkerEdgeColor','none','MarkerFaceColor',rgb('Gray'))%,'jitter','on','jitterAmount',.1)

errorbar(1,mean(cell2mat(mean_LR)),CI_LR(1),CI_LR(2))
errorbar(2,mean(cell2mat(mean_nLR)),CI_NLR(1),CI_NLR(2))

ylim([0 1])
set(gca,'YTick',ylim)
ylabel('lifetime sparseness')
ax = gca; box off; axis square;
ax.XTick = [1 2 3 4 5 6];
ax.XAxis.Limits = [.5 length(x)+.5];
