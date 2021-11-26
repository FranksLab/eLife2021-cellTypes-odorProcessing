clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng.txt'; % set to data and catalog directory
T = readtable(Catalog, 'Delimiter', ' ');
ROIfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Spontaneous FR

for R = 1:length(ROIfiles)
    clear efd
    clear SpikeTimes
    SpikeTimes = SpikeTimes_Beast(FindFilesKK(ROIfiles{R}));
    efd = EFDmaker_Beast(ROIfiles{R},'bhv');
    
    FVon = reshape(efd.ValveTimes.FVSwitchTimesOn,[],1);
    FVall = cat(2,FVon{:});
    
    for unit = 1:length(SpikeTimes.tsec)
        st = SpikeTimes.tsec{unit}(SpikeTimes.tsec{unit}>min(FVall) & SpikeTimes.tsec{unit}<max(FVall));
        for f = 1:length(FVall)            
            Toss = st>FVall(f)-2 & st<FVall(f)+4;
            st(Toss) = [];
        end
        TossTime = 6*sum(FVall>min(FVall) & FVall<max(FVall));
        spontRate{R}(unit) = length(st)/(max(FVall)-min(FVall)-TossTime);
    end
end

%% Group LR and nLR cells spontaneous FR for each ROI

for R = 1:length(ROIfiles)
    clear LRcells
    LRcells = LRcellPicker_chgPt(ROIfiles{R},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    
    spontRate_LR{R} = spontRate{R}(LR_idx{1});
    spontRate_nLR{R} = spontRate{R}(LR_idx{2});
    LR = cat(2,spontRate_LR{:});
    nLR = cat(2,spontRate_nLR{:});
end

%% Log spont FR

logLR = log10(LR);
lognLR = log10(nLR);

%% Plot

edges = -2:.1:2;
figure; hold on
colors = {rgb('DarkGoldenrod'),rgb('ForestGreen')};

subplot(2,2,1); hold on
medLR = median(logLR);
mednLR = median(lognLR);
plot(medLR,.15,'bx'); plot(mednLR,.15,'kx');
histogram(logLR,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor',colors{1},'LineWidth',1);
histogram(lognLR,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor',rgb('Gray'),'LineWidth',1);
ax = gca; axis square
ax.YAxis.Limits = [0 .12];
ax.XTick = [min(edges) 0 max(edges)];
ax.YTick = [0 .12];
ax.XTickLabel = {'0.01', '1', '100'};
ax.XLabel.String = 'spontaneous rate (Hz)';
ax.YLabel.String = 'fraction cells';

