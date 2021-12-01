clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\eLife2021_DryadData\ExperimentCatalog_Ntng.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Waveform peak-trough time and spont rate

for k = 1:length(KWIKfiles)
    clear efd
    clear SpikeTimes
    clear LRcells
    SpikeTimes = SpikeTimes_Beast(FindFilesKK(KWIKfiles{k}));
    efd = EFDmaker_Beast(KWIKfiles{k},'bhv');
    
    LRcells = LRcellPicker_chgPt(KWIKfiles{k},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    
    % peak-trough time
    [ypos{k},pttime{k},asym{k},hfw{k},bigwave{k},ypos_real{k}] = WaveStats_Beast(SpikeTimes.Wave);
    
    pttime_LR{k} = pttime{k}(LR_idx{1});
    pttime_nLR{k} = pttime{k}(LR_idx{2});
    
    wave_LR{k} = bigwave{k}(LR_idx{1},:);
    wave_nLR{k} = bigwave{k}(LR_idx{2},:);
    
    asym_LR{k} = asym{k}(LR_idx{1});
    asym_nLR{k} = asym{k}(LR_idx{2});
    
    hfw_LR{k} = hfw{k}(LR_idx{1});
    hfw_nLR{k} = hfw{k}(LR_idx{2});
    
    % spontaneous FR
    FVon = reshape(efd.ValveTimes.FVSwitchTimesOn,[],1);
    FVall = cat(2,FVon{:});
    for unit = 1:length(SpikeTimes.tsec)
        for f = 1:length(FVall)
            st = SpikeTimes.tsec{unit}(SpikeTimes.tsec{unit}>min(FVall) & SpikeTimes.tsec{unit}<max(FVall));
            Toss = st>FVall(f)-2 & st<FVall(f)+4;
            st(Toss) = [];
        end
        TossTime = 6*sum(FVall>min(FVall) & FVall<max(FVall));   
        spontRate{k}(unit) = length(st)/(max(FVall)-min(FVall)-TossTime);   
    end
    spontR_LR{k} = spontRate{k}(LR_idx{1});
    spontR_nLR{k} = spontRate{k}(LR_idx{2});
end

pt = cat(1,pttime{:});
ptnLR = cat(1,pttime_nLR{:,:});
ptLR = cat(1,pttime_LR{:,:});

spontR_nLR = cat(2,spontR_nLR{:,:})';
spontR_LR = cat(2,spontR_LR{:,:})';
log_nLR = log10(spontR_nLR);
log_LR = log10(spontR_LR);

wf_LR = cat(1,wave_LR{:});
wf_nLR = cat(1,wave_nLR{:});

as = cat(1,asym{:});
as_LR = cat(1,asym_LR{:});
as_nLR = cat(1,asym_nLR{:});

hw_LR = cat(1,hfw_LR{:});
hw_nLR = cat(1,hfw_nLR{:});

%% Plot spont FR vs. pttime overlapping

figure; subplot(2,2,1); hold on;
scatter(log_nLR,ptnLR,'MarkerFaceColor','k','MarkerEdgeColor','none');
scatter(log_LR,ptLR,'MarkerFaceColor','b','MarkerEdgeColor','none');
ax = gca; axis square
ax.YAxis.Limits = [0 1];
ax.XAxis.Limits = [-1 2];

subplot(2,2,2); hold on;
scatter(log_nLR,as_nLR,'MarkerFaceColor','k','MarkerEdgeColor','none');
scatter(log_LR,as_LR,'MarkerFaceColor','b','MarkerEdgeColor','none');
ax = gca; axis square
ax.YAxis.Limits = [-1 1];
ax.XAxis.Limits = [-1 2];

%% histogram pttime

figure; subplot(2,2,1); hold on
edges = 0:.03:1;
histogram(pt,edges,'Normalization','probability','DisplayStyle','stairs','LineWidth',1);

for k = 1:length(KWIKfiles)
    fracPT{k} = sum(pttime{k} < .4)/length(pttime{k});
end

subplot(2,2,2); hold on
scatter(y,cell2mat(fracPT))
errorbar(mean(cell2mat(fracPT)),std(cell2mat(fracPT)))
box off; axis square; ylim([0 .5])

%% histogram asym

subplot(2,2,3); hold on
edges = -1:.06:1;
histogram(as,edges,'Normalization','probability','DisplayStyle','stairs','LineWidth',1);

for k = 1:length(KWIKfiles)
    fracPT{k} = sum(asym{k} > .1)/length(asym{k});
end

subplot(2,2,4); hold on
scatter(y,cell2mat(fracPT))
errorbar(mean(cell2mat(fracPT)),std(cell2mat(fracPT)))
box off; axis square; ylim([0 .5])

