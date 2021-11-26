clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng.txt'; % set to data and catalog directory
T = readtable(Catalog, 'Delimiter', ' ');
ROIfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Burst Index

for R = 1:length(ROIfiles)
    clear SpikeTimes
    SpikeTimes = SpikeTimes_Beast(FindFilesKK(ROIfiles{R}));
    efd(R) = EFDmaker_Beast(ROIfiles{R},'bhv');
    
    FVon = reshape(efd(R).ValveTimes.FVSwitchTimesOn,[],1);
    FVall = cat(2,FVon{:});

    for unit = 1:length(SpikeTimes.tsec)
        st = SpikeTimes.tsec{unit}(SpikeTimes.tsec{unit}>min(FVall) & SpikeTimes.tsec{unit}<max(FVall));

        for f = 1:length(FVall)            
            Toss = st>FVall(f)-2 & st<FVall(f)+4;
            st(Toss) = [];
        end

        ISIspont{R}{unit} = diff(st);
    end
end

%%

for R = 1:length(ROIfiles)    
    
    LRcells = LRcellPicker_chgPt(ROIfiles{R},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    
    for lset = 1:length(LR_idx)
        ISIspont_celltype{lset}{R} = ISIspont{R}(LR_idx{lset});
    end
end

%%

for lset = 1:length(LR_idx)
    ISIall_cat{lset} = cat(2,ISIspont_celltype{lset}{:});
    for unit = 1:length(ISIall_cat{lset})
        medianISI{lset}(unit) = median(ISIall_cat{lset}{unit});
        meanISI{lset}(unit) = mean(ISIall_cat{lset}{unit});
        normISI{lset}(unit) = medianISI{lset}(unit)/meanISI{lset}(unit);
    end
end

figure; hold on
edges = 0:.03:1;

subplot(2,2,1); hold on
histogram(normISI{1},edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',1);
histogram(normISI{2},edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',1);
ax = gca; axis square; box off
% ax.YAxis.Limits = [0 .1];

