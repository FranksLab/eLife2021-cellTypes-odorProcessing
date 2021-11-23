clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng_Inh_Opto.txt';
T = readtable(Catalog, 'Delimiter', ' ');
ROIfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% summed spike counts for all odors all cells

VOI = 2:6;
Conc = 1;
PST = [-.6 -.1];
Trials = 1:2:20; 

for R = 1:length(ROIfiles)
    efd = EFDmaker_Beast(ROIfiles{R},'bhv');

    LRcells = LRcellPicker_chgPt(ROIfiles{R},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    for lset = 1:length(LR_idx)
        Raster = efd.ValveSpikes.RasterAlign(VOI,Conc,LR_idx{lset});
        [~, PSTHtrials, PSTHt] = PSTHmaker_Beast(Raster, PST, diff(PST), Trials);
        A = cell2mat(PSTHtrials);
        a = squeeze(A);
        b{lset}{R} = a./diff(PST);
        M = mean(b{lset}{R},2);
        allMean{R,lset} = mean(M,1);
    end
end

%%

figure;
for lset = 1:length(LR_idx)
    x{lset} = cat(2,b{lset}{:});
    xx{lset} = reshape(x{lset},[],size(x{lset},3));
    
    %bootstrap
    ci = bootci(1000,{@mean,xx{lset}},'Type','per');
    ci(1,:) = mean(xx{lset})-ci(1,:);
    ci(2,:) = ci(2,:)-mean(xx{lset});

    subplot(2,2,1); hold on
    boundedline(Trials',mean(xx{lset})',ci')
%     boundedline(Trials,mean(aM{lset}),sem(aM{lset}))
    ax=gca; ax.YAxis.Limits = [0 4];
end

    
    
    