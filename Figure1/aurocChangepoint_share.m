clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng.txt'; % set to data directory and catalog
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%%

PST = [-.1 .1];

for k = 1:length(KWIKfiles)
    [efd] = EFDmaker_Beast(KWIKfiles{k},'bhv');
    % auroc and ranksum p-value
    for unit = 1:size(efd.LaserSpikes.SpikesBeforeLaser,2)
        Control = efd.LaserSpikes.SpikesBeforeLaser{unit};
        Stimulus = efd.LaserSpikes.SpikesDuringLaser{unit};
        [auROC{k}(unit), p{k}(unit)] = RankSumROC(Control, Stimulus); % ranksum select laser responsive cells
        LR_p{k}(unit) =  auROC{k}(unit) < .5 & p{k}(unit) < .0000001;
    end
    % changepoint
    numTrials = length(efd.LaserSpikes.RasterAlign{1});    
    for unit = 1:size(efd.LaserSpikes.RasterAlign,2)
        [PSTH, ~, PSTHt] = PSTHmaker_Beast(efd.LaserSpikes.RasterAlign(unit), PST, 0.002, 1:numTrials);
        ipt = findchangepts(PSTH{:});
        changept{k}(unit) = PSTHt(ipt);
        LR_lat{k}(unit) = changept{k}(unit) > -.05 & changept{k}(unit) <= .01;
    end
   LR_idx{k} = find(LR_p{k} == 1 & LR_lat{k} == 1);     
end

%% plot auroc vs changept

figure; subplot(2,2,1); hold on
plot(cell2mat(changept),cell2mat(auROC),'ko',...
    'MarkerEdgeColor','none','MarkerFaceColor','k')
ax = gca; box off; ax.XAxis.Limits = [-.1 .1]; ax.YAxis.Limits = [0 1]; axis square
ax.XLabel.String = 'changepoint (s)';
ax.YLabel.String = 'auroc';
hold on
rectangle('Position',[-.05 0 .01 .5],'EdgeColor','r')

%% mark opto-tagged cells as red circles

for k = 1:length(KWIKfiles)
    plot(changept{k}(LR_idx{k}),auROC{k}(LR_idx{k}),'ro',...
    'MarkerEdgeColor','none','MarkerFaceColor','r')
end
    
    
