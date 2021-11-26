clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng_Inh_Opto.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

Params.PST = [0 .3];
Params.KS = 0.01;
Params.VOI = 2:6; % valve 1 is mineral oil blank stimulus
Params.Conc = 1;
Params.Cycle = 2; % first sniff
Trials{1} = 3:2:19; % light-off trial
Trials{2} = 2:2:19; % light-on trials

%% Heatmap one column per conc, all odors

for tset = 1:length(Trials)
    [~,~,~,~,Pk{tset},Pk_nLR{tset},mx{tset},mx_nLR{tset},mxt{tset},mxt_nLR{tset}] = LatDurPeak(KWIKfiles,Params,Trials{tset});
end

%%

for k = 1:length(KWIKfiles)
    LRcells = LRcellPicker_chgPt(KWIKfiles{k},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    
    for tset = 1:length(Trials)
        Scores = SCOmaker_Beast_PreInh(KWIKfiles{k},Trials(tset));
        for lset = 1:2 %LR/NLR cells
            auroc{tset,lset}{k} = Scores.auROC(Params.VOI,Params.Conc,LR_idx{lset},Params.Cycle);
            pVal{tset,lset}{k} = Scores.AURp(Params.VOI,Params.Conc,LR_idx{lset},Params.Cycle);
        end
    end
end

%% means for each experiment

PST = Params.PST;
for tset = 1:length(Trials)
    for k = 1:length(KWIKfiles)
        mx_active{tset}{k} = mx{tset}{k}(squeeze(auroc{1,1}{k}>.5));
        mx_nLR_active{tset}{k} = mx_nLR{tset}{k}(squeeze(auroc{1,2}{k}>.5));
        mxt_active{tset}{k} = mxt{tset}{k}(squeeze(auroc{1,1}{k}>.5));
        mxt_nLR_active{tset}{k} = mxt_nLR{tset}{k}(squeeze(auroc{1,2}{k}>.5));
        
        Nopks{tset}{k} = mxt_active{tset}{k} < (PST(1)+.001) | mxt_active{tset}{k} > (PST(2)-.001);
        Nopks_nLR{tset}{k} = mxt_nLR_active{tset}{k} < (PST(1)+.001) | mxt_nLR_active{tset}{k} > (PST(2)-.001);
        mx_active{tset}{k}(Nopks{tset}{k}) = nan;
        mx_nLR_active{tset}{k}(Nopks_nLR{tset}{k}) = nan;
        mx_active{tset}{k}(isnan(mx_active{tset}{k})) = 0;
        mx_nLR_active{tset}{k}(isnan(mx_nLR_active{tset}{k})) = 0;        
    end
end

for tset = 1:2
    for k = 1:length(KWIKfiles)
        mean_pk{tset}{k} = mean(mx_active{tset}{k});
        mean_pk_nLR{tset}{k} = mean(mx_nLR_active{tset}{k});        
    end
end

%% 

figure; 

subplot(2,2,1); hold on
x = 1:4;
scatter((repmat(x(1),length(KWIKfiles),1)),cell2mat(mean_pk{1}),'MarkerEdgeColor','none','MarkerFaceColor','r');%,'jitter','on','jitterAmount',.1)
scatter((repmat(x(2),length(KWIKfiles),1)),cell2mat(mean_pk{2}),'MarkerEdgeColor','none','MarkerFaceColor',rgb('Gray'));%,'jitter','on','jitterAmount',.1)

scatter((repmat(x(3),length(KWIKfiles),1)),cell2mat(mean_pk_nLR{1}),'MarkerEdgeColor','none','MarkerFaceColor','r');%,'jitter','on','jitterAmount',.1)
scatter((repmat(x(4),length(KWIKfiles),1)),cell2mat(mean_pk_nLR{2}),'MarkerEdgeColor','none','MarkerFaceColor',rgb('Gray'));%,'jitter','on','jitterAmount',.1)

% bootstrap
for tset = 1:2
    ci_LR{tset} = bootci(1000,{@mean,cell2mat(mean_pk{tset})},'Type','per');
    ci_LR{tset}(1) = mean(cell2mat(mean_pk{tset}))-ci_LR{tset}(1);
    ci_LR{tset}(2) = ci_LR{tset}(2)-mean(cell2mat(mean_pk{tset}));
    ci_NLR{tset} = bootci(1000,{@mean,cell2mat(mean_pk_nLR{tset})},'Type','per');
    ci_NLR{tset}(1) = mean(cell2mat(mean_pk_nLR{tset}))-ci_NLR{tset}(1);
    ci_NLR{tset}(2) = ci_NLR{tset}(2)-mean(cell2mat(mean_pk_nLR{tset}));
end

errorbar(1,mean(cell2mat(mean_pk{1})),ci_LR{1}(1),ci_LR{1}(2),'kx')
errorbar(2,mean(cell2mat(mean_pk{2})),ci_LR{2}(1),ci_LR{2}(2),'kx')
errorbar(3,mean(cell2mat(mean_pk_nLR{1})),ci_NLR{1}(1),ci_NLR{1}(2),'kx')
errorbar(4,mean(cell2mat(mean_pk_nLR{2})),ci_NLR{2}(1),ci_NLR{2}(2),'kx')

ylim([0 25])
box off; axis square;

for k = 1:length(KWIKfiles)
    plot([1 2],[mean_pk{1}{k} mean_pk{2}{k}])
    plot([3 4],[mean_pk_nLR{1}{k} mean_pk_nLR{2}{k}])
end

%% all cell-odor pairs combined

for tset = 1:length(Trials)
    for lset = 1:2 %LR/NLR cells
        aur{tset,lset} = permute(cat(3,auroc{tset,lset}{:}),[3,1,2]);
        p{tset,lset} = permute(cat(3,pVal{tset,lset}{:}),[3,1,2]);
        aur_string{tset,lset} = reshape(aur{tset,lset}',[],1);
        p_string{tset,lset} = reshape(p{tset,lset}',[],1);  
    end
end

for tset = 1:length(Trials)
    pk_active_LR{tset} = Pk{tset}(aur_string{1,1}>.5 & p_string{1,1}<.05);
    pk_active_nLR{tset} = Pk_nLR{tset}(aur_string{1,2}>.5 & p_string{1,2}<.05);
    pk_active_LR{tset}(isnan(pk_active_LR{tset})) = 0;
    pk_active_nLR{tset}(isnan(pk_active_nLR{tset})) = 0;
end

%%

figure; hold on
PSTaxis = [0 2.5];

subplot(2,2,1); hold on
plot(log10(pk_active_LR{1}+1),log10(pk_active_LR{2}+1),'k.');
ax = gca; ax.XAxis.Limits = PSTaxis; ax.YAxis.Limits = PSTaxis; box off; axis square;
plot(PSTaxis, PSTaxis,'r--');
plot(nanmedian(log10(pk_active_LR{1}+1)),nanmedian(log10(pk_active_LR{2}+1)),'rx');

subplot(2,2,2); hold on
plot(log10(pk_active_nLR{1}+1),log10(pk_active_nLR{2}+1),'k.');
ax = gca; ax.XAxis.Limits = PSTaxis; ax.YAxis.Limits = PSTaxis; box off; axis square;
plot(PSTaxis, PSTaxis,'r--');
plot(nanmedian(log10(pk_active_nLR{1}+1)),nanmedian(log10(pk_active_nLR{2}+1)),'rx');
