clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng_Inh_Opto.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

VOI = 2:6;
Conc = 1;
TrialSets{1} = 3:2:19;
TrialSets{2} = 2:2:19;
Cycle = 2; % first sniff

%% Heatmap one column per conc, all odors

for k = 1:length(KWIKfiles)
    for tset = 1:length(TrialSets)
        for lset = 1:2 %LR/NLR cells
            clear LRcells
            Scores{k} = SCOmaker_Beast_PreInh(KWIKfiles{k},TrialSets(tset));
            
            LRcells = LRcellPicker_chgPt(KWIKfiles{k},[-.1 .1]);
            LR_idx{1} = LRcells.primLR;
            LR_idx{2} = LRcells.nonLR;
            
            rawRate{tset,lset}{k} = Scores{k}.RawRate(VOI,Conc,LR_idx{lset},Cycle);
            
        end
    end
end

%% concatenate all experiments

for tset = 1:length(TrialSets)
    for lset = 1:2 %LR/NLR cells
        rr{tset,lset} = squeeze(cat(3,rawRate{tset,lset}{:}));
    end
end                                                 

%% Rank correlation

for lset = 1:2
    rho{lset} = corr(rr{1,lset}',rr{2,lset}','Type','Pearson');
    for v = 1:length(VOI)
        diag{lset}(v) = rho{lset}(v,v);
    end
end

% bootstrap CI
for lset = 1:2
    CI{lset} = bootci(1000,{@mean,diag{lset}},'Type','per');
    CI{lset}(1) = mean(diag{lset})-CI{lset}(1);
    CI{lset}(2) = CI{lset}(2)-mean(diag{lset});
end

%%

figure; 

subplot(2,2,1); hold on
for lset = 1:2  
    x = 1:2;    
    scatter(repmat(x(lset),length(diag{lset}),1)',diag{lset},'jitter','on','jitterAmount',0.1,'MarkerFaceColor','k','MarkerEdgeColor','none')
    ax = gca; box off; axis square; 
    ax.YAxis.Limits = [.2 1];
    ax.XAxis.Limits = [0 3];
    set(gca,'XTick',[1 2])
    ylabel('corr. coef.')
    hold on
    errorbar(lset,mean(diag{lset}),CI{lset}(1),CI{lset}(2),'rx')
end
