clear all
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
RNum = 1;
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

VOI = 2:11;
Conc = 3;
PST = [0 .5];
prePST = [-.6 -.1];
Trials = 2:11;

%% 

for k = 1:length(KWIKfiles)
    clear LRcells
    
    [efd(k)] = EFDmaker_Beast(KWIKfiles{k},'bhv');

    LRcells = LRcellPicker_chgPt(KWIKfiles{k},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    
    for lset = 1:length(LR_idx)     
        Raster{lset}{k} = efd(k).ValveSpikes.RasterAlign(VOI,Conc,LR_idx{lset});
        [~,data{lset}{k},~] = BinRearranger(Raster{lset}{k},PST,diff(PST),Trials);
        [~,dataPre{lset}{k},~] = BinRearranger(Raster{lset}{k},prePST,diff(PST),Trials);
                
        data_zsc{lset}{k} = (data{lset}{k}-mean(dataPre{lset}{k},1))./std(dataPre{lset}{k},1);
        data_zsc{lset}{k}(isinf(data_zsc{lset}{k})) = nan;
         
        rho{lset}{k} = corr(data_zsc{lset}{k}','Type','Pearson');
        rho_real{lset}{k} = rho{lset}{k} + diag(diag(nan(length(rho{lset}{k}))));
    end
end

%% Means of rho across experiments

figure;
for lset = 1:length(LR_idx)
    rho_layer{lset} = cat(3,rho_real{lset}{:});
    mean_rho{lset} = nanmean(rho_layer{lset},3);
    subplot(2,2,lset);
    imagesc(mean_rho{lset});
    colormap(pmkmp(128,'LinearL'));
    axis square; box off; axis off; ax=gca; 
%     caxis([min(min(mean_rho{lset})) max(max(mean_rho{lset}))]); 
    colorbar
%     caxis([.2 .8]);
end

%% mean within odor and mean across odor

subplot(2,2,3); hold on

for lset = 1:length(LR_idx)
    for v = 1:length(VOI)
        segment = length(Trials)*(v-1)+1:length(Trials)*v;
        meanWithin_trial{lset}{v} = nanmean(mean_rho{lset}(segment,segment));
        meanWithin{lset}(v) = nanmean(nanmean(mean_rho{lset}(segment,segment)));
        A = mean_rho{lset};
        A(segment,segment) = NaN;
        meanAcross_trial{lset}{v} = nanmean(A(segment,:)');
        meanAcross{lset}(v) = nanmean(nanmean(A(segment,:)));
    end
    
    subplot(2,2,3); hold on
    plot(lset,meanWithin{lset},'ko');
    plot(lset,meanAcross{lset},'kx');
    
    % bootstrap CI
    CI_within{lset} = bootci(1000,{@mean,meanWithin{lset}},'Type','per');
    CI_within{lset}(1) = mean(meanWithin{lset})-CI_within{lset}(1);
    CI_within{lset}(2) = CI_within{lset}(2)-mean(meanWithin{lset});
    
    CI_across{lset} = bootci(1000,{@mean,meanAcross{lset}},'Type','per');
    CI_across{lset}(1) = mean(meanAcross{lset})-CI_across{lset}(1);
    CI_across{lset}(2) = CI_across{lset}(2)-mean(meanAcross{lset});
    
    errorbar(lset+2,mean(meanWithin{lset}),CI_within{lset}(1),CI_within{lset}(2),'bx');
    errorbar(lset+2,mean(meanAcross{lset}),CI_across{lset}(1),CI_across{lset}(2),'rx');
    axis square; box off; ax=gca; ax.YAxis.Limits=[.1 .6];
    
end
