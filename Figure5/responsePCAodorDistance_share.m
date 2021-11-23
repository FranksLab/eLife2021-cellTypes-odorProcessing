clearvars
% close all
clc

%% Pick out files with 'kwik' in its name and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

VOI = 2:11;
PST = [0 .5];
prePST = [-.6 -.1];
BinSize = diff(PST);
pcaTrials = 2:11;
Conc = 3;

%% Rearrange spike counts trial by trial

for k = 1:length(KWIKfiles)
    efd = EFDmaker_Beast(KWIKfiles{k},'bhv');

    LRcells = LRcellPicker_chgPt(KWIKfiles{k},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    
    for lset = 1:2
        Raster{lset}{k} = efd.ValveSpikes.RasterAlign(VOI,Conc,LR_idx{lset});
        [label,pcadata{lset}{k},~] = BinRearranger(Raster{lset}{k},PST,BinSize,pcaTrials);
        [~,pcadataPre{lset}{k},~] = BinRearranger(Raster{lset}{k},prePST,BinSize,pcaTrials);
        pcadata_zsc{lset}{k} = (pcadata{lset}{k}-mean(pcadataPre{lset}{k},1))./std(pcadataPre{lset}{k},1);
        pcadata_zsc{lset}{k}(isinf(pcadata_zsc{lset}{k})) = nan;
    end
end

%% PCA

for lset = 1:2
    allData{lset} = cell2mat(pcadata_zsc{lset});
    [~,mappedX{lset},~,~,explained{lset}] = pca(allData{lset},'Algorithm','svd','NumComponents',10);
end 

%% color codes

CT{1}=rgb('LightSlateGray');
CT{2}=rgb('FireBrick');
CT{3}=rgb('PaleVioletRed');
CT{4}=rgb('DarkOrange');
CT{5}=rgb('DarkGreen');
CT{6}=rgb('BurlyWood');
CT{7}=rgb('SaddleBrown');
CT{8}=rgb('Gold');
CT{9}=rgb('SteelBlue');
CT{10}=rgb('Indigo');

%% plot on xyz axes

sz = 30;

for lset = 1:2
    figure;
    for o = 1:length(VOI)
        odorTrials = length(pcaTrials)*(o-1)+1:length(pcaTrials)*o;
        scatter3(mappedX{lset}(odorTrials,1),mappedX{lset}(odorTrials,2),mappedX{lset}(odorTrials,3),sz,CT{o},'filled')
%         line(mappedX{lset}(odorTrials,1),mappedX{lset}(odorTrials,2),mappedX{lset}(odorTrials,3))
        hold on
        [x,y,z] = ellipsoid(mean(mappedX{lset}(odorTrials,1)),mean(mappedX{lset}(odorTrials,2)),...
            mean(mappedX{lset}(odorTrials,3)),std(mappedX{lset}(odorTrials,1)),...
            std(mappedX{lset}(odorTrials,2)),std(mappedX{lset}(odorTrials,3)));
        surf(x,y,z,'FaceAlpha',0.5,'EdgeColor','none','FaceColor',CT{o},'LineStyle','--');
    end
    ax=gca; box off; axis square;
    ax.XAxis.Limits=[-20 80]; ax.YAxis.Limits=[-20 60]; ax.ZAxis.Limits=[-30 20];
%     view(-40,10)
end

%% Distance between mean of all trials

figure; hold on
subplot(2,2,1); hold on
for lset = 1:2
    clear mean_coord
    for o = 1:length(VOI)
        odorTrials = length(pcaTrials)*(o-1)+1:length(pcaTrials)*o;
        mean_coord{o} = [mean(mappedX{lset}(odorTrials,1)),mean(mappedX{lset}(odorTrials,2)),...
            mean(mappedX{lset}(odorTrials,3)),mean(mappedX{lset}(odorTrials,4)),mean(mappedX{lset}(odorTrials,5)),mean(mappedX{lset}(odorTrials,6))];
        std_coord{lset}{o} = [std(mappedX{lset}(odorTrials,1)),std(mappedX{lset}(odorTrials,2)),...
            std(mappedX{lset}(odorTrials,3)),std(mappedX{lset}(odorTrials,4)),std(mappedX{lset}(odorTrials,5)),std(mappedX{lset}(odorTrials,6))];
        std_odor{lset}{o} = mean(std_coord{lset}{o});
    end
    dist_odors{lset} = pdist(cell2mat(mean_coord'),'euclidean'); 
    bar(lset,mean(dist_odors{lset}),'BarWidth',.5);
    scatter(repmat(lset,1,length(dist_odors{lset})),dist_odors{lset})
    errorbar(mean(dist_odors{lset}),sem(dist_odors{lset}'))
    ax = gca; box off; ax.YAxis.Limits = [0 50]; axis square;
end

edges = 0:2:50;
subplot(2,2,2); hold on
histogram(dist_odors{1},edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','r')
histogram(dist_odors{2},edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k')
ax = gca; box off; ax.YAxis.Limits = [0 .3]; ax.XAxis.Limits = [0 50]; axis square;

%% distance confusion matrix

figure;
for lset = 1:2
    Z{lset} = squareform(dist_odors{lset});
    Z_real{lset} = Z{lset} + diag(diag(nan(length(Z{lset}))));
    subplot(2,2,lset); hold on; 
    imagesc(Z{lset})
    colormap(pmkmp(128,'LinearL'));
    caxis([0 50])
    axis square
end
