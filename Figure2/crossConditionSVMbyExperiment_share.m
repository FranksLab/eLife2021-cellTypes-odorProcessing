clearvars
% close all
clc

%% Pick out files with 'kwik' in its name and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng_Inh_Opto.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

VOI = 2:6;
BinSize = .5;
PST = [0 .5];
trainTrials = 3:2:19;
testTrials = 2:2:19;
Conc = 1;
numIter = 50;

numPerm{1} = 20;
numPerm{2} = 15;
numPerm{3} = 20;
numPerm{4} = 30;
numPerm{5} = 15;
numPerm{6} = 30;
numPerm{7} = 25;
numPerm{8} = 25;
numPerm{9} = 25;

%%
for k = 1:length(KWIKfiles)
    efd = EFDmaker_Beast(KWIKfiles{k},'bhv');
    
    LRcells = LRcellPicker_chgPt(KWIKfiles{k},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    
    for lset = 1:length(LR_idx)
        Raster = efd.ValveSpikes.RasterAlign(VOI,Conc,LR_idx{lset});
        [trainlabel,traindata{lset}{k},~] = BinRearranger(Raster,PST,BinSize,trainTrials);
        [~,testdata{lset}{k},~] = BinRearranger(Raster,PST,BinSize,testTrials);
       
        for iter = 1:numIter
            templabel = trainlabel;
            randset = randperm(size(traindata{lset}{k},2),numPerm{k});
            tempdata_train = traindata{lset}{k}(:,randset);
            tempdata_test = testdata{lset}{k}(:,randset);
            
            obsindex = 1:length(templabel);
            
            for o = obsindex % iteratively remove population activity vector one at a time
                trl = templabel(obsindex~=o);
                trd = tempdata_train(obsindex~=o,:);
%                 classes = unique(trl,'stable'); % maintain indexing when points are removed           

                model = train(trl,sparse(trd),'-s 4 -q');
                [predict_label{lset}{o}, ~, ~] = predict(1,sparse(tempdata_test(o,:)),model,'-q');
            end
            CM{lset} = confusionmat(trainlabel,cell2mat(predict_label{lset}));
            CM{lset} = bsxfun(@rdivide,CM{lset},sum(CM{lset},2));
            ACC{lset}{k}(iter) = mean(diag(CM{lset}));
        end
    end
end

%%

for lset = 1:2
    for k = 1:length(KWIKfiles)
        meanACC{lset}{k} = mean(ACC{lset}{k});
    end
    meanAcc{lset} = cell2mat(meanACC{lset});
end

%% 

figure; 

subplot(2,2,1); hold on
x = 1:4;
scatter((repmat(x(1),length(offOFF{1}),1)),offOFF{1},'MarkerEdgeColor','none','MarkerFaceColor','r');%,'jitter','on','jitterAmount',.1)
scatter((repmat(x(2),length(offON{1}),1)),offON{1},'MarkerEdgeColor','none','MarkerFaceColor',rgb('Gray'));%,'jitter','on','jitterAmount',.1)

scatter((repmat(x(3),length(offOFF{2}),1)),offOFF{2},'MarkerEdgeColor','none','MarkerFaceColor','r');%,'jitter','on','jitterAmount',.1)
scatter((repmat(x(4),length(offON{2}),1)),offON{2},'MarkerEdgeColor','none','MarkerFaceColor',rgb('Gray'));%,'jitter','on','jitterAmount',.1)

for lset = 1:2
    CI_off{lset} = bootci(1000,{@mean,offOFF{lset}},'Type','per');
    CI_off{lset}(1) = mean(offOFF{lset})-CI_off{lset}(1);
    CI_off{lset}(2) = CI_off{lset}(2)-mean(offOFF{lset});
    
    CI_on{lset} = bootci(1000,{@mean,offON{lset}},'Type','per');
    CI_on{lset}(1) = mean(offON{lset})-CI_on{lset}(1);
    CI_on{lset}(2) = CI_on{lset}(2)-mean(offON{lset});
end

errorbar(1,mean(offOFF{1}),CI_off{1}(1),CI_off{1}(2),'kx')
errorbar(2,mean(offON{1}),CI_on{1}(1),CI_on{1}(2),'kx')
errorbar(3,mean(offOFF{2}),CI_off{2}(1),CI_off{2}(2),'kx')
errorbar(4,mean(offON{2}),CI_on{2}(1),CI_on{2}(2),'kx')

ylim([0 1])
set(gca,'YTick',ylim)
ax = gca; box off; axis square;
ax.XTick = [0 1 2 3 4 5];
ax.XAxis.Limits = [.5 length(x)+.5];

for k = 1:length(KWIKfiles)
    plot([1 2],[offOFF{1}(k) offON{1}(k)])
    plot([3 4],[offOFF{2}(k) offON{2}(k)])
end

    
    