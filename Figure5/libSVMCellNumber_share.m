clearvars
% close all
clc

%% Pick out files with 'kwik' in its name and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_SL_inh_workingcopy.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

VOI = 2:6;
BinSize = .5;
PST = [0 .5];
PST_pre = [-.6 -.1];
Trials = 2:2:20;
Conc = 1;
numPerm = 150;
numIter = 100;

%% Rearrange data and run classifier

for k = 1:length(KWIKfiles)
    efd = EFDmaker_Beast(KWIKfiles{k},'bhv');
%     SpikeTimes = SpikeTimes_Beast(FindFilesKK(KWIKfiles{k}));
%     [~,pt,~,~,~,~] = WaveStats_Beast(SpikeTimes.Wave);
    
    LRcells = LRcellPicker_chgPt(KWIKfiles{k},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    for lset = 1:length(LR_idx)
%         pttime{lset} = pt(LR_idx{lset});
        Raster = efd.ValveSpikes.RasterAlign(VOI,Conc,LR_idx{lset});
%         Raster = Raster(:,:,pttime{lset} > .4);
%         RasterMO = efd.ValveSpikes.RasterAlign(1,Conc,LR_idx{lset});
%         RasterMO = RasterMO(:,:,pttime{lset} > .4);
        [trainlabel,data{lset}{k},~] = BinRearranger(Raster,PST,BinSize,Trials);
        [~,dataPre{lset}{k},~] = BinRearranger(Raster,PST_pre,BinSize,Trials);
%         [~,dataMO{lset}{k},~] = BinRearranger(RasterMO,PST,BinSize,Trials);
        data_zsc{lset}{k} = (data{lset}{k}-mean(dataPre{lset}{k},1))./std(dataPre{lset}{k},1);
%         data_center{lset}{k} = (data{lset}{k}-mean(data{lset}{k},2));
        data_zsc{lset}{k}(isinf(data_zsc{lset}{k})) = nan;
%         data_zsc{lset}{k}(isnan(data_zsc{lset}{k})) = 0;  
    end
end

%%  

for lset = 1:length(LR_idx)
    allData{lset} = cell2mat(data{lset});
    allData_zsc{lset} = cell2mat(data_zsc{lset});
    data_shuffle{lset} = allData_zsc{lset}(randperm(size(allData_zsc{lset},1)),:);
    
%     for ii = 1:size(allData{lset},2)
% %         for o = 1:length(VOI)
% %             odorTrials = length(Trials)*(o-1)+1:length(Trials)*o;
%             unit = allData{lset}(:,ii);
%             allData_norm{lset}(:,ii) = (unit-min(unit))/(max(unit)-min(unit));
%             allData_norm{lset}(isnan(allData_norm{lset})) = 0;
% %         end
%     end
end

%%

for lset = 1:length(LR_idx)
    for b = 50%1:numPerm        
        for iter = 1:numIter
%             clearvars randset_bin; clearvars rs_bin;
%             templabel = trainlabel;
%             randset = randperm((size(allData{lset},2)/16),b);
%             for rr = 1:length(randset)
%                 randset_bin{rr} = randset(rr)*16-16+1:randset(rr)*16;
%             end
%             rs_bin = cell2mat(randset_bin);
%             tempdata = allData{lset}(:,rs_bin);
            templabel = trainlabel;
            randset = randperm(size(allData{lset},2),b);
            tempdata = allData{lset}(:,randset);
%             randset = randperm(size(data_shuffle{lset},2),b);
%             tempdata = data_shuffle{lset}(:,randset);
            
            obsindex = 1:length(templabel);
            
            for o = obsindex % iteratively remove population activity vector one at a time
                trl = templabel(obsindex~=o);
                trd = tempdata(obsindex~=o,:);
                classes = unique(trl,'stable'); % maintain indexing when points are removed           
                
%                 model = svmtrain(trl,sparse(trd),'-s 0 -t 0 -q');
%                 [predict_label{lset}{o}, ~, ~] = svmpredict(1,sparse(tempdata(o,:)),model,'-q');      

                model = train(trl,sparse(trd),'-s 4 -q');
                [predict_label{lset}{o}, ~, ~] = predict(1,sparse(tempdata(o,:)),model,'-q');
            end
            CM{lset} = confusionmat(trainlabel,cell2mat(predict_label{lset}));
            CMat{lset}{iter} = bsxfun(@rdivide,CM{lset},sum(CM{lset},2));
            ACC{lset}(b,iter) = mean(diag(CMat{lset}{iter}));
        end
    end
end

%% plotting

figure; subplot(2,2,1); hold on
colors = {[0.7188 0.5234 0.0430; 0 0 0], [.3 .3 .3; 0 0 0]};

for lset = 1:length(LR_idx)
    x=1:numPerm;
    boundedline(x,mean(ACC{lset}(1:end,:),2),sem(ACC{lset}(1:end,:)'),'cmap',colors{lset})
    hold on
    box off; ylim([0 1]); xlim([0 400]);
    axis square;
end

% significance testing
% for b = 10:numPerm
%     [pVal(b),~] = ranksum(ACC{1}(b,:),ACC{2}(b,:));
%     if pVal(b) < .05
%         plot(b,.1,'r*')
%     end
% end
% 
% subplot(2,2,2); hold on
% for lset = 1:2
%     errorbar(mean(ACC{lset}(end,:)),std(ACC{lset}(end,:)),'Marker','x')
%     box off; ylim([0 1.2]); axis square;
% end

%%

figure;
for lset = 1:length(LR_idx)
    cm_iter{lset} = cat(3,CMat{lset}{:});
    meanCM{lset} = mean(cm_iter{lset},3);
    
    subplot(2,2,lset);
    imagesc(meanCM{lset}); colormap(jet); colorbar
    axis square; box off; axis on; caxis([0 1]);
end

    
    