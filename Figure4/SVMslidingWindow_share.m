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
BinSize = .03;
PST = [0 .5];
pcaTrials = 2:11;
Conc = 3;
numPerm = 400;
numIter = 100;

%% Window-making

WindowType = 'Sliding';

switch WindowType
    case 'Sliding'
        stepSize = BinSize/2;
        windowend = PST(1):stepSize:PST(2); 
        win_t = windowend+BinSize/2;
        
    case 'Expanding'
        windowend = (PST(1)+BinSize):BinSize:PST(2); 
        win_t = windowend-BinSize/2;
end

%% Rearrange data and run classifier

for w = 1:length(windowend)
    for k = 1:length(KWIKfiles)
        efd = EFDmaker_Beast(KWIKfiles{k},'bhv');
        
        LRcells = LRcellPicker_chgPt(KWIKfiles{k},[-.1 .1]);
        LR_idx{1} = LRcells.primLR;
        LR_idx{2} = LRcells.nonLR;
        for lset = 1:length(LR_idx)
            Raster = efd.ValveSpikes.RasterAlign(VOI,Conc,LR_idx{lset});

            switch WindowType
                case 'Sliding'
                    [trainlabel,pcadata{lset}{k},~] = BinRearranger(Raster,[0 BinSize]+windowend(w),BinSize,pcaTrials);
                case 'Expanding'
                    [trainlabel,pcadata{lset}{k},~] = BinRearranger(Raster,[0 windowend(w)],BinSize,pcaTrials);
            end
        end
    end
    
    for lset = 1:length(LR_idx)
        allData{lset} = cell2mat(pcadata{lset}); % combine cells from all experiments
        
        for iter = 1:numIter
        templabel = trainlabel;
        randset = randperm(size(allData{lset},2),numPerm);
        tempdata = allData{lset}(:,randset);
        
        obsindex = 1:length(templabel);
        
        for o = obsindex % iteratively remove population activity vector one at a time
            trl = templabel(obsindex~=o);
            trd = tempdata(obsindex~=o,:);
            clslist = unique(trl,'stable'); % maintain indexing when points are removed
            
            model = train(trl,sparse(trd),'-s 4 -q');
            [predict_label{lset}{o}, ~, ~] = predict(1,sparse(tempdata(o,:)),model,'-q');             
        end
        
        CM{lset} = confusionmat(trainlabel,cell2mat(predict_label{lset}));
        CM{lset} = bsxfun(@rdivide,CM{lset},sum(CM{lset},2));
        ACC{lset}(w,iter) = mean(diag(CM{lset}));
        end
    end
end

%% plotting

figure; subplot(2,2,1); hold on
colors = {[0.7188 0.5234 0.0430; 0 0 0], [.3 .3 .3; 0 0 0]};
figure;

for lset = 1:length(LR_idx)
    for w = 1:length(windowend)
        CI{lset}{w} = bootci(1000,{@mean,ACC{lset}(w,:)},'Type','per');
        CI{lset}{w}(1) = mean(ACC{lset}(w,:))-CI{lset}{w}(1);
        CI{lset}{w}(2) = CI{lset}{w}(2)-mean(ACC{lset}(w,:));
    end
end
    
for lset = 1:length(LR_idx)
    boundedline(win_t,mean(ACC{lset},2)',cell2mat(CI{lset})','cmap',colors{lset})
    hold on
    box off; ylim([0 1]); xlim([0 .1]);
    axis square;
end