clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng.txt'; % set to data and catalog directory
T = readtable(Catalog, 'Delimiter', ' ');
ROIfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

Params.PST = [-.5 .5];
Params.KS = 0.002; % kernel size
Params.TOI = 1:1000; % trials
% TOI for Figure1-S5c: 1:1000
% TOI for Figure1-S5e (top): 141:1140
% TOI for Figure1-S5c (bottom): all
% TOI for Figure1-S5c: 1:1000

%% Response index

for R = 1:length(ROIfiles)
    clear efd
    efd = EFDmaker_Beast(ROIfiles{R},'bhv');
    
    LRcells = LRcellPicker_chgPt(ROIfiles{R},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    LR_idx{2} = LRcells.nonLR;
    
    [KDF{R}, ~, ~] = PSTHmaker_Beast(efd.LaserSpikes.RasterAlign(LR_idx{1}), Params.PST, Params.KS, Params.TOI);    
    [KDF_nLR{R}, ~, KDFt] = PSTHmaker_Beast(efd.LaserSpikes.RasterAlign(LR_idx{2}), Params.PST, Params.KS, Params.TOI);    
end

% realPST = KDFt>=Params.PST(1) & KDFt<=Params.PST(2);
% KDFt = KDFt(realPST);    
    
LR = cat(2,KDF{:});
nLR = cat(2,KDF_nLR{:});

lr = cell2mat(LR');
nlr = cell2mat(nLR');

%%

figure;
imagesc(KDFt,[],lr);
HT = hot(32);
HT = HT(1:end-4,:);
colormap(HT)
caxis([0 15])
        
figure;
imagesc(KDFt,[],nlr);
HT = hot(32);
HT = HT(1:end-4,:);
colormap(HT)
caxis([0 10])

