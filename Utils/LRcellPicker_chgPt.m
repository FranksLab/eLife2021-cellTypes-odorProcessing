function [LR] = LRcellPicker_chgPt(KWIKfile,PST)
% function [LR,auROC,changept,p] = LRcellPicker_chgPt(KWIKfile,PST)

efd = EFDmaker_Beast(KWIKfile,'bhv');
%% Load file

[a,b] = fileparts(KWIKfile);
LRfile = [a,filesep,b,'.lr'];
if exist(LRfile,'file')
    load(LRfile,'-mat')
else
    
%% Identify LR and nonLR cells

for unit = 1:size(efd.LaserSpikes.SpikesBeforeLaser,2)
    Control = efd.LaserSpikes.SpikesBeforeLaser{unit};
    Stimulus = efd.LaserSpikes.SpikesDuringLaser{unit};
    [auROC(unit), p(unit)] = RankSumROC(Control, Stimulus); % ranksum select laser responsive cells
    LR_p(unit) =  auROC(unit) < .5 & p(unit) < .00000001;
end
 
%% Change-point analysis

numTrials = length(efd.LaserSpikes.RasterAlign{1});

for unit = 1:size(efd.LaserSpikes.RasterAlign,2)
    [PSTH, ~, PSTHt] = PSTHmaker_Beast(efd.LaserSpikes.RasterAlign(unit), PST, 0.002, 141:numTrials);
    ipt = findchangepts(PSTH{:});
    changept(unit) = PSTHt(ipt);
    LR_lat(unit) = changept(unit) > -.05 & changept(unit) <= .01;
end

%%

LR_idx = find(LR_p == 1 & LR_lat == 1);
nonLR_idx = find(LR_p == 0 | LR_lat == 0);


LR.primLR = LR_idx;
LR.nonLR = nonLR_idx;

save(LRfile,'LR')

end