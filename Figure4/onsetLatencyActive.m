function [O,OnLR,OV,OV_nLR,D,DnLR] = onsetLatencyActive(KWIKfiles,Params)

VOI = Params.VOI;
Conc = Params.Conc;
PST = [-Params.BaseTime Params.ResponseTime];
KernelSize = Params.KS;
TOI = Params.TOI;
Cycle = Params.Cycle;

for R = 1:length(KWIKfiles)
    clear efd
    clear KDF
    efd = EFDmaker_Beast(KWIKfiles{R});
    
    Scores = SCOmaker_Beast_PreInh(KWIKfiles{R},{TOI});
    
    LRcells = LRcellPicker_chgPt(KWIKfiles{R},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    if numel(fieldnames(LRcells)) > 2
        LR_idx{2} = sort([LRcells.nonLR,LRcells.secLR]);
    else
        LR_idx{2} = LRcells.nonLR;
    end
        
    [KDF,~,KDFt,~] = KDFmaker_Beast(efd.ValveSpikes.RasterAlign(VOI,Conc,LR_idx{1}),PST,KernelSize);
    [KDF_nLR,~,~,~] = KDFmaker_Beast(efd.ValveSpikes.RasterAlign(VOI,Conc,LR_idx{2}),PST,KernelSize);
    
    aurocLR = Scores.auROC(:,Conc,LR_idx{1},Cycle);
    for valve = 1:(size(KDF,1))
        for unit = 1:size(KDF,3)
            isActiveLR(valve,unit) = aurocLR(valve,:,unit) > 0.5;
        end
    end
    aurocNLR = Scores.auROC(:,Conc,LR_idx{2},Cycle);
    for valve = 1:(size(KDF_nLR,1))
        for unit = 1:size(KDF_nLR,3)
            isActiveNLR(valve,unit) = aurocNLR(valve,:,unit) > 0.5;
        end
    end
    
    realPST = KDFt >= PST(1) & KDFt <= PST(2);
    KDFt = KDFt(realPST);
    
    basePST = KDFt >= KDFt(1) & KDFt <= KDFt(find(KDFt > 0,1));
    stimPST = KDFt > KDFt(find(KDFt > 0,1));
    stimTime = KDFt(stimPST);
    
    for valve = 1:(size(KDF,1))
        for unit = 1:size(KDF,3)
            if  isActiveLR(valve,unit) == 1
                KDF{valve,1,unit} = KDF{valve,1,unit}(realPST);
                meanBase = nanmean(KDF{valve,1,unit}(basePST));
                stdBase = nanstd(KDF{valve,1,unit}(basePST));
                KDFstim = KDF{valve,1,unit}(stimPST);                
                onT_idx = find(KDFstim > 2*stdBase+meanBase,1);
                
                if ~isempty(onT_idx)
                    onT_LR{R}{valve,unit} = stimTime(onT_idx);
                    onV_LR{R}{valve,unit} = KDFstim(onT_idx);
                    offT_idx = find(KDFstim(onT_idx:end) < 1*stdBase+meanBase,1);
                    offT_LR{R}{valve,unit} = stimTime(offT_idx+onT_idx-1);
                    offV_LR{R}{valve,unit} = KDFstim(offT_idx+onT_idx-1);
                    dur_LR{R}{valve,unit} = offT_LR{R}{valve,unit} - onT_LR{R}{valve,unit};
                else
                    onT_LR{R}{valve,unit} = NaN;
                    onV_LR{R}{valve,unit} = NaN;
                    offT_LR{R}{valve,unit} = NaN;
                    offV_LR{R}{valve,unit} = NaN;
                    dur_LR{R}{valve,unit} = NaN;
                end
                
                else
                    onT_LR{R}{valve,unit} = NaN;
                    onV_LR{R}{valve,unit} = NaN;
                    offT_LR{R}{valve,unit} = NaN;
                    dur_LR{R}{valve,unit} = NaN;
                    
            end
        end
    end

    
    for valve = 1:(size(KDF_nLR,1))
        for unit = 1:size(KDF_nLR,3)
            if  isActiveNLR(valve,unit) == 1
                KDF_nLR{valve,1,unit} = KDF_nLR{valve,1,unit}(realPST);
                meanBase = nanmean(KDF_nLR{valve,1,unit}(basePST));
                stdBase = nanstd(KDF_nLR{valve,1,unit}(basePST));
                KDFstim_nLR = KDF_nLR{valve,1,unit}(stimPST);               
                onT_idx_nLR = find(KDFstim_nLR > 2*stdBase+meanBase,1);
                
                if ~isempty(onT_idx_nLR)
                    onT_nLR{R}{valve,unit} = stimTime(onT_idx_nLR);
                    onV_nLR{R}{valve,unit} = KDFstim_nLR(onT_idx_nLR);
                    offT_idx_nLR = find(KDFstim_nLR(onT_idx_nLR:end) < 1*stdBase+meanBase,1);
                    offT_nLR{R}{valve,unit} = stimTime(offT_idx_nLR + onT_idx_nLR - 1);
                    dur_nLR{R}{valve,unit} = offT_nLR{R}{valve,unit} - onT_nLR{R}{valve,unit};
                else
                    onT_nLR{R}{valve,unit} = NaN;
                    onV_nLR{R}{valve,unit} = NaN;
                    offT_nLR{R}{valve,unit} = NaN;
                    dur_nLR{R}{valve,unit} = NaN;
                end
                
            else
                onT_nLR{R}{valve,unit} = NaN;
                onV_nLR{R}{valve,unit} = NaN;
                offT_nLR{R}{valve,unit} = NaN;
                dur_nLR{R}{valve,unit} = NaN;
                               
            end
        end
    end
end    

%%

O = cat(2,onT_LR{:}); O = cell2mat(O(:)');
OnLR = cat(2,onT_nLR{:}); OnLR = cell2mat(OnLR(:)');
Nopks = O < (stimTime(1)+.001) | O > (PST(2)-.001);
Nopks_nLR = OnLR < (stimTime(1)+.001) | OnLR > (PST(2)-.001);
O(Nopks) = nan;
OnLR(Nopks_nLR) = nan;

OV = cat(2,onV_LR{:}); OV = cell2mat(OV(:)');
OV(Nopks) = nan;
OV_nLR = cat(2,onV_nLR{:}); OV_nLR = cell2mat(OV_nLR(:)');
OV_nLR(Nopks_nLR) = nan;

D = cat(2,dur_LR{:}); D = cell2mat(D(:)');
D(Nopks) = nan;
DnLR = cat(2,dur_nLR{:}); DnLR = cell2mat(DnLR(:)');
DnLR(Nopks_nLR) = nan;
