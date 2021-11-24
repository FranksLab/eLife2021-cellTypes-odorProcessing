function [La,LanLR,D,DnLR,Pk,PknLR,mx,mx_nLR,mxt,mxt_nLR] = LatDurPeak(KWIKfiles,Params,Trials)

VOI = Params.VOI;
Conc = Params.Conc;
PST = Params.PST;
KernelSize = Params.KS;
TOI = Trials;

for R = 1:length(KWIKfiles)
    clear efd
    clear KDF
    efd = EFDmaker_Beast(KWIKfiles{R});
%     Scores = SCOmaker_Beast(KWIKfiles{R},{TOI});
    LRcells = LRcellPicker_chgPt(KWIKfiles{R},[-.1 .1]);
    LR_idx{1} = LRcells.primLR;
    if numel(fieldnames(LRcells)) > 2
        LR_idx{2} = sort([LRcells.nonLR,LRcells.secLR]);
    else
        LR_idx{2} = LRcells.nonLR;
    end
    [KDF,~,KDFt,~] = KDFmaker_Beast(efd.ValveSpikes.RasterAlign(VOI,Conc,LR_idx{1}),PST,KernelSize,TOI);
    [KDF_nLR,~,~,~] = KDFmaker_Beast(efd.ValveSpikes.RasterAlign(VOI,Conc,LR_idx{2}),PST,KernelSize,TOI);
    
    realPST = KDFt>=PST(1) & KDFt<=PST(2);
    KDFt = KDFt(realPST);

    
    for valve = 1:size(KDF,1)
        for unit = 1:size(KDF,3)
            if ~isempty(KDF{valve,1,unit})
                KDF{valve,1,unit} = KDF{valve,1,unit}(realPST);
                [mx{R}(valve,unit),mxt_temp] = max(KDF{valve,1,unit});
                mxt{R}(valve,unit) = KDFt(mxt_temp);
                mxd{R}(valve,unit) = FWHM(KDFt,KDF{valve,1,unit},.5);
            else
                KDF{valve,unit}  = nan(1,length(realPST));
                mx{R}(valve,unit) = nan;
                mxt{R}(valve,unit) = nan;
                mxd{R}(valve,unit) = nan;
            end
        end
    end
    
    for valve = 1:size(KDF_nLR,1)
        for unit = 1:size(KDF_nLR,3)
            if ~isempty(KDF_nLR{valve,1,unit})
                KDF_nLR{valve,1,unit} = KDF_nLR{valve,1,unit}(realPST);
                [mx_nLR{R}(valve,unit),mxt_temp] = max(KDF_nLR{valve,1,unit});
                mxt_nLR{R}(valve,unit) = KDFt(mxt_temp);
                mxd_nLR{R}(valve,unit) = FWHM(KDFt,KDF_nLR{valve,1,unit},.5);
            else
                KDF_nLR{valve,unit}  = nan(1,length(realPST));
                mx_nLR{R}(valve,unit) = nan;
                mxt_nLR{R}(valve,unit) = nan;
                mxd_nLR{R}(valve,unit) = nan;
            end
        end
    end
end

%%

L = cat(2,mxt{:}); La = L(:);
LnLR = cat(2,mxt_nLR{:}); LanLR = LnLR(:);
Nopks = La < (PST(1)+.001) | La > (PST(2)-.001);
Nopks_nLR = LanLR < (PST(1)+.001) | LanLR > (PST(2)-.001);
La(Nopks) = nan;
LanLR(Nopks_nLR) = nan;

P = cat(2,mx{:}); Pk = P(:);
PnLR = cat(2,mx_nLR{:}); PknLR = PnLR(:);
Pk(Nopks) = nan;
PknLR(Nopks_nLR) = nan;

D = cat(2,mxd{:}); D = D(:);
DnLR = cat(2,mxd_nLR{:}); DnLR = DnLR(:);
D(Nopks) = nan;
DnLR(Nopks_nLR) = nan;
