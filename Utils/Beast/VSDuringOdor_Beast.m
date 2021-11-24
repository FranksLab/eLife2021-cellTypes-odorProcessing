function [SpikesDuringOdor] = VSDuringOdor_Beast(ValveTimes,SpikeTimes)

SpikesDuringOdor = cell(size(ValveTimes.PREXTimes,1),size(ValveTimes.PREXTimes,2),size(SpikeTimes.tsec,1));

for i = 1:size(ValveTimes.FVSwitchTimesOff,1)
    a(i) = size(ValveTimes.FVSwitchTimesOff{i},2);
end
maxa = max(a);


for Unit = 1:size(SpikeTimes.tsec,1)
    st = SpikeTimes.tsec{Unit};
    for Valve = 1:size(ValveTimes.PREXTimes,1) 
        for Conc = 1:size(ValveTimes.PREXTimes,2)
            if ~isempty(ValveTimes.PREXTimes{Valve})
                Opening = ValveTimes.FVSwitchTimesOn{Valve,Conc}(:);
                Closing = ValveTimes.FVSwitchTimesOff{Valve,Conc}(:);
                fvl = min(length(Opening),length(Closing));
                Opening = Opening(1:fvl);
                Closing = Closing(1:fvl);

                x = bsxfun(@gt,st,Opening');
                x2 = bsxfun(@lt,st,Closing');
                x3 = x+x2-1;

                SpikesDuringOdor{Valve,Conc,Unit} = sum(x3==1);
                SpikesDuringOdor{Valve,Conc,Unit}(maxa+(a(Valve)-maxa+1):maxa) = NaN;
            end
        end
    end
end