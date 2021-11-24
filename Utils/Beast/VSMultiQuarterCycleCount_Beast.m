function [MultiCycleSpikeCount,MultiCycleSpikeRate,MultiCycleBreathPeriod] = VSMultiQuarterCycleCount_Beast(ValveTimes,SpikeTimes,PREX,CyclestoCheck)
if ~iscell(CyclestoCheck) == 1
    CycleList = 1:CyclestoCheck;
else
    CycleList = cell2mat(CyclestoCheck);
    CyclestoCheck = length(CycleList);
end

MultiCycleSpikeCount = cell(size(ValveTimes.PREXIndex,1),size(ValveTimes.PREXIndex,2),size(SpikeTimes.tsec,1),CyclestoCheck);
MultiCycleSpikeRate = cell(size(ValveTimes.PREXIndex,1),size(ValveTimes.PREXIndex,2),size(SpikeTimes.tsec,1),CyclestoCheck);
MultiCycleBreathPeriod = cell(CyclestoCheck,1);

for i = 1:size(ValveTimes.PREXIndex,1)
    a(i) = size(ValveTimes.PREXIndex{i},2);
end
maxa = max(a);

for Unit = 1:size(SpikeTimes.tsec,1)
    st = SpikeTimes.tsec{Unit};
    
    for Valve = 1:size(ValveTimes.PREXIndex,1)
        if ~isempty(ValveTimes.PREXTimes{Valve})
            for Conc = 1:size(ValveTimes.PREXIndex,2)                
                for Cycle = 1:length(CycleList)
                    Adder = CycleList(Cycle)-1;
                    ValveTimes.PREXIndex{Valve,Conc}((ValveTimes.PREXIndex{Valve,Conc}(:)+1+Adder)>length(PREX)) = [];
                    Beginning = PREX(ValveTimes.PREXIndex{Valve,Conc}(:)+Adder);
                    EndofCycle = PREX(ValveTimes.PREXIndex{Valve,Conc}(:)+1+Adder);
                    MultiCycleBreathPeriod{Valve,Conc,Cycle} = EndofCycle-Beginning;
                    QuarterPeriod = (EndofCycle-Beginning)/4;
                    try
                        x = bsxfun(@gt,st,Beginning);
                        x2 = bsxfun(@lt,st,(Beginning+QuarterPeriod));
                        x3 = x+x2-1;
                        
                        MultiCycleSpikeCount{Valve,Conc,Unit,Cycle} = sum(x3==1);
                        MultiCycleSpikeCount{Valve,Conc,Unit,Cycle}(maxa+(a(Valve)-maxa+1):maxa) = NaN;
                        MultiCycleSpikeRate{Valve,Conc,Unit,Cycle} = MultiCycleSpikeCount{Valve,Conc,Unit,Cycle}./(EndofCycle-Beginning);
                    catch
                        for Trial = 1:size(ValveTimes.PREXIndex{Valve},2)
                            % Loop through the trials and if the cycle is too far
                            % ... or just use PREX and CyclestoCheck to figure it
                            % out.. Then make those NaNs.
                            MultiCycleSpikeCount{Valve,Conc,Unit,Cycle}(Trial) = sum(st>Beginning(Trial) & st<EndofCycle(Trial));
                            MultiCycleSpikeRate{Valve,Conc,Unit,Cycle}(Trial) = sum(st>Beginning(Trial) & st<EndofCycle(Trial))/(EndofCycle(Trial)-Beginning(Trial));
                            
                        end
                        MultiCycleSpikeCount{Valve,Conc,Unit,Cycle}(maxa+(a(Valve)-maxa+1):maxa) = NaN;
                        MultiCycleSpikeRate{Valve,Conc,Unit,Cycle}(maxa+(a(Valve)-maxa+1):maxa) = NaN;
                        
                    end
                end
            end
        end
    end
    
    
end