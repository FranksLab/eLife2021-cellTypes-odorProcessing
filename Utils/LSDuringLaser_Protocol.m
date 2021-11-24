function [SpikesDuringLaser, SpikesBeforeLaser, SpikesDuringLaserLate] = LSDuringLaser_Protocol(LaserTimes,SpikeTimes)

% SpikesDuringLaser = cell(size(LaserTimes.LaserOff,2),size(SpikeTimes.tsec,1));

for Unit = 1:size(SpikeTimes.tsec,1)
    st = SpikeTimes.tsec{Unit};
    
    Opening_all = LaserTimes.LaserOn{1}(:);
%     Closing = LaserTimes.LastPulseOff{1}(:);
    Closing_all = LaserTimes.LaserOff{1}(:);
    fvl = min(length(Opening_all),length(Closing_all));
    Opening_all = Opening_all(1:fvl);
    Closing_all = Closing_all(1:fvl);
    Diff = Closing_all-Opening_all;
    Protocol_idx = find(Diff<1,1); % Looking for the start of the tagging protocol (pulses<1s) in experiments that use laser not just for tagging. 
    if ~isempty(Protocol_idx)
        Opening = Opening_all(Protocol_idx:end);
        Closing = Closing_all(Protocol_idx:end);
    else
        Opening = Opening_all;
        Closing = Closing_all;
    end
    PreOpening = Opening - mean(Closing-Opening);
    LateOpening = Opening + mean(Closing-Opening)/2;
    
    %% Count spikes during laser
    x = bsxfun(@gt,st,Opening'); % logical True spikes greater than On time 
    x2 = bsxfun(@lt,st,Closing'); % logical True spikes less than Off time 
    x3 = x+x2-1; % True + True -1 = 1; True + False -1 = 0;
    SpikesDuringLaser{Unit} = sum(x3==1);
    
    %% Count spikes during laser second half (because some cells only transiently activate)
    x = bsxfun(@gt,st,LateOpening'); % logical True spikes greater than On time 
    x2 = bsxfun(@lt,st,Closing'); % logical True spikes less than Off time 
    x3 = x+x2-1; % True + True -1 = 1; True + False -1 = 0;
    SpikesDuringLaserLate{Unit} = sum(x3==1);
    
    %% Count spikes before laser
    x = bsxfun(@gt,st,PreOpening'); % logical True spikes greater than On time 
    x2 = bsxfun(@lt,st,Opening'); % logical True spikes less than Off time 
    x3 = x+x2-1; % True + True -1 = 1; True + False -1 = 0;
    SpikesBeforeLaser{Unit} = sum(x3==1);

end