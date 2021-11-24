function [ValveSpikes,LaserSpikes,LVSpikes] = VSmaker_Beast(ValveTimes,LaserTimes,LVTimes,SpikeTimes,PREX)

%% Odor processing
if isfield(ValveTimes,'PREXTimes')
    %% Aligned Raster
    [ValveSpikes.RasterAlign] = VSRasterAlign_Beast(ValveTimes,SpikeTimes);
    %% Spikes in Multi Cycles
    [ValveSpikes.MultiCycleSpikeCount,ValveSpikes.MultiCycleSpikeRate,ValveSpikes.MultiCycleBreathPeriod] = VSMultiCycleCount_Beast(ValveTimes,SpikeTimes,PREX,{0:2});
    %% Spikes During Odor
    ValveSpikes.SpikesDuringOdor = VSDuringOdor_Beast(ValveTimes,SpikeTimes);
else
    ValveSpikes = [];
end

%% Laser processing
if isfield(LaserTimes,'LaserOn')
    %% Laser Aligned Raster
    [LaserSpikes.RasterAlign] = VSRasterAlignLaser(LaserTimes,SpikeTimes);
    %% Spikes During Laser
    [LaserSpikes.SpikesDuringLaser, LaserSpikes.SpikesBeforeLaser, LaserSpikes.SpikesDuringLaserLate] = LSDuringLaser_Protocol(LaserTimes,SpikeTimes);
else
    LaserSpikes=[];
end

%% Odor and Laser processing together
% if ~isstr(LVTimes)
%     for LS = 1:2
%         %% Laser Aligned Raster
%         [LVSpikes{LS}.RasterAlign] = VSRasterAlign(LVTimes{LS},SpikeTimes);
%     end
% else
    LVSpikes = [];
% end

end
