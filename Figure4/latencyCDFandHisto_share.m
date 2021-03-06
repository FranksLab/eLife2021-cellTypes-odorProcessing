clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\catalog\ExperimentCatalog_Ntng.txt';
T = readtable(Catalog, 'Delimiter', ' ');
KWIKfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

Params.PST = [0 .3];
Params.KS = 0.01;
Params.VOI = 2:11;
Params.Conc = 3;
Params.TOI = 2:11;
Params.Cycle = 3;

Params.BaseTime = .3;
Params.ResponseTime = .3;

%% Determine latency, duration and peak

% [Lat,Lat_nLR,~,~,~,~] = peakLatencyActive(KWIKfiles,Params);
[Lat,Lat_nLR,~,~,~,~] = onsetLatencyActive(KWIKfiles,Params);

Lat = Lat(~isnan(Lat));
Lat_nLR = Lat_nLR(~isnan(Lat_nLR));

%% Plotting CDF

edges = 0:.002:.3;

[NLat,edges] = histcounts(Lat,edges);
[NLat_nLR,edges] = histcounts(Lat_nLR,edges);
cdfLat = cumsum(NLat)/sum(NLat);
cdfLat_nLR = cumsum(NLat_nLR)/sum(NLat_nLR);

colors = {rgb('ForestGreen'),rgb('DarkGoldenrod')};
figure; hold on

subplot(2,2,1); hold on
plot(edges(1:end-1),cdfLat,'Color',colors{1})
plot(edges(1:end-1),cdfLat_nLR,'Color',rgb('Gray'))
xlim([0 .3])
ylim([0 1])
set(gca,'XTick',xlim)
set(gca,'YTick',ylim)
xlabel('latency (s)')
axis square; box off;

%% Plotting histograms

subplot(2,2,2); hold on
edges = 0:.002:.1;
histogram(Lat,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k','LineWidth',1);
histogram(Lat_nLR,edges,'Normalization','probability','DisplayStyle','stairs','EdgeColor',colors{1},'LineWidth',1);
xlim([0 .1])
set(gca,'XTick',xlim)
set(gca,'YTick',ylim)
xlabel('latency (s)')
ylabel('fraction responses')
axis square; box off;
