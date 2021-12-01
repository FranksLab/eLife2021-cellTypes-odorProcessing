clearvars
% close all
clc

%% Pick out recording files and put each in one cell

ROI = {'SL'};
Catalog = 'B:\Expt_Sets\eLife2021_DryadData\ExperimentCatalog_Ntng_Inh_Opto.txt';
T = readtable(Catalog, 'Delimiter', ' ');
ROIfiles = T.kwikfile(logical(T.include) & strcmp(T.ROI,ROI));

%% Params

VOI = 2:6;
Cycle = 2;
Conc = 1;
TrialSet{1} = 3:2:19;
TrialSet{2} = 2:2:19;

%% 

for R = 1:length(ROIfiles)
    for tset = 1:length(TrialSet)
        clear Scores
        clear LRcells
        Scores = SCOmaker_Beast_PreInh(ROIfiles{R},TrialSet(tset));

        LRcells = LRcellPicker_chgPt(ROIfiles{R},[-.1 .1]);
        LR_idx{1} = LRcells.primLR;
        LR_idx{2} = LRcells.nonLR;

        sig_act{R,tset} = Scores.auROC(VOI,Conc,LR_idx{1},Cycle) > 0.5 & Scores.AURp(VOI,Conc,LR_idx{1},Cycle) < 0.05;
        sig_sup{R,tset} = Scores.auROC(VOI,Conc,LR_idx{1},Cycle) < 0.5 & Scores.AURp(VOI,Conc,LR_idx{1},Cycle) < 0.05;
        ratechange_LR{R,tset} = Scores.RawRate(VOI,Conc,LR_idx{1},Cycle);

        sig_act_nLR{R,tset} = Scores.auROC(VOI,Conc,LR_idx{2},Cycle) > 0.5 & Scores.AURp(VOI,Conc,LR_idx{2},Cycle) < 0.05;
        sig_sup_nLR{R,tset} = Scores.auROC(VOI,Conc,LR_idx{2},Cycle) < 0.5 & Scores.AURp(VOI,Conc,LR_idx{2},Cycle) < 0.05;
        ratechange_nLR{R,tset} = Scores.RawRate(VOI,Conc,LR_idx{2},Cycle);
    end
end

%% Reshape and stack

for tset = 1:length(TrialSet)
    ratechg_LR{tset} = reshape(cat(3,ratechange_LR{:,tset}),[],1);
    ratechg_nLR{tset} = reshape(cat(3,ratechange_nLR{:,tset}),[],1);
    act_LR{tset} = reshape(cat(3,sig_act{:,tset}),[],1);
    sup_LR{tset} = reshape(cat(3,sig_sup{:,tset}),[],1);
    act_nLR{tset} = reshape(cat(3,sig_act_nLR{:,tset}),[],1);
    sup_nLR{tset} = reshape(cat(3,sig_sup_nLR{:,tset}),[],1);
end

ratechg_LR = cell2mat(ratechg_LR);
ratechg_nLR = cell2mat(ratechg_nLR);
act_LR = cell2mat(act_LR);
sup_LR = cell2mat(sup_LR);
act_nLR = cell2mat(act_nLR);
sup_nLR = cell2mat(sup_nLR);

%% Plotting LR

for m = 1:length(ratechg_LR)
    if any(act_LR(m,1))
        ups{m} = ratechg_LR(m,:);
    elseif any(sup_LR(m,1))
        downs{m} = ratechg_LR(m,:);
    else
        stat{m} = ratechg_LR(m,:);
    end
end

%% Plotting NLR

for m = 1:length(ratechg_nLR)
    if any(act_nLR(m,1))
        ups_nLR{m} = ratechg_nLR(m,:);
    elseif any(sup_nLR(m,1))
        downs_nLR{m} = ratechg_nLR(m,:);
    else
        stat_nLR{m} = ratechg_nLR(m,:);
    end
end

%%

for tset = 1:2
    CI_nLR_ups{tset} = bootci(1000,{@mean,ups_nLR(:,tset)},'Type','per');
    CI_nLR_ups{tset}(1) = mean(ups_nLR(:,tset))-CI_nLR_ups{tset}(1);
    CI_nLR_ups{tset}(2) = CI_nLR_ups{tset}(2)- mean(ups_nLR(:,tset));
    
    CI_nLR_downs{tset} = bootci(1000,{@mean,downs_nLR(:,tset)},'Type','per');
    CI_nLR_downs{tset}(1) = mean(downs_nLR(:,tset))-CI_nLR_downs{tset}(1);
    CI_nLR_downs{tset}(2) = CI_nLR_downs{tset}(2)-mean(downs_nLR(:,tset));
    
    CI_nLR_stat{tset} = bootci(1000,{@mean,stat_nLR(:,tset)},'Type','per');
    CI_nLR_stat{tset}(1) = mean(stat_nLR(:,tset))-CI_nLR_stat{tset}(1);
    CI_nLR_stat{tset}(2) = CI_nLR_stat{tset}(2)- mean(stat_nLR(:,tset));
    
    CI_LR_ups{tset} = bootci(1000,{@mean,ups(:,tset)},'Type','per');
    CI_LR_ups{tset}(1) = mean(ups(:,tset))-CI_LR_ups{tset}(1);
    CI_LR_ups{tset}(2) = CI_LR_ups{tset}(2)- mean(ups(:,tset));
    
    CI_LR_downs{tset} = bootci(1000,{@mean,downs(:,tset)},'Type','per');
    CI_LR_downs{tset}(1) = mean(downs(:,tset))-CI_LR_downs{tset}(1);
    CI_LR_downs{tset}(2) = CI_LR_downs{tset}(2)-mean(downs(:,tset));
    
    CI_LR_stat{tset} = bootci(1000,{@mean,stat(:,tset)},'Type','per');
    CI_LR_stat{tset}(1) = mean(stat(:,tset))-CI_LR_stat{tset}(1);
    CI_LR_stat{tset}(2) = CI_LR_stat{tset}(2)- mean(stat(:,tset));
end

figure; subplot(2,2,1); hold on;
for tset = 1:2
errorbar(tset,mean(ups(:,tset)),CI_LR_ups{tset}(1),CI_LR_ups{tset}(2),'rx')
errorbar(tset,mean(downs(:,tset)),CI_LR_downs{tset}(1),CI_LR_downs{tset}(2),'bx')
errorbar(tset,mean(stat(:,tset)),CI_LR_stat{tset}(1),CI_LR_stat{tset}(2),'kx')
end
ax = gca; ax.YAxis.Limits = [0 12];

subplot(2,2,2); hold on;
for tset = 1:2
errorbar(tset,mean(ups_nLR(:,tset)),CI_nLR_ups{tset}(1),CI_nLR_ups{tset}(2),'rx')
errorbar(tset,mean(downs_nLR(:,tset)),CI_nLR_downs{tset}(1),CI_nLR_downs{tset}(2),'bx')
errorbar(tset,mean(stat_nLR(:,tset)),CI_nLR_stat{tset}(1),CI_nLR_stat{tset}(2),'kx')
end
ax = gca; ax.YAxis.Limits = [0 12];



