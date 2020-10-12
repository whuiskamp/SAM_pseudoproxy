%% ENSO regression maps
% This script, based on the ENSO indices calculated in 'ENSO_stuff.m',
% regresses these indices against the running correlation between SAM and
% SAT/Precip at each grid cell. 
%
% May, 2020

% Setup
clear
load('DataFiles/model_output.mat','lat');
load('DataFiles/filtered_n34.mat');
land = ncread('DataFiles/sftlf_A1.static.nc','sftlf')';
load('DataFiles/model_output.mat','sat_detr','SAM','precip_detr');
lat_2 = 45; % change this from 45 to get some sub-set of the SH
n34_SAT_regr = nan(3,lat_2,144); n34_pre_regr = nan(3,lat_2,144);
n34 = nan(3,499); n34(1,:) = filt_ann_30; n34(2,:) = filt_ann_60;
n34(3,:) = filt_ann_90;


for windowsize = [31 61 91]
    load(['DataFiles/runcorr',num2str(windowsize),'yrwdw.mat']);
    for i = 1:size(sat_runcorr,1)
        sat_runcorr(i,isnan(land)) = nan;
        precip_runcorr(i,isnan(land)) = nan;
    end
    range = ceil(windowsize/2):500-floor(windowsize/2);
    sat_runcorr = sat_runcorr(range,1:lat_2,:); % we only need the SH and eliminate the NaNs from the running correlation at either end.
    pre_runcorr = precip_runcorr(range,1:lat_2,:);
    
    for i = 1:size(sat_runcorr,2)
        for j = 1:size(sat_runcorr,3)
            if ~isnan(sat_runcorr(1,i,j))
               p = polyfit(zscore(n34(floor(windowsize/30),range))',zscore(sat_runcorr(:,i,j)),1); n34_SAT_regr(floor(windowsize/30),i,j) = p(1);
               p = polyfit(zscore(n34(floor(windowsize/30),range))',zscore(pre_runcorr(:,i,j)),1); n34_pre_regr(floor(windowsize/30),i,j) = p(1);
            end
        end
    end
end

SAM_SAT_regr = nan(90,144);
SAM_precip_regr = nan(90,144);
% Calculate regression slopes for SAM-SAT/precip
for i = 1:size(sat_detr,2)
    for j = 1:size(sat_detr,3)
        if ~isnan(land(i,j))
            p = polyfit(zscore(SAM),zscore(sat_detr(:,i,j)),1); SAM_SAT_regr(i,j) = p(1);
            p = polyfit(zscore(SAM),zscore(precip_detr(:,i,j)),1); SAM_precip_regr(i,j) = p(1);
        end
    end
end

% Create a scatter plot of ENSO regression slopes vs SAM regressions slopes
save('DataFiles/ENSO_regression.mat','n34_SAT_regr','n34_pre_regr')
save('DataFiles/SAM_regression.mat','SAM_SAT_regr','SAM_precip_regr')

%% Once the above has been run, simply start here

clear
load('DataFiles/ENSO_corrs.mat','SAT_corr_sig','pre_corr_sig')
load('DataFiles/ENSO_regression.mat','n34_SAT_regr','n34_pre_regr')
load('DataFiles/SAM_regression.mat','SAM_SAT_regr','SAM_precip_regr')
land = ncread('DataFiles/sftlf_A1.static.nc','sftlf')';

figure(1)
% SAT
for w = [1 2 3]
    subplot(2,3,w)
    hold on
    xlabel('SAM-SAT regr. coeff.')
    ylabel('n34-SAM/SAT running corr regr. coeff')
    for i = 1:size(SAT_corr_sig,2)
        for j = 1:size(SAT_corr_sig,3)
            if isnan(SAT_corr_sig(w,i,j)) && ~isnan(land(i,j))  % Plot non-significant correlations first
                scatter(abs(SAM_SAT_regr(i,j)),abs(n34_SAT_regr(w,i,j)),'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
            end           
        end
    end
end
% To make sure significant corrs are not covered, we need a separate
% loop...

for w = [1 2 3]
    subplot(2,3,w)
    ylim([0 1]); xlim([0 1])
    for i = 1:size(SAT_corr_sig,2)
        for j = 1:size(SAT_corr_sig,3)
            if ~isnan(SAT_corr_sig(w,i,j)) && ~isnan(land(i,j))  % Plot non-significant correlations first
                scatter(abs(SAM_SAT_regr(i,j)),abs(n34_SAT_regr(w,i,j)),'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
            end           
        end
    end
end

% Add 0 lines
% subplot(2,3,1)
% line([0 0],[-1 1.5],'Color',[17 17 17]/255,'LineStyle','--')
% line([-6 4],[0 0],'Color',[17 17 17]/255,'LineStyle','--')
% ylim([-1 1.5]); xlim([-6 4])
% subplot(2,3,2)
% line([0 0],[-1 1.5],'Color',[17 17 17]/255,'LineStyle','--')
% line([-6 4],[0 0],'Color',[17 17 17]/255,'LineStyle','--')
% xlim([-6 4])
% subplot(2,3,3)
% line([0 0],[-4 8],'Color',[17 17 17]/255,'LineStyle','--')
% line([-6 4],[0 0],'Color',[17 17 17]/255,'LineStyle','--')
% xlim([-6 4])

% Precip
for w = [1 2 3]
    subplot(2,3,w+3)
    ylim([0 1]); xlim([0 1])
    hold on
    xlabel('SAM-precip regr. coeff.')
    ylabel('n34-SAM/precip running corr regr. coeff')
    for i = 1:size(pre_corr_sig,2)
        for j = 1:size(pre_corr_sig,3)
            if isnan(pre_corr_sig(w,i,j)) && ~isnan(land(i,j))
                scatter(abs(SAM_precip_regr(i,j)),abs(n34_pre_regr(w,i,j)),'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
            end           
        end
    end
end
for w = [1 2 3]
    subplot(2,3,w+3)
    hold on
    for i = 1:size(pre_corr_sig,2)
        for j = 1:size(pre_corr_sig,3)
            if ~isnan(pre_corr_sig(w,i,j)) && ~isnan(land(i,j))
                scatter(abs(SAM_precip_regr(i,j)),abs(n34_pre_regr(w,i,j)),'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
            end           
        end
    end
end

% Add 0 lines
% subplot(2,3,4)
% line([0 0],[-1.5 1.5],'Color',[17 17 17]/255,'LineStyle','--')
% line([-1e-4 1.5e-4],[0 0],'Color',[17 17 17]/255,'LineStyle','--')
% xlim([-1e-4 1.5e-4])
% subplot(2,3,5)
% line([0 0],[-1.5 1.5],'Color',[17 17 17]/255,'LineStyle','--')
% line([-1e-4 1.5e-4],[0 0],'Color',[17 17 17]/255,'LineStyle','--')
% xlim([-1e-4 1.5e-4])
% subplot(2,3,6)
% line([0 0],[-6 6],'Color',[17 17 17]/255,'LineStyle','--')
% line([-1e-4 1.5e-4],[0 0],'Color',[17 17 17]/255,'LineStyle','--')
% xlim([-1e-4 1.5e-4])

h=gcf;
set(h,'PaperOrientation','landscape');
print(gcf,'figures/zscore_ENSO_absregr_scatter','-dpdf','-bestfit')


%% Alternatively, we produce a scatter plot comparing ENSO-SAM/SAT corrs with
%  corr. variance

load('DataFiles/ENSO_corrs.mat','SAT_corr_sig','pre_corr_sig')
sat_std = ncread('DataFiles/corr_std.nc','std_sat_wdw');
pre_std = ncread('DataFiles/corr_std.nc','std_precip_wdw');
land = ncread('DataFiles/sftlf_A1.static.nc','sftlf')';

figure(2)

for w = [1 2 3]
    subplot(2,3,w)
    hold on
    xlabel('SAM-SAT std. corr.')
    ylabel('n34-SAM/SAT running corr r')
    xlim([0 0.8])
    ylim([-1 1])
    for i = 1:size(SAT_corr_sig,2)
        for j = 1:size(SAT_corr_sig,3)
            if ~isnan(SAT_corr_sig(w,i,j)) && ~isnan(land(i,j))
                scatter(sat_std(i,j,w),SAT_corr_sig(w,i,j),'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
            end           
        end
    end
end

for w = [1 2 3]
    subplot(2,3,w+3)
    hold on
    xlabel('SAM-precip std. corr.')
    ylabel('n34-SAM/precip running corr r')
    xlim([0 0.8])
    ylim([-1 1])
    for i = 1:size(pre_corr_sig,2)
        for j = 1:size(pre_corr_sig,3)
            if ~isnan(pre_corr_sig(w,i,j)) && ~isnan(land(i,j))
                scatter(pre_std(i,j,w),pre_corr_sig(w,i,j),'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
            end           
        end
    end
end

h=gcf;
set(h,'PaperOrientation','landscape');
print(gcf,'figures/ENSO_corr_SAM_std','-dpdf','-bestfit')


%% Orrrr, ENSO regressions vs SAM-proxy correlations

clear
load('DataFiles/ENSO_corrs.mat','SAT_corr_sig','pre_corr_sig')
load('DataFiles/ENSO_regression.mat','n34_SAT_regr','n34_pre_regr')
load('DataFiles/SAM_corrs.mat','model_pre_corr_land','model_SAT_corr_land')
land = ncread('DataFiles/sftlf_A1.static.nc','sftlf')';

figure(1)
% SAT
nonsig_count = zeros(3);
for w = [1 2 3]
    subplot(2,3,w)
    hold on
    xlabel('SAM-SAT r')
    ylabel('n34-SAM/SAT running corr regr. coeff')
    for i = 1:size(SAT_corr_sig,2)
        for j = 1:size(SAT_corr_sig,3)
            if isnan(SAT_corr_sig(w,i,j)) && ~isnan(land(i,j))  % Plot non-significant correlations first
                scatter(abs(model_SAT_corr_land(i,j)),abs(n34_SAT_regr(w,i,j)),'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
                nonsig_count(w,1) = nonsig_count(w,1) + 1;
            end           
        end
    end
end
% To make sure significant corrs are not covered, we need a separate
% loop...
sig_count = zeros(3);
for w = [1 2 3]
    subplot(2,3,w)
    for i = 1:size(SAT_corr_sig,2)
        for j = 1:size(SAT_corr_sig,3)
            if ~isnan(SAT_corr_sig(w,i,j)) && ~isnan(land(i,j))  % Plot non-significant correlations first
                scatter(abs(model_SAT_corr_land(i,j)),abs(n34_SAT_regr(w,i,j)),'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
                sig_count(w,1) = sig_count(w,1) + 1;
            end           
        end
    end
end
pcnt_sig_SAT = sig_count./(nonsig_count+sig_count);


% Precip
nonsig_count = zeros(3);
for w = [1 2 3]
    subplot(2,3,w+3)
    hold on
    xlabel('SAM-precip r')
    ylabel('n34-SAM/precip running corr regr. coeff')
    for i = 1:size(pre_corr_sig,2)
        for j = 1:size(pre_corr_sig,3)
            if isnan(pre_corr_sig(w,i,j)) && ~isnan(land(i,j))
                scatter(abs(model_pre_corr_land(i,j)),abs(n34_pre_regr(w,i,j)),'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
                nonsig_count(w,1) = nonsig_count(w,1) + 1;
            end           
        end
    end
end
sig_count = zeros(3);
for w = [1 2 3]
    subplot(2,3,w+3)
    hold on
    for i = 1:size(pre_corr_sig,2)
        for j = 1:size(pre_corr_sig,3)
            if ~isnan(pre_corr_sig(w,i,j)) && ~isnan(land(i,j))
                scatter(abs(model_pre_corr_land(i,j)),abs(n34_pre_regr(w,i,j)),'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
                sig_count(w,1) = sig_count(w,1) + 1;
            end           
        end
    end
end
pcnt_sig_precip = sig_count./(nonsig_count+sig_count);

h=gcf;
set(h,'PaperOrientation','landscape');
print(gcf,'figures/zscoreENSO_regr_SAM_r_scatter_2','-dpdf','-painters','-bestfit')










