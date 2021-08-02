%% Making new Figure 1/ Comparison of model and ERA data
% Create new running correlations
clear all
irange = 144;
jrange = 90;
windowsize = 36;

% Data pre-processing
% Do CM2.1 first
load('DataFiles/model_output.mat')

tic
for windowsize = [31 61 91]
    for i= 1:irange
        for j= 1:jrange
            [sat_CM_runcorr(:,i,j), ~, sat_CM_pvals(:,i,j)] = movingCorrelation([squeeze(sat_detr(:,j,i)),SAM],windowsize,2);
            [pre_CM_runcorr(:,i,j), ~, pre_CM_pvals(:,i,j)] = movingCorrelation([squeeze(precip_detr(:,j,i)),SAM],windowsize,2);
        end
        i
    end
    save(['DataFiles/CM2_',num2str(windowsize),'yr_runcorr.mat'],'sat_CM_runcorr','sat_CM_pvals','pre_CM_runcorr','pre_CM_pvals');
    clear sat_CM_runcorr pre_CM_runcorr sat_CM_pvals pre_CM_pvals
    windowsize
end
toc % This operation takes ~1hr22 minutes on a desktop

%save('DataFiles/CM2_ERA_runcorr.mat','sat_CM_runcorr','sat_CM_pvals','pre_CM_runcorr','pre_CM_pvals','-append'); % This is for the 36yr ERA window
clear sat_detr precip_detr SAM

% Now lets do the GISS data
load('/media/huiskamp/My Book/CM2.1/SAM/GISS_data/model_output.mat')

tic
for i= 1:irange
    for j= 1:jrange
        [sat_GISS_runcorr(:,i,j), ~ , sat_GISS_pvals(:,i,j)] = movingCorrelation([squeeze(satG_detr(:,j,i)),SAM_G],windowsize,2);
        [pre_GISS_runcorr(:,i,j), ~ , pre_GISS_pvals(:,i,j)] = movingCorrelation([squeeze(precipG_detr(:,j,i)),SAM_G],windowsize,2);
    end
    i
end
toc % This operation takes ~53 minutes on a desktop


save('DataFiles/CM2_GISS_runcorr.mat','sat_GISS_runcorr','sat_GISS_pvals','pre_GISS_runcorr','pre_GISS_pvals');

% 1. Does ERA data fall within the range of model variance?
clear all

load('DataFiles/CM2_ERA_runcorr.mat')
%load('/media/huiskamp/My Book/CM2.1/SAM/GISS_data/model_output.mat')
load('DataFiles/CM2_GISS_runcorr.mat','sat_GISS_runcorr','pre_GISS_runcorr');

era_sat = ncread('DataFiles/ERA_data_regrid.nc','t2m');
era_precip = ncread('DataFiles/ERA_data_regrid.nc','tp');
load('marshall_SAM.mat') % Ann, Aut, Win, Spr, Sum

for i = 1:size(era_sat,1)
    for j = 1:size(era_sat,2)
        era_sat_corr(i,j) = corr(Marshall_SAM(23:58,2),squeeze(era_sat(i,j,1:36)));
        era_precip_corr(i,j) = corr(Marshall_SAM(23:58,2),squeeze(era_precip(i,j,1:36)));
    end
end

%lat_era = double(ncread('DataFiles/ERA_data_regrid.nc','lat'));
%lon_era = double(ncread('DataFiles/ERA_data_regrid.nc','lon'));

irange = 144;
jrange = 45;

for i = 1:irange
    for j = 1:jrange
        if era_sat_corr(i,j) >= min(sat_CM_runcorr(:,i,j)) && era_sat_corr(i,j) <= max(sat_CM_runcorr(:,i,j))   
            CM2_fit_sat(i,j) = 1;
            else CM2_fit_sat(i,j) = 0;
        end
        if era_precip_corr(i,j) >= min(pre_CM_runcorr(:,i,j)) && era_precip_corr(i,j) <= max(pre_CM_runcorr(:,i,j))   
            CM2_fit_pre(i,j) = 1;
            else CM2_fit_pre(i,j) = 0;
        end
        if era_sat_corr(i,j) >= min(sat_GISS_runcorr(:,i,j)) && era_sat_corr(i,j) <= max(sat_GISS_runcorr(:,i,j))   
            GISS_fit_sat(i,j) = 1;
            else GISS_fit_sat(i,j) = 0;
        end
        if era_precip_corr(i,j) >= min(pre_GISS_runcorr(:,i,j)) && era_precip_corr(i,j) <= max(pre_GISS_runcorr(:,i,j))   
            GISS_fit_pre(i,j) = 1;
            else GISS_fit_pre(i,j) = 0;
        end
    end
end

save('model_fit.mat','CM2_fit_pre','CM2_fit_sat','GISS_fit_pre','GISS_fit_sat')

%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%
clear
load('DataFiles/model_output.mat','lat','lon')
load('DataFiles/SAM_corrs.mat','model_SAT_corr','model_pre_corr')
model_SAT_corr = double(model_SAT_corr); model_pre_corr = double(model_pre_corr);
load('DataFiles/GISS_corrs.mat')
load('DataFiles/GISS_model_output.mat')
load('model_fit.mat')
irange = 144;
jrange = 45;
% land = ncread('DataFiles/sftlf_A1.static.nc','sftlf');
for i = 1:4
    figure(i)
    axesm('MapProjection','pcarree','MapLatLimit',[-90 -0],'maplonlim',[0 360])
    framem;
    gridm;
    mlabel;
    plabel;
    tightmap
end

levels = -1:0.1:1;
levels_n = -1:0.1:-0.3;
levels_p = 0.3:0.1:1;

% model SAT
figure(1)
contourfm(lat(1:45),lon,model_SAT_corr(1:45,:),levels,'linestyle','none')
hold on;
colormap(b2r(-1,1))
for j = 1:jrange
    for i = 1:irange
        if CM2_fit_sat(i,j) ~= 1
            scatterm(lat(j),lon(i),6,CM2_fit_sat(i,j),'k','filled')
        %else
        %    scatterm(lat(j),lon(i),6,CM2_fit_sat(i,j),'g','filled')
        end
    end
end

% model Precip
figure(2)
contourfm(lat(1:45),lon,model_pre_corr(1:45,:),levels,'linestyle','none')
hold on;
colormap(b2r(-1,1))
for i = 1:irange
    for j = 1:jrange
        if CM2_fit_pre(i,j) ~= 1
            scatterm(lat(j),lon(i),6,CM2_fit_pre(i,j),'k','filled')
        %else
        %    scatterm(lat(j),lon(i),6,CM2_fit_pre(i,j),'g','filled')
        end
    end
end

% GISS SAT
figure(3)
contourfm(latG(1:45),lonG,GISS_SAT_corr(1:45,:),levels,'linestyle','none')
hold on;
colormap(b2r(-1,1))
for i = 1:irange
    for j = 1:jrange
        if GISS_fit_sat(i,j) ~= 1
            scatterm(latG(j),lonG(i),6,GISS_fit_sat(i,j),'k','filled')
        %else
        %    scatterm(latG(j),lonG(i),6,GISS_fit_sat(i,j),'g','filled')
        end
    end
end

% GISS precip
figure(4)
contourfm(latG(1:45),lonG,GISS_pre_corr(1:45,:),levels,'linestyle','none')
hold on;
colormap(b2r(-1,1))
for i = 1:irange
    for j = 1:jrange
        if GISS_fit_pre(i,j) ~= 1
            scatterm(latG(j),lonG(i),6,GISS_fit_pre(i,j),'k','filled')
        %else
        %    scatterm(latG(j),lonG(i),6,GISS_fit_pre(i,j),'g','filled')
        end
    end
end

% Add continents

for i = 1:4
    figure(i)
    load coast
    plotm(lat,long,'k','linewidth',2)
end

figure(1)
print('CM2p1_sat_fit.pdf','-dpdf','-painters')%,'-bestfit')
figure(2)
print('CM2p1_precip_fit.pdf','-dpdf','-painters')%,'-bestfit')
figure(3)
print('GISS_sat_fit.pdf','-dpdf','-painters')%,'-bestfit')
figure(4)
print('GISS_precip_fit.pdf','-dpdf','-painters')%,'-bestfit')

% What % of sites fit within the model range?
% Total # of SH land cells  = 2241 
num_land = sum(sum(~isnan(land(:,1:45))));

for i = 1:144
    for j = 1:45
        if isnan(land(i,j))
            CM2_fit_sat(i,j) = 0;
            CM2_fit_pre(i,j) = 0;
        end
    end
end

sat_fit_pcnt = sum(sum(CM2_fit_sat))/num_land*100;
pre_fit_pcnt = sum(sum(CM2_fit_pre))/num_land*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make scatter plots of r vs pval for each of the running correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For land cells only...
clear all
land = ncread('DataFiles/sftlf_A1.static.nc','sftlf')';
load('DataFiles/model_output.mat')
irange = 144;
jrange = 45; % Only do the SH to save time!
nsim = 1000;

% Create 1k synthetic, red noise time-series for each proxy.
% Correlate each synthetic proxy with each SAM index
%

tic % This takes ~20 mins for only land cells, only has to be done once
for i = 1:irange
    for j = 1:jrange
        if ~exist(['noise_proxies/noise_prx_i',num2str(i),'_j',num2str(j),'.mat'],'file') %~isnan(land(j,i))
            value = rand(1); 
            noise_sat(:,:) = ebisuzaki(sat_detr(:,j,i),nsim,value); % Create 1,000 synthetic time-series' of each proxy
            noise_pre(:,:) = ebisuzaki(pre_detr(:,j,i),nsim,value);
            save(['noise_proxies/noise_prx_i',num2str(i),'_j',num2str(j),'.mat'],'noise_sat','noise_pre')
%            save(['noise_proxies/noise_prx_i',num2str(i),'_j',num2str(j),'.mat'],'noise_pre','-append')
            clear noise_pre noise_sat
        end
    end
end
toc

% Calculate the r values for noise proxies ...
% To do this, run the r_noise_calc_* scripts
% The output is then used in the next section...

% Calculate P vals for running corrs

for windowsize = [91]%[31 61 91]
    p_sat_noise = nan(500,irange,jrange);
    p_pre_noise = nan(500,irange,jrange);
    load(['DataFiles/CM2_',num2str(windowsize),'yr_runcorr.mat'],'sat_CM_runcorr','pre_CM_runcorr');
    for i = 1:irange
        for j = 1:jrange
            if ~isnan(land(j,i))
                load(['noise_proxies/noise_rvals',num2str(windowsize),'yrWDW_i',num2str(i),'_j',num2str(j),'.mat'],'r_sat','r_pre')
                for k = floor(windowsize/2):500-floor(windowsize/2)
                    if sat_CM_runcorr(k,i,j) > 0
                        p_sat_noise(k,i,j) = 1-(sum(sat_CM_runcorr(k,i,j) > r_sat(k,:)))/nsim;
                    elseif sat_CM_runcorr(k,i,j) < 0
                        p_sat_noise(k,i,j) = 1-(sum(sat_CM_runcorr(k,i,j) < r_sat(k,:)))/nsim;
                    end
                    if pre_CM_runcorr(k,i,j) > 0
                        p_pre_noise(k,i,j) = 1-(sum(pre_CM_runcorr(k,i,j) > r_pre(k,:)))/nsim;
                    elseif pre_CM_runcorr(k,i,j) < 0
                        p_pre_noise(k,i,j) = 1-(sum(pre_CM_runcorr(k,i,j) < r_pre(k,:)))/nsim;
                    end
                end
            end
        end
    i    
    end
    save(['DataFiles/CM2_',num2str(windowsize),'yr_runcorr.mat'],'p_sat_noise','p_pre_noise','-append')
    windowsize
end

% Calculate P vals for 500 year window
clear all
irange = 144;
jrange = 45;
nsim = 1000;
p_sat_noise = nan(irange,jrange);
p_pre_noise = nan(irange,jrange);
load('DataFiles/model_output.mat','SAM','sat_detr','precip_detr','lat','lon')
% Calculate 500 year correlation fields
for i = 1:size(sat_detr,2)
    for j = 1:size(sat_detr,3)
        model_SAT_corr(j,i) = double(corr(SAM,squeeze(sat_detr(:,i,j))));
        model_pre_corr(j,i) = double(corr(SAM,squeeze(precip_detr(:,i,j))));
    end
end
% Calculate P vals
for i = 1:irange
        for j = 1:jrange
            load(['noise_proxies/noise_rvals500yrWDW_i',num2str(i),'_j',num2str(j),'.mat'],'r_sat','r_pre')
            if model_SAT_corr(i,j) > 0
                p_sat_noise(i,j) = 1-(sum(model_SAT_corr(i,j) > r_sat(:)))/nsim;
            elseif model_SAT_corr(i,j) < 0
                p_sat_noise(i,j) = 1-(sum(model_SAT_corr(i,j) < r_sat(:)))/nsim;
            end
            if model_pre_corr(i,j) > 0
                p_pre_noise(i,j) = 1-(sum(model_pre_corr(i,j) > r_pre(:)))/nsim;
            elseif model_pre_corr(i,j) < 0
                p_pre_noise(i,j) = 1-(sum(model_pre_corr(i,j) < r_pre(:)))/nsim;
            end
        end
end
save('DataFiles/pvals_500yr_corrs.mat','p_sat_noise','p_pre_noise','model_SAT_corr','model_pre_corr')

figure(1)
% For SAT
for windowsize = [31 61 91]
    load(['DataFiles/CM2_',num2str(windowsize),'yr_runcorr.mat'])
    sat_CM_runcorr(:,isnan(land)') = nan;
    subplot(2,3,floor(windowsize/30))
    hold on
    for i = floor(windowsize/2):500-floor(windowsize/2)
        tmp_r = abs(reshape(sat_CM_runcorr(i,:,1:jrange),1,(irange*jrange)));
        tmp_p = reshape(p_sat_noise(i,:,:),1,(irange*jrange));
        plot(tmp_r,tmp_p,'o','linestyle','none','Color','k','MarkerFaceColor','r','MarkerSize',0.5)
    end
    xlim([0 1])
    ylim([0 0.1])
    line([0.3 0.3],[0 0.1],'linestyle','--','linewidth',1.5)
    %xlabel('r val'); ylabel('p val')
end
% For precip
for windowsize = [31 61 91]
    load(['DataFiles/CM2_',num2str(windowsize),'yr_runcorr.mat'])
    pre_CM_runcorr(:,isnan(land)') = nan;
    subplot(2,3,floor(windowsize/30)+3)
    hold on
    for i = floor(windowsize/2):500-floor(windowsize/2)
        tmp_r = abs(reshape(pre_CM_runcorr(i,:,1:jrange),1,(irange*jrange)));
        tmp_p = reshape(p_pre_noise(i,:,:),1,(irange*jrange));
        plot(tmp_r,tmp_p,'o','linestyle','none','Color','k','MarkerFaceColor','r','MarkerSize',0.5)
    end
    xlim([0 1])
    ylim([0 0.1])
    line([0.3 0.3],[0 0.1],'linestyle','--','linewidth',1.5)
    %xlabel('r val'); ylabel('p val')
end

subplot(2,3,1); title('31yr window'); ylabel('p val')
subplot(2,3,2); title('61yr window')
subplot(2,3,3); title('91yr window')
subplot(2,3,4); ylabel('p val')
subplot(2,3,5); xlabel('r val')

% print('Rs_vs_Ps.pdf','-painters','-dpdf') % Don't do this - it will lock
% up
print('Rs_vs_Ps.png','-dpng','-r300');  % this takes several minutes, be patient.

% Plotting r's vs p's for 500yr correlations
clear all
load('DataFiles/pvals_500yr_corrs.mat')
figure(2)
subplot(1,2,1)
hold on
% All correlations > 0.07 are significant at p < 0.1
tmp_r = abs(reshape(model_SAT_corr(:,1:45),1,(144*45)));
tmp_p = reshape(p_sat_noise,1,(144*45));
plot(tmp_r,tmp_p,'o','linestyle','none','Color','k','MarkerFaceColor','r','MarkerSize',0.5)
xlim([0 1])
ylim([0 0.1])
grid on
title('SAT'); ylabel('p val'); xlabel('r val')

subplot(1,2,2)
hold on
% All correlations > 0.07-0.08 are significant at p < 0.1
tmp_r = abs(reshape(model_pre_corr(:,1:45),1,(144*45)));
tmp_p = reshape(p_pre_noise,1,(144*45));
plot(tmp_r,tmp_p,'o','linestyle','none','Color','k','MarkerFaceColor','r','MarkerSize',0.5)
xlim([0 1])
ylim([0 0.1])
grid on
title('Precip.'); xlabel('r val')

print('Rs_vs_Ps_500.png','-dpng','-r300');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% # Non-stat proxies vs. Recon skill %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
NUM_YRS = 500;
NUM_TRIALS = 1000;
for window = [31 61 91]
    load(['DataFiles/nonstat_map_',num2str(window),'yrwdw.mat'])
    i=1; NUM_CAL_WDW = 10; 
    clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
    for c=0:9
        CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window);
    end
    for group_size=[2,5,10,30,50,70]
        nstat_nstns = nan(NUM_CAL_WDW,NUM_TRIALS);
        for c=1:size(CAL_WDW,1)
            DIR_NAME = ['Proxies/NoResample/',num2str(window),'yrWindow/'];
            load([DIR_NAME,'CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/',num2str(group_size),'stns_1000prox.mat'], ...
             'sat_lat','sat_lon','stn_lon','stn_lat');

            nstat_yrs_sat = nan(size(sat_lat));
            nstat_yrs_pre = nan(size(stn_lat));
            for m=1:group_size
                for tr=1:NUM_TRIALS
                    nstat_yrs_sat(tr,m) = nonstat_satmap(sat_lat(tr,m),sat_lon(tr,m));
                    nstat_yrs_pre(tr,m) = nonstat_precipmap(stn_lat(tr,m),stn_lon(tr,m));
                end
            end

            nstat_nstns_sat(c,:) = sum(nstat_yrs_sat > ceil(0.1*(NUM_YRS-window)),2);
            nstat_nstns_pre(c,:) = sum(nstat_yrs_pre > ceil(0.1*(NUM_YRS-window)),2);
            save(['DataFiles/num_nstat',num2str(window),'yrWdw_',num2str(group_size),'prox.mat'],'nstat_nstns_sat','nstat_nstns_pre')
        end
        %xbins = (0:1.0/group_size:1)-0.0000001; xbins(1)=0; xbins(end)=1;
        %h=histc(nstat_nstns(:),xbins)/100;
        %hold on; Hnd(i) = plot(xbins,h,'Color',cmap(i,:),'LineWidth',2);
        %h(~h) = NaN;
        %plot(xbins,h,'o','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:),'MarkerSize',7);
        %i=i+1;
    end
end

% Plotting
clear all
NUM_YRS = 500;
NUM_GROUPS=6;
NUM_CAL_WDW = 10;
cmap = hsv(NUM_GROUPS);
slopes = nan(2,6,3); % 2 variables, 6 groups, 3 windows

for windowsize = [31 61 91]
    for g = 1:NUM_GROUPS
        if g == 1
            group_size = 2;
        elseif g == 2
            group_size = 5;
        elseif g == 3
            group_size = 10; 
        elseif g == 4
            group_size = 30;
        elseif g == 5
            group_size = 50;
        elseif g == 6
            groupsize = 70;
        end
        clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
        for c = 0:9
            CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
        end
        for c = 1:size(CAL_WDW,1)
            load(['Proxies/NoResample/',num2str(windowsize),'yrWindow/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
                'all_sat_corr_CPS','all_stn_corr_CPS')
            sat_rvals(c,:) = all_sat_corr_CPS(group_size,:);
            stn_rvals(c,:) = all_stn_corr_CPS(group_size,:);
        end
        load(['DataFiles/num_nstat',num2str(windowsize),'yrWdw_',num2str(group_size),'prox.mat'])
        sat_r = reshape(sat_rvals,1,10*1000);
        pre_r = reshape(stn_rvals,1,10*1000);
        nstat_sat = reshape(nstat_nstns_sat,1,10*1000);
        nstat_pre = reshape(nstat_nstns_pre,1,10*1000);
        slopes(1,g,floor(windowsize/30)) = (corr2(sat_r,nstat_sat))^2;
        slopes(2,g,floor(windowsize/30)) = (corr2(pre_r,nstat_pre))^2;
        
        % Plotting time
        figure(1)
        idx = (floor(windowsize/30)-1)*12;
        subplot(6,6,idx+g)
        Fit = polyfit(sat_r,nstat_sat,1);
        x1 = linspace(0, 1, 1000);
        y1 = polyval(Fit,x1);
        plot(sat_r,nstat_sat,'o','linestyle','none','color',cmap(g,:),'MarkerFaceColor',cmap(g,:));
        hold on
        plot(x1,y1,'Color','k','linewidth',2)
        xlim([0 1]); 
        if g == 1
            ylim([0 2])
        elseif g == 2
            ylim([0 5])
        elseif g == 3 
            ylim([0 10])
        elseif g == 4 || g == 5 || g == 6
            ylim([0 20])
        end
        % Now for Precip
        subplot(6,6,idx+g+6)
        Fit = polyfit(pre_r,nstat_pre,1);
        x1 = linspace(0, 1, 1000);
        y1 = polyval(Fit,x1);
        plot(pre_r,nstat_pre,'o','linestyle','none','color',cmap(g,:),'MarkerFaceColor',cmap(g,:));
        hold on
        plot(x1,y1,'Color','k','linewidth',2)
        xlim([0 1]);
        if g == 1
            ylim([0 2])
        elseif g == 2
            ylim([0 5])
        elseif g == 3 
            ylim([0 10])
        elseif g == 4 || g == 5 || g == 6
            ylim([0 20])
        end
    end
end
            
print('Rs_vs_nstat.pdf','-painters','-dpdf')
        
   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Alternative figure - linegraph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
NUM_YRS = 500;
NUM_TRIALS = 1000;
num_prox = 70;
NUM_CAL_WDW = 10;

for window = [31 61 91]
    DIR_NAME = ['Proxies/NoResample/',num2str(window),'yrWindow/'];
    load(['DataFiles/nonstat_map_',num2str(window),'yrwdw.mat'])
    clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
    for c=0:9
        CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window);
    end
    for group_size= 2:num_prox
        for c=1:NUM_CAL_WDW
            load([DIR_NAME,'CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/',num2str(group_size),'stns_1000prox.mat'], ...
             'sat_lat','sat_lon','stn_lon','stn_lat');
            nstat_yrs_sat = nan(size(sat_lat));
            nstat_yrs_pre = nan(size(stn_lat));
            for m=1:group_size
                for tr=1:NUM_TRIALS
                    nstat_yrs_sat(tr,m) = nonstat_satmap(sat_lat(tr,m),sat_lon(tr,m));
                    nstat_yrs_pre(tr,m) = nonstat_precipmap(stn_lat(tr,m),stn_lon(tr,m));
                end
            end

            nstat_nstns_sat(c,:) = sum(nstat_yrs_sat > ceil(0.1*(NUM_YRS-window)),2);
            nstat_nstns_pre(c,:) = sum(nstat_yrs_pre > ceil(0.1*(NUM_YRS-window)),2);
            save(['DataFiles/num_nstat',num2str(window),'yrWdw_',num2str(group_size),'prox.mat'],'nstat_nstns_sat','nstat_nstns_pre')
        end
    end
end

%%%% calculate correlations between recon skill and nstat %. Consolidate into windowsizes

clear all
NUM_YRS = 500;
num_prox = 70;
NUM_CAL_WDW = 10;
nstat_corrs = nan(3,2,69); % windows, vars, network size
regr_slopes = nan(3,2,69);
for windowsize = [31 61 91]
    for group_size = 2:num_prox
        clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
        for c = 0:9
            CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
        end
        for c = 1:size(CAL_WDW,1)
            load(['Proxies/NoResample/',num2str(windowsize),'yrWindow/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
                'all_sat_corr_CPS','all_stn_corr_CPS')
            sat_rvals(c,:) = all_sat_corr_CPS(group_size,:);
            stn_rvals(c,:) = all_stn_corr_CPS(group_size,:);
        end
        load(['DataFiles/num_nstat',num2str(windowsize),'yrWdw_',num2str(group_size),'prox.mat'])
        sat_r = reshape(sat_rvals,1,10*1000);
        pre_r = reshape(stn_rvals,1,10*1000);
        nstat_sat = reshape(nstat_nstns_sat,1,10*1000)./group_size;
        nstat_pre = reshape(nstat_nstns_pre,1,10*1000)./group_size;
        [Rs1 Ps1] = corrcoef(sat_r,nstat_sat);
        [Rs2 Ps2] = corrcoef(pre_r,nstat_pre);
        nstat_corrs(floor(windowsize/30),1,group_size) = Rs1(1,2); nstat_Ps(floor(windowsize/30),1,group_size) = Ps1(1,2);
        nstat_corrs(floor(windowsize/30),2,group_size) = Rs2(1,2); nstat_Ps(floor(windowsize/30),2,group_size) = Ps2(1,2);
        p_sat = polyfit(sat_r,nstat_sat,1); p_pre = polyfit(pre_r,nstat_pre,1);
        regr_slopes(floor(windowsize/30),1,group_size) = p_sat(1);
        regr_slopes(floor(windowsize/30),2,group_size) = p_pre(1);
    end
end
        
% Plotting - I think we're better off going for a line graph here...

figure(1)
subplot(1,2,1)
hold on 
grid on
% SAT
plot(squeeze(nstat_corrs(1,1,:)),'o','color',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(nstat_corrs(2,1,:)),'o','color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(nstat_corrs(3,1,:)),'o','color',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'markersize',3,'linestyle','-','linewidth',2)     
% Precip
plot(squeeze(nstat_corrs(1,2,:)),'o','color',[0.4940, 0.1840, 0.5560],'MarkerFaceColor',[0.4940, 0.1840, 0.5560],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(nstat_corrs(2,2,:)),'o','color',[0.4660, 0.6740, 0.1880],'MarkerFaceColor',[0.4660, 0.6740, 0.1880],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(nstat_corrs(3,2,:)),'o','color',[0.3010, 0.7450, 0.9330],'MarkerFaceColor',[0.3010, 0.7450, 0.9330],'markersize',3,'linestyle','-','linewidth',2)
line([0 70],[0 0],'linestyle','--','color','k') % 0 line
line([0 70],[-0.02 -0.02],'color','r'); line([0 70],[0.027 0.027],'color','r') % region where corrs are not significant at p < 0.05
legend('SAT-31yrWdw','SAT-61yrWdw','SAT-91yrWdw','Pre-31yrWdw','Pre-61yrWdw','Pre-91yrWdw')
xlabel('Network Size'); ylabel('r');

% Plot regression slopes...
subplot(1,2,2)
hold on 
grid on
% SAT
plot(squeeze(regr_slopes(1,1,:)),'o','color',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(regr_slopes(2,1,:)),'o','color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(regr_slopes(3,1,:)),'o','color',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'markersize',3,'linestyle','-','linewidth',2)     
% Precip
plot(squeeze(regr_slopes(1,2,:)),'o','color',[0.4940, 0.1840, 0.5560],'MarkerFaceColor',[0.4940, 0.1840, 0.5560],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(regr_slopes(2,2,:)),'o','color',[0.4660, 0.6740, 0.1880],'MarkerFaceColor',[0.4660, 0.6740, 0.1880],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(regr_slopes(3,2,:)),'o','color',[0.3010, 0.7450, 0.9330],'MarkerFaceColor',[0.3010, 0.7450, 0.9330],'markersize',3,'linestyle','-','linewidth',2)
line([0 70],[0 0],'linestyle','--','color','k') % 0 line
line([0 70],[-0.02 -0.02],'color','r'); line([0 70],[0.027 0.027],'color','r') % region where corrs are not significant at p < 0.05
ylabel('regression slope');


print('nonstat_impact.pdf','-dpdf','-painters')%,'-bestfit')

%% Calculate the % of land points that are non-stationary
num_nstat = nan(2,3); pct_nstat = nan(2,3);
land = ncread('DataFiles/sftlf_A1.static.nc','sftlf')';
for window = [31 61 91]
    load(['DataFiles/nonstat_map_',num2str(window),'yrwdw.mat'])

    nonstat_precipmap(nonstat_precipmap<((500-(window-1))/10)) = nan;
    nonstat_satmap(nonstat_satmap<((500-(window-1))/10)) = nan;
    for i = 1:144
        for j = 1:45
            if isnan(land(j,i))
                nonstat_precipmap(j,i) = nan;
                nonstat_satmap(j,i) = nan;
            end
        end
    end
    num_nstat(1,floor(window/30)) =  sum(sum(~isnan(nonstat_satmap)));
    num_nstat(2,floor(window/30)) =  sum(sum(~isnan(nonstat_precipmap)));
    pct_nstat(1,floor(window/30)) =  sum(sum(~isnan(nonstat_satmap)))/2241*100;
    pct_nstat(2,floor(window/30)) =  sum(sum(~isnan(nonstat_precipmap)))/2241*100;
end







