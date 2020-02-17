%% ENSO correlation maps
% This script, based on the ENSO indices calculated in 'ENSO_stuff.m',
% correlates these indices with the running correlation between SAM and
% SAT/Precip at each grid cell. 
%
% December, 2019

% Setup
clear
load('DataFiles/model_output.mat','lat');
load('DataFiles/filtered_n34.mat');
land = ncread('DataFiles/sftlf_A1.static.nc','sftlf')';
lat_2 = 45; % change this from 45 to get some sub-set of the SH
n34_SAM_SAT = nan(3,lat_2,144); n34_SAM_pre = nan(3,lat_2,144);
n34 = nan(3,499); n34(1,:) = filt_ann_30; n34(2,:) = filt_ann_60;
n34(3,:) = filt_ann_90;

    
for windowsize = [31 61 91]
    load(['DataFiles/runcorr',num2str(windowsize),'yrwdw.mat']);
%     for i = 1:size(sat_runcorr,1)
%         sat_runcorr(i,isnan(land)) = nan;
%         precip_runcorr(i,isnan(land)) = nan;
%     end
    range = ceil(windowsize/2):500-floor(windowsize/2);
    sat_runcorr = sat_runcorr(range,1:lat_2,:); % we only need the SH and eliminate the NaNs from the running correlation at either end.
    pre_runcorr = precip_runcorr(range,1:lat_2,:);
    
    for i = 1:size(sat_runcorr,2)
        for j = 1:size(sat_runcorr,3)
            if ~isnan(sat_runcorr(1,i,j))
                [n34_SAM_SAT(floor(windowsize/30),i,j),corr_sig_SAT(floor(windowsize/30),i,j),~] = ...
                    calc_statsig(sat_runcorr(:,i,j),n34(floor(windowsize/30),range)');
                [n34_SAM_pre(floor(windowsize/30),i,j),corr_sig_pre(floor(windowsize/30),i,j),~] = ...
                    calc_statsig(pre_runcorr(:,i,j),n34(floor(windowsize/30),range)');
            end
        end
    end
    
end
    
% This plots the first 3 EOFs of the running correlation between SAT and SAM
levels = -1:0.05:1;
levels2 = -1:0.01:1;
levels_p = 0.005:0.1:1;
levels_sig = 0.95:0.011:1;
levels_n = -1:0.1:0;
reg_lat = 1:45;
reg_lon = 1:144;

% Mask out regions of significant correlation
SAT_corr_sig = nan(size(n34_SAM_SAT));
pre_corr_sig = nan(size(n34_SAM_pre));
for i = 1:3
    [r, c] = find(squeeze(corr_sig_SAT(i,:,:))>90); % Ensure we're getting p = 95 or higher
    [r2,c2]= find(squeeze(corr_sig_pre(i,:,:))>90);
    for j = 1:size(c,1)
        SAT_corr_sig(i,r(j),c(j)) = n34_SAM_SAT(i,r(j),c(j));
    end
    for k= 1:size(c2,1)
        pre_corr_sig(i,r2(k),c2(k)) = n34_SAM_pre(i,r2(k),c2(k));
    end
end    

% Plotting in matlab is a pile of ass, so lets save the output and do it in
% python.

save('DataFiles/ENSO_corrs.mat','lat','lon','n34_SAM_SAT','n34_SAM_pre','SAT_corr_sig','pre_corr_sig')

% Plot for SAT
subplot1(3, 1, 'Gap', [.01 .03], 'XTickL', 'Margin', 'YTickL', 'Margin');
for i = 1:3
    subplot1(i)
    axesm('MapProjection','eqdcylin','origin',[0,180],'MapLatLimit',[-90 0],'meridianlabel','on','parallellabel','on','PLabelMeridian',0,'Mlabelparallel','south')
    framem
    gridm %If you want gridlinesload coast
    load('DataFiles/model_output.mat','lat','lon')
    contourfm(lat(reg_lat),lon(reg_lon),squeeze(n34_SAM_SAT(i,:,:)),levels,'linestyle','none') 
    %contourm(lat(reg_lat),lon(reg_lon),squeeze(n34_SAM_SAT(i,:,:)),levels_p,'k'); contourm(lat(reg_lat),lon(reg_lon),squeeze(n34_SAM_SAT(i,:,:)),levels_n,'k','linestyle','--');
    colormap(b2r(-1,1))
    contourm(lat(reg_lat),lon(reg_lon),squeeze(corr_sig_SAT(i,:,:)),levels2,'g','linewidth',1.5)
    load coast
    plotm(lat,long,'k','linewidth',2)
    %title(['SAT-SAM running corr for',num2str(i*10+1),'yr Window - ',num2str(round(sat_var_exp(i,:))),'% var. exp.'])
    clear lat long
end
    
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('figures/n34_SAM_SAT_corrs2.pdf','-dpdf','-besfit')
    
% Plot for Precip

subplot1(3, 1, 'Gap', [.01 .03], 'XTickL', 'Margin', 'YTickL', 'Margin');
for i = 1:3
    subplot1(i)
    axesm('MapProjection','eqdcylin','origin',[0,180],'MapLatLimit',[-90 0],'meridianlabel','on','parallellabel','on','PLabelMeridian',0,'Mlabelparallel','south')
    framem
    gridm %If you want gridlinesload coast
    load('DataFiles/model_output.mat','lat','lon')
    contourfm(lat(reg_lat),lon(reg_lon),squeeze(n34_SAM_pre(i,:,:)),levels,'linestyle','none') 
    %contourm(lat(reg_lat),lon(reg_lon),squeeze(n34_SAM_pre(i,:,:)),levels_p,'k'); contourm(lat(reg_lat),lon(reg_lon),squeeze(n34_SAM_pre(i,:,:)),levels_n,'k','linestyle','--');
    colormap(b2r(-1,1))
    contourm(lat(reg_lat),lon(reg_lon),squeeze(corr_sig_pre(i,:,:)),levels2,'g','linewidth',1.5)
    load coast
    plotm(lat,long,'k','linewidth',2)
    %title(['SAT-SAM running corr for',num2str(i*10+1),'yr Window - ',num2str(round(sat_var_exp(i,:))),'% var. exp.'])
    clear lat long
end
    
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('figures/n34_SAM_pre_corrs2.pdf','-dpdf','-besfit')
    
    
    
%% Alternatively, we can check significance with Nerilie's 2014 method
% Create 1000 series' with same spectral properties and compare their
% correlation to our model series.
% NOTE: Not sure of the appropriate method to apply here: Would we create 
% the simulated TS for both the n34 AND the running corr records? They both
% have the same issue of reduced degrees of freedom, but at that point
% we're just correlating noise with noise - does this make sense?
clear
load('DataFiles/model_output.mat','lat');
load('DataFiles/filtered_n34.mat');
land = ncread('DataFiles/sftlf_A1.static.nc','sftlf')';
lat_2 = 45; % change this from 45 to get some sub-set of the SH
n34_SAM_SAT = nan(3,lat_2,144); n34_SAM_pre = nan(3,lat_2,144);
n34 = nan(3,499); n34(1,:) = filt_ann_30; n34(2,:) = filt_ann_60;
n34(3,:) = filt_ann_90;
    
    
    
    
    
nsim = 10000;
noise_prox = nan(995,52,nsim);
r_noise_prox = nan(3,52,nsim);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    