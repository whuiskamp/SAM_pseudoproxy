% Original scripts written by Ryan Batehup, 2014
% Edited for use with SAM by Willem Huiskamp, 2015

% This script will load in data and calculate the anomalous temperature and
% precipitation data and store it in a .mat file.

ua_file = '/Media/My Book/CM2.1/ua_A0.0001-0500.nc'; pr_file = '/Media/My Book/CM2.1/pr_A0.0001-0500.nc'; 
slp_file = '/Media/My Book/CM2.1/psl_A0.0001-0500.nc'; sat_file = 'DataFiles/ts_A0.0001-0500.nc';

lat = ncread(ua_file,'lat');
lon = ncread(ua_file,'lon');
time = ncread(ua_file,'time'); % Assumes both files use the same time
ua = ncread(ua_file,'ua');
sat = ncread(sat_file,'ts');
precip = ncread(pr_file,'pr');
slp = ncread(slp_file,'psl');

% wind = squeeze(ua(:,3,:,:)); %for 850hPa

%% Calculating the SAM index

% defined in Gallant et al. 2013
% as the difference between the normalized zonally
% averaged sea level pressure anomalies at 40S and 60S,

slp_40s = squeeze(slp(:,26,:)); % calculate SLP for 40S and 60S
slp_60s = squeeze(slp(:,16,:)); 
slp_40_mean = detrend(nanmean(slp_40s,2)); % take the zonal mean and detrend
slp_60_mean = detrend(nanmean(slp_60s,2));
SAM1 = slp_40_mean - slp_60_mean;

for i = 1:length(SAM1)
   SAM2(i) = (SAM1(i)-min(SAM1))/(max(SAM1)-min(SAM1)); % normalise
end

for i = 1:length(SAM)
    SAM(i) = SAM2(i) - mean(SAM2); % calculate anomalies
end

%%
% Removing trends

for i = 1:90
    for j = 1:144
       precip_detr(:,i,j) = detrend(squeeze(precip(j,i,:),'linear');
       sat_detr(:,i,j) = detrend(squeeze(sat(j,i,:)),'linear');
    end
end
% for i = 1:90
%     for j = 1:144
%        wind_detr(:,i,j) = detrend(wind(:,i,j),'linear');
%     end
% end

%save('model_output.mat','SAM','precip_detr','wind_detr','lat','lon');
save('DataFiles/model_output.mat','sat_detr','-append');

