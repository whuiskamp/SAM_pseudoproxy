% This script will load in data and calculate the anomalous temperature and
% precipitation data and store it in a .mat file.

slp_file = '/p/projects/climber3/huiskamp/Other_projects/SAM_model/GISS_data/slp_ann.nc'; sat_file = '/p/projects/climber3/huiskamp/Other_projects/SAM_model/GISS_data/sat_ann.nc'; 
pr_file = '/p/projects/climber3/huiskamp/Other_projects/SAM_model/GISS_data/pr_ann.nc';

latG = ncread(slp_file,'lat');
lonG = ncread(slp_file,'lon');
time = ncread(slp_file,'time'); % Assumes both files use the same time

sat = ncread(sat_file,'tas');
precip = ncread(pr_file,'pr');
slp = ncread(slp_file,'psl');


%% Calculating the SAM index

% defined in Gallant et al. 2013
% as the difference between the normalized zonally
% averaged sea level pressure anomalies at 40S and 60S,

slp_40s = squeeze(slp(:,26,:)); % calculate SLP for 40S and 60S
slp_60s = squeeze(slp(:,16,:)); 
slp_40_mean = detrend(nanmean(slp_40s,1)); % take the zonal mean and detrend
slp_60_mean = detrend(nanmean(slp_60s,1));
SAM1 = slp_40_mean' - slp_60_mean';

for i = 1:length(SAM1)
   SAM2(i) = (SAM1(i)-min(SAM1))/(max(SAM1)-min(SAM1)); % normalise
end

SAM_G = zeros(851,1);
for i = 1:length(SAM_G)
    SAM_G(i) = SAM2(i) - mean(SAM2); % calculate anomalies
end

%%
% Removing trends

for i = 1:90
    for j = 1:144
       precipG_detr(:,i,j) = detrend(squeeze(precip(j,i,:)),'linear');
       satG_detr(:,i,j) = detrend(squeeze(sat(j,i,:)),'linear');
    end
end

save('GISS_data/model_output.mat','SAM_G','satG_detr','precipG_detr','latG','lonG');
