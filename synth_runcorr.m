% This script will calculate running correlations between synthetic temp/prec
% and the SAM index, and saves them.
%% Setup

load DataFiles/model_output.mat

windowsize = 31; % The running window in years

%% Calculating running correlations from SYNTHETIC data
tic;
% Limits of box to calculate corr coefs
S_lat = -90; N_lat = 0; W_lon = 0; E_lon = 360;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));

for n=1:1000
n
tic;
load(['../../../../../../media/My Book/CM2.1/Synth_Data/run',num2str(n),'syn.mat'])

% USE IF OLD FILES EXIST ALREADY
    if  ~exist(['../../../../../../media/My Book/CM2.1/Synth_runcorr/',num2str(windowsize),'yrWindow/run',num2str(n),'syncorr.mat'],'file')
        sat_synruncorr=NaN(size(nu_sat),'single');
        precip_synruncorr=NaN(size(nu_precip),'single');
    else
        load(['../../../../../../media/My Book/CM2.1/Synth_runcorr/',num2str(windowsize),'yrWindow/run',num2str(n),'syncorr.mat'])
    end

% Running Correlation of Temperature
    for i=S_bound:N_bound
        for j=W_bound:E_bound
            sat_synruncorr(:,i,j)=movingCorrelation([squeeze(nu_sat(:,i,j)),SAM],windowsize,2);
        end
    end

save(['../../../../../../media/My Book/CM2.1/Synth_runcorr/',num2str(windowsize),'yrWindow/run',num2str(n),'syncorr.mat'],'sat_synruncorr','precip_synruncorr','windowsize');
toc;
end

for n=1:1000
n
tic;
load(['../../../../../../media/My Book/CM2.1/Synth_Data/run',num2str(n),'syn.mat'])

% USE IF OLD FILES EXIST ALREADY
    if  ~exist(['../../../../../../media/My Book/CM2.1/Synth_runcorr/',num2str(windowsize),'yrWindow/run',num2str(n),'syncorr.mat'],'file')
        sat_synruncorr=NaN(size(nu_sat),'single');
        precip_synruncorr=NaN(size(nu_precip),'single');
    else
        load(['../../../../../../media/My Book/CM2.1/Synth_runcorr/',num2str(windowsize),'yrWindow/run',num2str(n),'syncorr.mat'])
    end

% Running Correlation of Precipitation
    for i=S_bound:N_bound
        for j=W_bound:E_bound
            precip_synruncorr(:,i,j)=movingCorrelation([squeeze(nu_precip(:,i,j)),SAM],windowsize,2);
        end
    end

save(['../../../../../../media/My Book/CM2.1/Synth_runcorr/',num2str(windowsize),'yrWindow/run',num2str(n),'syncorr.mat'],'precip_synruncorr','-append'); %make sure the other one runs first
toc;
end
