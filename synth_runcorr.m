% This script will calculate running correlations between synthetic temp/prec
% and the Nino3.4 index, and saves them
% This script needs mexcdf to be installed in matlab to run, the
% function 'b2r', 'plotworld', and the folder DataFiles
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
        wind_synruncorr=NaN(size(nu_wind),'single');
        precip_synruncorr=NaN(size(nu_precip),'single');
    else
        load(['../../../../../../media/My Book/CM2.1/Synth_runcorr/',num2str(windowsize),'yrWindow/run',num2str(n),'syncorr.mat'])
    end

% Running Correlation of Temperature
    for i=S_bound:N_bound
        for j=W_bound:E_bound
            wind_synruncorr(:,i,j)=movingCorrelation([squeeze(nu_wind(:,i,j)),SAM],windowsize,2);
        end
    end

save(['../../../../../../media/My Book/CM2.1/Synth_runcorr/',num2str(windowsize),'yrWindow/run',num2str(n),'syncorr.mat'],'wind_synruncorr','precip_synruncorr','windowsize');
toc;
end

for n=1:1000
n
tic;
load(['../../../../../../media/My Book/CM2.1/Synth_Data/run',num2str(n),'syn.mat'])

% USE IF OLD FILES EXIST ALREADY
    if  ~exist(['../../../../../../media/My Book/CM2.1/Synth_runcorr/',num2str(windowsize),'yrWindow/run',num2str(n),'syncorr.mat'],'file')
        wind_synruncorr=NaN(size(nu_wind),'single');
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


% %% Gathering Data for PDFs
% 
% x_lon = 160; [~,x_ind]= min(abs(lon-x_lon));
% y_lat = -20; [~,y_ind]= min(abs(lat-y_lat));
% 
% load(['Synth_corr_',num2str(windowsize),'yrwin/run',num2str(1000),'syncorr.mat']);
% if sum(~isnan(pr_synruncorr(:,y_ind,x_ind)),1)==0
%     error('Selected spot has no data');
% end
% 
% spot_ts=zeros(1000,499); % 1000 runs, 499 years of corr data
% spot_pr=zeros(1000,499); % 1000 runs, 499 years of corr data
% 
% for n=1:1000 % This will take about 2.5 minutes
%     load(['Synth_corr_',num2str(windowsize),'yrwin/run',num2str(n),'syncorr.mat']);
%     spot_ts(n,:) = squeeze(ts_synruncorr(:,y_ind,x_ind));
%     spot_pr(n,:) = squeeze(pr_synruncorr(:,y_ind,x_ind));
% end
% 
% % Translating to Fisher Z scores
% 
% spot_ts = 0.5*log( (1+spot_ts)./(1-spot_ts) );
% spot_pr = 0.5*log( (1+spot_pr)./(1-spot_pr) );
% pr_runcorrz = 0.5*log( (1+pr_runcorr)./(1-pr_runcorr) );
% ts_runcorrz = 0.5*log( (1+ts_runcorr)./(1-ts_runcorr) );
% 
% %% Fitting Distributions and Obtaining percentiles
% 
% for n=33:499
%     [pr_muhat pr_sighat] = normfit(spot_pr(:,n));
%     if kstest2(spot_pr(:,n),normrnd(pr_muhat,pr_sighat,1,10000))
%         disp(['Warning: Normal Series fit is inappropriate at n=',num2str(n)]);
%     end
% end
% % Percentiles
% pr_pc_spot = prctile(spot_pr,[2.5,97.5]);
% ts_pc_spot = prctile(spot_ts,[2.5,97.5]);
% plot(ts_pc_spot'); hold on;
% plot(squeeze(ts_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3);
% hold off;
% 
% length(find(squeeze(pr_runcorrz(:,y_ind,x_ind))<pr_pc_spot(1,:)'|squeeze(pr_runcorrz(:,y_ind,x_ind))>pr_pc_spot(2,:)'))
% length(find(squeeze(ts_runcorrz(:,y_ind,x_ind))<ts_pc_spot(1,:)'|squeeze(ts_runcorrz(:,y_ind,x_ind))>ts_pc_spot(2,:)'))
% %% Plot of Z scores for several time series
% 
% plot(squeeze(pr_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3)
% hold on;
% plot(squeeze(spot_pr(1,:)),'r')
% plot(squeeze(spot_pr(2,:)),'y')
% plot(squeeze(spot_pr(3,:)),'b')
% plot(squeeze(spot_pr(4,:)),'g')
% plot(squeeze(spot_pr(5,:)),'m')
% hold off
% title(['Fisher Z scores of Runing correlations of GFDL data vs synthetic data at ', num2str(x_lon),'E, ',num2str(y_lat),'N']);
% ylabel('Fisher Z Correlation Coefficients');
% xlabel('Year after the final year of the running correlation period');
% legend('GFDL','Synth1','Synth2','Synth3','Synth4','Synth5');
% 
% plot(squeeze(ts_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3)
% hold on;
% plot(squeeze(spot_ts(1,:)),'r')
% plot(squeeze(spot_ts(2,:)),'y')
% plot(squeeze(spot_ts(3,:)),'b')
% plot(squeeze(spot_ts(4,:)),'g')
% plot(squeeze(spot_ts(5,:)),'m')
% hold off
% title(['Fisher Z scores of Runing correlations of GFDL data vs synthetic data at ', num2str(x_lon),'E, ',num2str(y_lat),'N']);
% ylabel('Fisher Z Correlation Coefficients');
% xlabel('Year after the final year of the running correlation period');
% legend('GFDL','Synth1','Synth2','Synth3','Synth4','Synth5');
% 
% %% Plotting Area with data
% 
% load Synth_corr_31yrwin/run1000syncorr.mat
% pcolor(lon,lat,squeeze(ts_synruncorr(33,:,:)))
% title('Correlation Coefficients of last synthetic run, first correlation period');
% shading flat
% plotworld
