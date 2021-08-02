% This script will use existing running correlations between synthetic temp/prec
% and the SAM index (from synth_corr), modifying the format so each
% grid point can be analysed for percentiles and nonstationarities

%% Setup
load DataFiles/model_output.mat

window = 31; % The running window in years

% Loading 
% load('DataFiles/runcorr.mat');
% pr_runcorrz = pr_runcorr;
% pr_runcorrz = 0.5*log( (1+pr_runcorrz)./(1-pr_runcorrz) );
% ts_runcorrz = ts_runcorr;
% ts_runcorrz = 0.5*log( (1+ts_runcorrz)./(1-ts_runcorrz) );

%% Changing data organisation to spatial-point form
tic;
% Limits of box
S_lat = -90; N_lat = 0; W_lon = 0; E_lon = 360;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));

spot_sat = zeros(1000,500);
spot_precip = zeros(1000,500);
for i=W_bound:E_bound
    for j=S_bound:N_bound
        if exist(['Synth_runcorr/',num2str(window),'yrWindow/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'],'file')
            load(['Synth_runcorr/',num2str(window),'yrWindow/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'])
            if size(spot_precip,1)<1000
                recalc=true;
            else 
                recalc=false;
            end
        else
            recalc=true;
        end
        if recalc==true
            for n=1:1000 % This will take about 6 minutes
                load(['Synth_runcorr/run',num2str(n),'syncorr.mat']);
                spot_sat(n,:) = ts_synruncorr(:,j,i);
                spot_precip(n,:) = pr_synruncorr(:,j,i);
            end
            save(['Synth_runcorr/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'],...
                'spot_sat','spot_precip');
        end
        toc;
    end
end
toc;
