% This script wil plot the number of nonstationarities in an area of the
% map from existing data of synthetic running correlations.
% This requires the synthetic correlations produced by synth_runcorr.m, but
% transformed by synth_pointform.m in order to produce a map of nonstationarities
% This script needs mexcdf to be installed in matlab to run, the
% function 'b2r', 'plotworld', and the folder DataFiles

%% Setup
clear
load DataFiles/model_output.mat

window = 31; % The running window in years

%% Fitting Distributions and Obtaining percentiles

nonstat_precipmap=nan(length(lat),length(lon));
nonstat_satmap=nan(length(lat),length(lon));
running_nonstat_precipmap=nan(499,length(lat),length(lon));
running_nonstat_satmap=nan(499,length(lat),length(lon));
nonstat_satmaprecord=zeros(length(n34_ind),length(lat),length(lon));
nonstat_precipmaprecord=zeros(length(n34_ind),length(lat),length(lon));
pr_pc = nan(2,length(n34_ind),length(lat),length(lon));
ts_pc = nan(2,length(n34_ind),length(lat),length(lon));
run_wdw = 30; % Years
load('DataFiles/runcorr',num2str(window),'yr.mat');
runcorr_wdw = 31; % Running correlation window in the runcorr.mat file

% Limits of box
S_lat = -90; N_lat = 90; W_lon = 0; E_lon = 360;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));

for i=W_bound:E_bound
    for j=S_bound:N_bound
        if exist(['/srv/ccrc/data34/z3372730/MATLAB/Synth_runcorr/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'],'file')
            % This takes 0.3 seconds per point
            load(['/srv/ccrc/data34/z3372730/MATLAB/Synth_runcorr/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat']);
            spotz_sat = 0.5*log( (1+spot_sat)./(1-spot_sat) ); % Fishers Z Score
            spotz_precip = 0.5*log( (1+spot_precip)./(1-spot_precip) );
            pr_runcorrz = 0.5*log( (1+pr_runcorr)./(1-pr_runcorr) );
            ts_runcorrz = 0.5*log( (1+ts_runcorr)./(1-ts_runcorr) );
            pr_pc_spotz = prctile(spotz_precip,[2.5,97.5]);
            ts_pc_spotz = prctile(spotz_sat,[2.5,97.5]);
            pr_pc(:,:,j,i) = prctile(spot_precip,[2.5,97.5]);
            ts_pc(:,:,j,i) = prctile(spot_sat,[2.5,97.5]);
            nonstat_precipmap(j,i)=length( find(squeeze(pr_runcorrz(:,j,i))<pr_pc_spotz(1,:)'|...
                                          squeeze(pr_runcorrz(:,j,i))>pr_pc_spotz(2,:)')   );
            nonstat_satmap(j,i)=length( find(squeeze(ts_runcorrz(:,j,i))<ts_pc_spotz(1,:)'|...
                                          squeeze(ts_runcorrz(:,j,i))>ts_pc_spotz(2,:)')   );
            for n=1:(length(ts_runcorrz(:,1,1))+1-run_wdw-runcorr_wdw)
                running_nonstat_precipmap(n,j,i) = length( find(squeeze( ...
                    pr_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))<pr_pc_spotz(1,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))'|...
                    squeeze(pr_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))>pr_pc_spotz(2,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))')   );
                running_nonstat_satmap(n,j,i) = length( find(squeeze( ...
                    ts_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))<ts_pc_spotz(1,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))'|...
                    squeeze(ts_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))>ts_pc_spotz(2,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))')   );
            end
            nonstat_satmaprecord(find(squeeze(ts_runcorrz(:,j,i))<ts_pc_spotz(1,:)'|...
                                          squeeze(ts_runcorrz(:,j,i))>ts_pc_spotz(2,:)'),j,i) = 1;
            nonstat_precipmaprecord(find(squeeze(pr_runcorrz(:,j,i))<pr_pc_spotz(1,:)'|...
                                          squeeze(pr_runcorrz(:,j,i))>pr_pc_spotz(2,:)'),j,i) = 1;
                                        
        else
            disp(['Data at ',num2str(lon(i)),'E',num2str(lat(j)),'N does not exist']);
        end
    end
end

save('DataFiles/nonstat_map.mat','nonstat_precipmap','nonstat_satmap', ...
     'pr_pc','ts_pc','running_nonstat_precipmap','running_nonstat_satmap',...
     'nonstat_satmaprecord','nonstat_precipmaprecord');
 
%% Plotting Map
clear
load('DataFiles/nonstat_map_31yrwdw.mat')
load('DataFiles/model_output.mat')
NUM_CONTOURS = 10;

subplot(2,1,1);
contourf(lon,lat,nonstat_precipmap,NUM_CONTOURS);
plotworld;
colorbar;
caxis([0, 200]);
colormap(flipud(hot(NUM_CONTOURS)));
title(['Num of nonstationary stations for prec, rcor window:',num2str(window),'yrs'])
subplot(2,1,2);
contourf(lon,lat,nonstat_satmap,NUM_CONTOURS);
plotworld;
caxis([0, 200]);
colorbar
colormap(flipud(hot(NUM_CONTOURS)));
title(['Num of nonstationary stations for temp, rcor window:',num2str(window),'yrs'])

% pcolor(lon,lat,double(nonstat_windmap)); plotworld; colorbar; caxis([0,200]); colormap(flipud(hot));
