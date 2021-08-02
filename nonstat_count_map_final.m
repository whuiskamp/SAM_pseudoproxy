% This script will plot the number of nonstationarities in an area of the
% map from existing data of synthetic running correlations.
% This requires the synthetic correlations produced by synth_runcorr.m, but
% transformed by synth_pointform.m in order to produce a map of nonstationarities

%% Setup
clear
load DataFiles/model_output.mat

window = 61; % The running window in years

%% Fitting Distributions and Obtaining percentiles

nonstat_precipmap=nan(length(lat),length(lon));
nonstat_satmap=nan(length(lat),length(lon));
running_nonstat_precipmap=nan(500,length(lat),length(lon));
running_nonstat_satmap=nan(500,length(lat),length(lon));
nonstat_satmaprecord=zeros(length(SAM),length(lat),length(lon));
nonstat_precipmaprecord=zeros(length(SAM),length(lat),length(lon));
pre_pc = nan(2,length(SAM),length(lat),length(lon));
sat_pc = nan(2,length(SAM),length(lat),length(lon));
run_wdw = 60; % Years
load(['DataFiles/runcorr',num2str(window),'yrwdw.mat']);
runcorr_wdw = 61; % Running correlation window in the runcorr.mat file

% Limits of box
S_lat = -90; N_lat = 0; W_lon = 0; E_lon = 360;
[~,S_bound]= min(abs(lat-S_lat));
[~,N_bound]= min(abs(lat-N_lat));
[~,W_bound]= min(abs(lon-W_lon));
[~,E_bound]= min(abs(lon-E_lon));

for i=W_bound:E_bound
    for j=S_bound:N_bound
        if exist(['Synth_runcorr/',num2str(window),'yrWindow/Point/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat'],'file')
            % This takes 0.3 seconds per point
            load(['Synth_runcorr/',num2str(window),'yrWindow/Point/',num2str(lon(i)),'E',num2str(lat(j)),'N_syncorr.mat']);
            spotz_sat = 0.5*log( (1+spot_sat)./(1-spot_sat) ); % Fishers Z Score
            spotz_precip = 0.5*log( (1+spot_precip)./(1-spot_precip) );
            pre_runcorrz = 0.5*log( (1+precip_runcorr)./(1-precip_runcorr) );
            sat_runcorrz = 0.5*log( (1+sat_runcorr)./(1-sat_runcorr) );
            pre_pc_spotz = prctile(spotz_precip,[2.5,97.5]);
            sat_pc_spotz = prctile(spotz_sat,[2.5,97.5]);
            pre_pc(:,:,j,i) = prctile(spot_precip,[2.5,97.5]);
            sat_pc(:,:,j,i) = prctile(spot_sat,[2.5,97.5]);
            nonstat_precipmap(j,i)=length( find(squeeze(pre_runcorrz(:,j,i))<pre_pc_spotz(1,:)'|...
                                          squeeze(pre_runcorrz(:,j,i))>pre_pc_spotz(2,:)')   );
            nonstat_satmap(j,i)=length( find(squeeze(sat_runcorrz(:,j,i))<sat_pc_spotz(1,:)'|...
                                          squeeze(sat_runcorrz(:,j,i))>sat_pc_spotz(2,:)')   );
            for n=1:(length(sat_runcorrz(:,1,1))+1-run_wdw-runcorr_wdw)
                running_nonstat_precipmap(n,j,i) = length( find(squeeze( ...
                    pre_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))<pre_pc_spotz(1,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))'|...
                    squeeze(pre_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))>pre_pc_spotz(2,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))')   );
                running_nonstat_satmap(n,j,i) = length( find(squeeze( ...
                    sat_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))<sat_pc_spotz(1,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))'|...
                    squeeze(sat_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))>sat_pc_spotz(2,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))')   );
            end
            nonstat_satmaprecord(find(squeeze(sat_runcorrz(:,j,i))<sat_pc_spotz(1,:)'|...
                                          squeeze(sat_runcorrz(:,j,i))>sat_pc_spotz(2,:)'),j,i) = 1;
            nonstat_precipmaprecord(find(squeeze(pre_runcorrz(:,j,i))<pre_pc_spotz(1,:)'|...
                                          squeeze(pre_runcorrz(:,j,i))>pre_pc_spotz(2,:)'),j,i) = 1;
                                        
        else
            disp(['Data at ',num2str(lon(i)),'E',num2str(lat(j)),'N does not exist']);
        end
    end
    i
end

save(['DataFiles/nonstat_map_',num2str(runcorr_wdw),'yrwdw.mat'],'nonstat_precipmap','nonstat_satmap', ...
     'pre_pc','sat_pc','running_nonstat_precipmap','running_nonstat_satmap',...
     'nonstat_satmaprecord','nonstat_precipmaprecord');
 

%% Plotting Map
% clear
% load('DataFiles/nonstat_map_31yrwdw.mat')
% load('DataFiles/model_output.mat')
% NUM_CONTOURS = 10;

% subplot(2,1,1);
% contourf(lon,lat,nonstat_precipmap,NUM_CONTOURS);
% plotworld;
% colorbar;
% caxis([0, 200]);
% colormap(flipud(hot(NUM_CONTOURS)));
% title(['Num of nonstationary stations for prec, rcor window:',num2str(window),'yrs'])
% subplot(2,1,2);
% contourf(lon,lat,nonstat_satmap,NUM_CONTOURS);
% plotworld;
% caxis([0, 200]);
% colorbar
% colormap(flipud(hot(NUM_CONTOURS)));
% title(['Num of nonstationary stations for temp, rcor window:',num2str(window),'yrs'])

% pcolor(lon,lat,double(nonstat_windmap)); plotworld; colorbar; caxis([0,200]); colormap(flipud(hot));
