% This script wil plot the number of nonstationarities in an area of the
% map from existing data of synthetic running correlations.
% This requires the synthetic correlations produced by synth_runcorr.m, but
% transformed by synth_pointform.m in order to produce a map of nonstationarities
% This script needs mexcdf to be installed in matlab to run, the
% function 'b2r', 'plotworld', and the folder DataFiles

%% Setup

load DataFiles/model_output.mat

window = 31; % The running window in years

%% Fitting Distributions and Obtaining percentiles

nonstat_precipmap=nan(length(lat),length(lon));
nonstat_windmap=nan(length(lat),length(lon));
running_nonstat_precipmap=nan(499,length(lat),length(lon));
running_nonstat_windmap=nan(499,length(lat),length(lon));
nonstat_windmaprecord=zeros(length(n34_ind),length(lat),length(lon));
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
            spotz_wind = 0.5*log( (1+spot_wind)./(1-spot_wind) ); % Fishers Z Score
            spotz_precip = 0.5*log( (1+spot_precip)./(1-spot_precip) );
            pr_runcorrz = 0.5*log( (1+pr_runcorr)./(1-pr_runcorr) );
            ts_runcorrz = 0.5*log( (1+ts_runcorr)./(1-ts_runcorr) );
            pr_pc_spotz = prctile(spotz_precip,[2.5,97.5]);
            ts_pc_spotz = prctile(spotz_wind,[2.5,97.5]);
            pr_pc(:,:,j,i) = prctile(spot_precip,[2.5,97.5]);
            ts_pc(:,:,j,i) = prctile(spot_wind,[2.5,97.5]);
            nonstat_precipmap(j,i)=length( find(squeeze(pr_runcorrz(:,j,i))<pr_pc_spotz(1,:)'|...
                                          squeeze(pr_runcorrz(:,j,i))>pr_pc_spotz(2,:)')   );
            nonstat_windmap(j,i)=length( find(squeeze(ts_runcorrz(:,j,i))<ts_pc_spotz(1,:)'|...
                                          squeeze(ts_runcorrz(:,j,i))>ts_pc_spotz(2,:)')   );
            for n=1:(length(ts_runcorrz(:,1,1))+1-run_wdw-runcorr_wdw)
                running_nonstat_precipmap(n,j,i) = length( find(squeeze( ...
                    pr_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))<pr_pc_spotz(1,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))'|...
                    squeeze(pr_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))>pr_pc_spotz(2,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))')   );
                running_nonstat_windmap(n,j,i) = length( find(squeeze( ...
                    ts_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))<ts_pc_spotz(1,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))'|...
                    squeeze(ts_runcorrz(n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw),j,i))>ts_pc_spotz(2,n+runcorr_wdw:(n+run_wdw-1+runcorr_wdw))')   );
            end
            nonstat_windmaprecord(find(squeeze(ts_runcorrz(:,j,i))<ts_pc_spotz(1,:)'|...
                                          squeeze(ts_runcorrz(:,j,i))>ts_pc_spotz(2,:)'),j,i) = 1;
            nonstat_precipmaprecord(find(squeeze(pr_runcorrz(:,j,i))<pr_pc_spotz(1,:)'|...
                                          squeeze(pr_runcorrz(:,j,i))>pr_pc_spotz(2,:)'),j,i) = 1;
                                        
        else
            disp(['Data at ',num2str(lon(i)),'E',num2str(lat(j)),'N does not exist']);
        end
    end
end

save('DataFiles/nonstat_map.mat','nonstat_precipmap','nonstat_windmap', ...
     'pr_pc','ts_pc','running_nonstat_precipmap','running_nonstat_windmap',...
     'nonstat_windmaprecord','nonstat_precipmaprecord');
 
%% Plotting Map
load /srv/ccrc/data34/z3372730/Katana_Data/MATLAB/DataFiles/nonstat_map31yrwdw.mat
NUM_CONTOURS = 10;

subplot(2,1,1);
contourf(lon,lat,nonstat_precipmap,NUM_CONTOURS);
plotworld;
colorbar;
caxis([0, 200]);
colormap(flipud(hot(NUM_CONTOURS)));
title(['Num of nonstationary stations for prec, rcor window:',num2str(window),'yrs'])
subplot(2,1,2);
contourf(lon,lat,nonstat_windmap,NUM_CONTOURS);
plotworld;
caxis([0, 200]);
colorbar
colormap(flipud(hot(NUM_CONTOURS)));
title(['Num of nonstationary stations for temp, rcor window:',num2str(window),'yrs'])

% pcolor(lon,lat,double(nonstat_windmap)); plotworld; colorbar; caxis([0,200]); colormap(flipud(hot));

%% Adding Difference in Percentiles Overlay

subplot(2,1,1)
hold on;
% [c,h]=contour(lon,lat,squeeze(max(abs(squeeze(pr_pc(1,:,:,:)-pr_pc(2,:,:,:))))),[0.1:0.1:0.9],'b');
[c,h]=contour(lon,lat,double(squeeze(abs(mean(pr_pc(1,(window+2):end,:,:)-pr_pc(2,(window+2):end,:,:),2)))),[0.1:0.1:0.9],'b');
clabel(c,h,'Color','b')
hold off;
subplot(2,1,2)
hold on
% [c,h]=contour(lon,lat,squeeze(max(abs(squeeze(ts_pc(1,:,:,:)-ts_pc(2,:,:,:))))),[0.1:0.1:0.9],'b');
[c,h]=contour(lon,lat,double(squeeze(abs(mean(ts_pc(1,(window+2):end,:,:)-ts_pc(2,(window+2):end,:,:),2)))),[0.1:0.1:0.9],'b');
clabel(c,h,'Color','b')
hold off;
set(gcf, 'PaperPosition', [0 0 19 28]); %x_width=19cm y_width=28cm
saveas(gcf,['Plots/nonstatmap+pc_rcor',num2str(window),'yr.jpg'])

% %% Removing insignificant teleconnections for Temperature
% 
% sig_chang         e=ones(size(ts_pc,3),size(ts_pc,4));
% for i=1:size(ts_pc,3)
%     for j=1:size(ts_pc,4)
%         if abs(max(ts_pc(1,:,i,j)-ts_pc(2,:,i,j)))<0.1
%             sig_change(i,j)=nan; % No significant change
%         end
%     end
% end

%% EOFS + Nonstat_wind_map

figure;
pcolor(lon,lat,squeeze(rc31_eof_wind_fm(1,:,:)));
plotworld; hold on;

[c,h]=contour(lon,lat,nonstat_windmap,0:25:100,'m'); hold off;
clabel(c,h,'Color','m')
colormap(b2r(-0.03,0.03));
%% Examining Certain Stations plotting Z-scores

% x_lon = 226; [~,x_ind]= min(abs(lon-x_lon));
% y_lat = 1; [~,y_ind]= min(abs(lat-y_lat));
% 
% load(['Synth_runcorr/',num2str(lon(x_ind)),'E',num2str(lat(y_ind)),'N_syncorr.mat']);
% load('DataFiles/runcorr.mat');
% 
% spotz_wind = 0.5*log( (1+spot_wind)./(1-spot_wind) );
% spotz_precip = 0.5*log( (1+spot_precip)./(1-spot_precip) );
% pr_runcorrz = 0.5*log( (1+pr_runcorr)./(1-pr_runcorr) );
% ts_runcorrz = 0.5*log( (1+ts_runcorr)./(1-ts_runcorr) );
% pr_pc_spotz = prctile(spotz_precip,[2.5,97.5]);
% ts_pc_spotz = prctile(spotz_wind,[2.5,97.5]);
% 
% % Plot of percentiles and GFDL runcorrs
% subplot(2,1,1)
% plot(pr_pc_spotz'); hold on;
% plot(squeeze(pr_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3);
% hold off;
% title(['Percentiles of Z scores of Runing correlations of pr at ', num2str(x_lon),'E, ',num2str(y_lat),'N']);
% subplot(2,1,2)
% plot(ts_pc_spotz'); hold on;
% plot(squeeze(ts_runcorrz(:,y_ind,x_ind)),'k','LineWidth',3);
% hold off;
% title(['Percentiles of Z scores of Runing correlations of ts at ', num2str(x_lon),'E, ',num2str(y_lat),'N']);
% 

% %% Plotting map of running nonstat
% clf;
% for i=1:(length(ts_runcorr(:,1,1))+1-run_wdw-runcorr_wdw)
%     pcolor(lon,lat,squeeze(running_nonstat_windmap(i,:,:))); plotworld; colorbar; caxis([0, 30]);
%     title(['Running nonstationary years in 30 year windows (temperature)',num2str(i)])
%     M(i)=getframe(gcf);
% end
% movie2avi(M,'Plots/running_nonstat_wind.avi','FPS',8)
% 
% %% Plotting map of nonstat time series
% clf;
% for i=1:(length(ts_runcorr(:,1,1)))
%     pcolor(lon,lat,squeeze(nonstat_windmaprecord(i,:,:))); plotworld; colorbar; caxis([0 1]);
%     colormap(flipud(gray)); title(['Nonstationary Areas (temperature)',num2str(i)])
%     M(i)=getframe(gcf);
% end
% movie2avi(M,'Plots/nonstat_wind.avi','FPS',8)

%% Plotting Different groups of nonstationarities

load DataFiles/nonstat_map.mat
load DataFiles/runcorr.mat

corr_precip = zeros(size(apr,2),size(apr,3));
for i=1:size(apr,2)
    for j=1:size(apr,3)
        corr_precip(i,j) = corr(n34_ind,apr(:,i,j));
    end
end

corr_wind = zeros(size(ats,2),size(ats,3));
for i=1:size(ats,2)
    for j=1:size(ats,3)
        corr_wind(i,j) = corr(n34_ind,ats(:,i,j));
    end
end
clear i j

run_corr_low = zeros(size(corr_wind));
mean_corr_low = zeros(size(corr_wind));
for i=1:size(apr,2)
    for j=1:size(apr,3)
        if (~isempty(find(abs(ts_runcorr(:,i,j)) < 0.3)) & ...
                                 nonstat_windmap(i,j) > 50 & ...
                                 abs(corr_wind(i,j)) > 0.3       )
            run_corr_low(i,j) = 1;
        elseif (~isempty(find(abs(ts_runcorr(:,i,j)) > 0.3)) & ...
                                 nonstat_windmap(i,j) > 50 & ...
                                 abs(corr_wind(i,j)) < 0.3       )
             mean_corr_low(i,j) = 1;
        end
    end
end

[is js] = ind2sub(size(corr_wind),find(run_corr_low));
[is2 js2] = ind2sub(size(corr_wind),find(mean_corr_low));
plotworld; hold on; run_hnd=scatter(lon(js), lat(is),'r.'); mean_hnd=scatter(lon(js2),lat(is2),'b.');
xlim([0,360]); ylim([-90,90]); hold off;
legend([run_hnd, mean_hnd],'abs(runcorr) > 0.3','corr\_wind > 0.3')
