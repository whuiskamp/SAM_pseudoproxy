%% Figure out which regions would provide best proxies
% Willem Huiskamp, Sept. 2017
%
% This script seeks to highlight regions with the strongest relationship
% between SAM and temp/ SAT fields.
% 
% Setup
load('DataFiles/model_output.mat','lat');
NUM_EOFs = 10; % Number of EOFs we want to retain later
lat_2 = 31; % change this from 45 to get some sub-set of the SH
lats = zeros(144*lat_2,1);
tick = 1:lat_2:(144*lat_2);
for i = tick
    lats(i:i+(lat_2-1),1) = sqrt(cos(lat(1:lat_2)*3.141/180)); % we are weighting by the sqrt cos of latitude (in radians, of course)
end
weight = diag(lats);

% 1: Import running correlations between SAM and detrended SAT and precip
% anomalies
% 

for windowsize = [31 61 91]
    load(['DataFiles/runcorr',num2str(windowsize),'yrwdw.mat']);
    % SAT - remove NaNs and reshape
    sat_runcorr = sat_runcorr(ceil(windowsize/2):500-floor(windowsize/2),1:lat_2,:); % we only need the SH and eliminate the NaNs from the running correlation at either end.
    sat_runcorr_2 = reshape(sat_runcorr,size(sat_runcorr,1),size(sat_runcorr,2)*size(sat_runcorr,3)); % reshape so that we have a 2D matrix of time X cells
    precip_runcorr = precip_runcorr(ceil(windowsize/2):500-floor(windowsize/2),1:lat_2,:);
    precip_runcorr_2 = reshape(precip_runcorr,size(precip_runcorr,1),size(precip_runcorr,2)*size(precip_runcorr,3));
    
    % Create anomalies
    [n,p] = size(sat_runcorr_2); [q,w] = size(precip_runcorr_2);
    Xbar = mean(sat_runcorr_2,1); Xbar_p = mean(sat_runcorr_2,1);
    sat_anom = sat_runcorr_2-ones(n,1)*Xbar; precip_anom = precip_runcorr_2-ones(q,1)*Xbar_p;
    % Weight by latitude
    sat_weighted = weight*sat_anom'; precip_weighted = weight*precip_anom';
    
    [sat_EOFs,sat_PC,~,~,sat_var_exp] = pca(sat_weighted'); % Calculate EOFs
    [precip_EOFs,precip_PC,~,~,precip_var_exp] = pca(precip_weighted');
    rotEOFs_sat = rotatefactors(sat_EOFs(:,1:NUM_EOFs),'Method','varimax','maxit',1000); % Rotate first five EOFs to maximise variance. Need to increase the maxit to 1000 from 250, otherwise it breaks. Less than 10 is best when using Varimax rotation
    rotEOFs_precip = rotatefactors(precip_EOFs(:,1:NUM_EOFs),'Method','varimax','maxit',1000);
    
    sat_rot_EOFs = reshape(rotEOFs_sat',NUM_EOFs,lat_2,144); % Reshape back into 3D matrix (PCs, lat, lon)
    precip_rot_EOFs = reshape(rotEOFs_precip',NUM_EOFs,lat_2,144);
    sat_EOFs = reshape(sat_EOFs',size(sat_EOFs,2),lat_2,144);
    precip_EOFs = reshape(precip_EOFs',size(precip_EOFs,2),lat_2,144);
    
    save(['DataFiles/runcorr',num2str(windowsize),'yrwdw_30S.mat'],'sat_EOFs','sat_rot_EOFs','sat_PC','sat_var_exp','precip_EOFs','precip_rot_EOFs','precip_PC','precip_var_exp')%,'-append')
    clear sat_runcorr precip_runcorr sat_anom precip_anom sat_weighted precip_weighted sat_EOFs sat_PC sat_var_exp precip_EOFs precip_PC precip_var_exp rotEOFs_sat rotEOFs_precip sat_rot_EOFs precip_rot_EOFs
end

% Normalise the PCs and plot them (first 3)
windowsize = 31;

for i = 1:(size(sat_PC,2))
    zsat_PC(:,i) = zscore(sat_PC(:,i));
    zprecip_PC(:,i) = zscore(precip_PC(:,i));
    
end

for w = [31 61 91]
    load(['DataFiles/runcorr',num2str(w),'yrwdw_30S.mat'],'sat_PC','precip_PC','sat_EOFs','precip_EOFs')
    for i = 1:(size(sat_PC,2))
        zsat_PC(:,i) = zscore(sat_PC(:,i));
        zprecip_PC(:,i) = zscore(precip_PC(:,i));
        std_sat_PC(floor(w/30),1,i) = std(sat_PC(:,i));
        std_precip_PC(floor(w/30),1,i) = std(precip_PC(:,i));
    end
    
    figure(floor(w/30))
    for i = 1:3
        subplot(3,1,i)
        plot(zsat_PC(:,i),'color','black','linewidth',2)
        title(['PC',num2str(i),' SAT - Black, Precip - Red'])
        hold
        plot(zprecip_PC(:,i),'color','red','linewidth',2)
    end
        h=gcf;
        set(h,'PaperPositionMode','auto');         
        set(h,'PaperOrientation','landscape');
        set(h,'PaperUnits','normalized');
        set(h,'PaperPosition', [0 0 1 1]);
    %print(gcf,'-dpdf','-painters',['3PCs_runncorr_',num2str(w),'yrWdW.pdf'])
    clear zsat_PC zprecip_PC
end
save(['DataFiles/runcorr31yrwdw_30S.mat'],'std_sat_PC','std_precip_PC','-append')
% Regress first N PCs onto SAT & Precip field
N = 3; % select number of PCs to ivestigate


for windowsize = [31 61 91]
    load(['DataFiles/model_output.mat'],'sat_detr','precip_detr');
    load(['DataFiles/runcorr',num2str(windowsize),'yrwdw_30S.mat'],'sat_PC','precip_PC')
    sat_detr = sat_detr(ceil(windowsize/2):500-floor(windowsize/2),:,:);
    precip_detr = precip_detr(ceil(windowsize/2):500-floor(windowsize/2),:,:);
    for n = 1:N
        for i = 1:90
            for j = 1:144
                SAM_SAT_regr(n,i,j) = regress(sat_PC(:,n),squeeze(sat_detr(:,i,j)));
                SAM_precip_regr(n,i,j) = regress(precip_PC(:,n),squeeze(precip_detr(:,i,j)));
            end
        end
    end
    save(['DataFiles/runcorr',num2str(windowsize),'yrwdw_30S.mat'],'SAM_SAT_regr','SAM_precip_regr','-append')
    clear sat_detr precip_detr SAM_SAT_regr SAM_precip_regr sat_PC precip_PC
end

% Correlate global average surface air temp with the first PC

load('DataFiles/model_output.mat','sat_detr')
test = squeeze(mean(sat_detr,3));
sat_glob_mean = squeeze(mean(test,2));

for wdw = [31 61 91]
    for i = 1:3
        load(['DataFiles/runcorr',num2str(wdw),'yrwdw_30S.mat'],'sat_PC')
        PC_SAT_CORR(floor(wdw/30),i) = corr(sat_glob_mean(1:size(sat_PC,1),:),sat_PC(:,i));
    end
end
        

%% Plotting 
% Setup
clear
windowsize = 91; % change this for other plots.
lat_2 = 31;
load(['DataFiles/runcorr',num2str(windowsize),'yrwdw_30S.mat'],'sat_EOFs','sat_rot_EOFs','sat_PC','sat_var_exp','precip_EOFs','precip_rot_EOFs','precip_PC','precip_var_exp');
levels = -0.1:0.01:0.1;
levels_p = 0.005:0.02:0.1;
levels_n = -0.1:0.02:0;


% This plots the first 3 EOFs of the running correlation between SAT and
% SAM
load('DataFiles/runcorr31yrwdw_30S.mat','std_sat_PC','std_precip_PC')
subplot1(3, 1, 'Gap', [.01 .03], 'XTickL', 'Margin', 'YTickL', 'Margin');
for i = 1:3;
    subplot1(i)
    axesm('MapProjection','eqdcylin','origin',[0,180],'MapLatLimit',[-90 -30],'meridianlabel','on','parallellabel','on','PLabelMeridian',0,'Mlabelparallel','south')
    framem
    gridm %If you want gridlinesload coast
    load('DataFiles/model_output.mat','lat','lon')
    contourfm(lat(1:lat_2),lon,squeeze(sat_EOFs(i,:,:)*(squeeze(std_sat_PC(floor(windowsize/30),:,i)))),levels,'linestyle','none') 
    contourm(lat(1:lat_2),lon,squeeze(sat_EOFs(i,:,:)*(squeeze(std_sat_PC(floor(windowsize/30),:,i)))),levels_p,'k');contourm(lat(1:lat_2),lon,squeeze(sat_EOFs(i,:,:)),levels_n,'k','linestyle','--');
    colormap(b2r(-0.1,0.1))
    load coast
    plotm(lat,long,'k','linewidth',2)
    title(['SAT EOF',num2str(i),' for',num2str(windowsize),'yr Window - ',num2str(round(sat_var_exp(i,:))),'% var. exp.'])
    colorbar
    clear lat long
end
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf,'-dpdf','-painters',['SAT_EOFs_',num2str(windowsize),'yr_30S'])


figure(2)
subplot1(3, 1, 'Gap', [.01 .03], 'XTickL', 'Margin', 'YTickL', 'Margin');
for i = 1:3; % Same plots for precip
    subplot1(i)
    axesm('MapProjection','eqdcylin','origin',[0,180],'MapLatLimit',[-90 -30],'meridianlabel','on','parallellabel','on','PLabelMeridian',0,'Mlabelparallel','south')
    framem
    gridm %If you want gridlinesload coast
    load('DataFiles/model_output.mat','lat','lon')
    contourfm(lat(1:lat_2),lon,squeeze(precip_EOFs(i,:,:)*(squeeze(std_precip_PC(floor(windowsize/30),:,i)))),levels,'linestyle','none')
    contourm(lat(1:lat_2),lon,squeeze(precip_EOFs(i,:,:)*(squeeze(std_precip_PC(floor(windowsize/30),:,i)))),levels_p,'k');contourm(lat(1:lat_2),lon,squeeze(precip_EOFs(i,:,:)*(squeeze(std_precip_PC(floor(windowsize/30),:,i)))),levels_n,'k','linestyle','--');
    colormap(b2r(-0.1,0.1))
    load coast
    plotm(lat,long,'k','linewidth',2)
    title(['Precip EOF',num2str(i),' for ',num2str(windowsize),'yr Window - ',num2str(round(precip_var_exp(i,:))),'% var. exp.'])
    colorbar
    clear lat long
end

h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf,'-dpdf','-painters',['Precip_EOFs_',num2str(windowsize),'yr_30S'])


% Now do the same for rotated EOFs
levels = -0.075:0.005:0.075;
levels_p = 0.005:0.005:0.075;
levels_n = -0.075:0.005:0;

for i = 1:3;
    subplot(3,1,i)
    axesm('MapProjection','eqdcylin','origin',[0,180],'MapLatLimit',[-90 -0],'meridianlabel','on','parallellabel','on','PLabelMeridian',0,'Mlabelparallel','south')
    framem
    gridm %If you want gridlinesload coast
    load('DataFiles/model_output.mat','lat','lon')
    contourfm(lat(1:lat_2),lon,squeeze(sat_rot_EOFs(i,:,:)),levels,'linestyle','none')
    contourm(lat(1:lat_2),lon,squeeze(sat_rot_EOFs(i,:,:)),levels_p,'k');contourm(lat(1:lat_2),lon,squeeze(sat_rot_EOFs(i,:,:)),levels_n,'k','linestyle','--');
    colormap(b2r(-0.04,0.04))
    load coast
    plotm(lat,long,'k','linewidth',2)
    title(['Rot. SAT EOF',num2str(i),' for',num2str(windowsize),'yr Window - ',num2str(round(sat_var_exp(i,:))),'% var. exp.'])
    colorbar
    clear lat long
end

print(['SAT_rot_EOFs_',num2str(windowsize),'yr'],'-dpng','-r500')

figure(2)
for i = 1:3; % Same plots for precip
    subplot(3,1,i)
    axesm('MapProjection','eqdcylin','origin',[0,180],'MapLatLimit',[-90 -0],'meridianlabel','on','parallellabel','on','PLabelMeridian',0,'Mlabelparallel','south')
    framem
    gridm %If you want gridlinesload coast
    load('DataFiles/model_output.mat','lat','lon')
    contourfm(lat(1:lat_2),lon,squeeze(precip_rot_EOFs(i,:,:)),levels,'linestyle','none')
    contourm(lat(1:lat_2),lon,squeeze(precip_rot_EOFs(i,:,:)),levels_p,'k');contourm(lat(1:lat_2),lon,squeeze(precip_rot_EOFs(i,:,:)),levels_n,'k','linestyle','--');
    colormap(b2r(-0.04,0.04))
    load coast
    plotm(lat,long,'k','linewidth',2)
    title(['Rot. Precip EOF',num2str(i),' for ',num2str(windowsize),'yr Window - ',num2str(round(precip_var_exp(i,:))),'% var. exp.'])
    colorbar
    clear lat long
end

print(['precip_rot_EOFs_',num2str(windowsize),'yr'],'-dpng','-r500')

figure(3)
levels_p = 0.1:0.2:2;
levels_n = -0.1:-0.2:-2;
windowsize = 61;
load(['DataFiles/runcorr',num2str(windowsize),'yrwdw_30S.mat'],'SAM_SAT_regr','SAM_precip_regr')

for i = 1:3;
    subplot(3,1,i)
    load('DataFiles/model_output.mat','lat','lon')
    contourf(lon,lat,squeeze(SAM_SAT_regr(i,:,:)),'linestyle','none'); hold on;
    contour(lon,lat,squeeze(SAM_SAT_regr(i,:,:)),levels_p,'k');contour(lon,lat,squeeze(SAM_SAT_regr(i,:,:)),levels_n,'k','linestyle','--');
    colormap(b2r(-2,2))
    plotworld
    title(['SAM-SAT regres. EOF',num2str(i),' for',num2str(windowsize),'yr Window'])
    colorbar
end

print(['SAM_SAT_regr_',num2str(windowsize),'yr'],'-dpng','-r500')

% Plot of explained variance for each Variable/ windowsize

for i = [31 61 91]
    load(['DataFiles/runcorr',num2str(i),'yrwdw.mat'],'sat_var_exp','precip_var_exp');
    var_exp_sat(1:size(sat_var_exp,1),floor(i/30)) = sat_var_exp(:,1);
    var_exp_precip(1:size(precip_var_exp,1),floor(i/30)) = precip_var_exp(:,1);
    clear sat_var_exp precip_var_exp
end

figure
subplot(2,1,1)














