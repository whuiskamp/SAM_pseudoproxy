% Extract correlations for all the calibration windows and look at the std.
% deviation of correlations.
% Updated May 7 2020

clear
precip_corrwindows = nan(3,10,90,144); % 3 windowsizes, 10 calibration windows, lat, lon
sat_corrwindows = nan(3,10,90,144);
NUM_YRS = 500;
NUM_CAL_WDW = 10;

for windowsize = [31 61 91]
  
  if windowsize > 31
      clear CAL_WDW
  end
  
  if windowsize == 31
      i = 1;
  elseif windowsize == 61
      i = 2;
  elseif windowsize == 91
      i = 3;
  end
  
  for c=0:9
       overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
       CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
  end
  
  
  for c = 1:10
    load(['Proxies/NoResample/',num2str(windowsize),'yrWindow/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/10stns_1000prox.mat']);
    precip_corrwindows(i,c,:,:) = corr_precip(:,:);
    sat_corrwindows(i,c,:,:) = corr_sat(:,:);
  end
end

std_precip = nan(3,90,144);
std_sat = nan(3,90,144);


for w = 1:3 % window
    for i = 1:90
        for j = 1:144
            std_precip(w,i,j) = squeeze(nanstd(precip_corrwindows(w,:,i,j)));
            std_sat(w,i,j) = squeeze(nanstd(sat_corrwindows(w,:,i,j)));
        end
    end
end

figure(1)
precip = tight_subplot(1,3,[0.05 0.01],[0.10 0.01],[0.1 0.01]);
figure(2)
sat = tight_subplot(1,3,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

for i = [precip sat]
    for j = 1:3
        axes(i(j));
        axesm('MapProjection','stereo','origin',[-90,0],'MapLatLimit',[-90 -20])
        framem
        load coast
        plotm(lat,long,'k')
        gridm %If you want gridlines
    end
end

load(['DataFiles/model_output.mat']);

figure(1)
for j = 1:3
    axes(precip(j));
    contourfm(lat,lon,squeeze(std_precip(j,:,:)),'linestyle','none')
end
        
figure(2)
for j = 1:3
    axes(sat(j));
    contourfm(lat,lon,squeeze(std_sat(j,:,:)),'linestyle','none')
end

std_precip = permute(std_precip,[2 3 1]);
std_sat = permute(std_sat,[2 3 1]);

nccreate('corr_std.nc','std_precip_wdw','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144,'time',3},'Format','classic')
ncwrite('corr_std.nc','std_precip_wdw',std_precip(1:45,1:144,1:3))
nccreate('corr_std.nc','std_sat_wdw','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144,'time',3},'Format','classic')
ncwrite('corr_std.nc','std_sat_wdw',std_sat(1:45,1:144,1:3))

%% How many times does each point on the map meet the correlation criteria?

for windowsize = 1:3
    for i = 1:90
        for j = 1:144
            corr_precip_count(windowsize,i,j) = sum(~isnan(precip_corrwindows(windowsize,:,i,j)));
            corr_sat_count(windowsize,i,j) = sum(~isnan(sat_corrwindows(windowsize,:,i,j)));
        end
    end
end


%% Same thing but with the full fields of the 500yr running correlation

clear
precip_fullcorr = nan(3,500,90,144); % 3 windowsizes, 500yrs, lat, lon
sat_fullcorr = nan(3,500,90,144);
NUM_YRS = 500;

load(['DataFiles/runcorr31yrwdw.mat'])
precip_fullcorr(1,:,:,:) = precip_runcorr;
sat_fullcorr(1,:,:,:) = sat_runcorr;
load(['DataFiles/runcorr61yrwdw.mat'])
precip_fullcorr(2,:,:,:) = precip_runcorr;
sat_fullcorr(2,:,:,:) = sat_runcorr;
load(['DataFiles/runcorr91yrwdw.mat'])
precip_fullcorr(3,:,:,:) = precip_runcorr;
sat_fullcorr(3,:,:,:) = sat_runcorr;

for w = 1:3 % window
    for i = 1:90
        for j = 1:144
            std_precip(w,i,j) = squeeze(nanstd(precip_fullcorr(w,:,i,j)));
            std_sat(w,i,j) = squeeze(nanstd(sat_fullcorr(w,:,i,j)));
        end  
    end
    w
end

figure(1)
precip = tight_subplot(1,3,[0.05 0.01],[0.10 0.01],[0.1 0.01]);
figure(2)
sat = tight_subplot(1,3,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

for i = [precip sat]
    for j = 1:3
        axes(i(j));
        set(gca, 'CLim', [0, 0.25])
        axesm('MapProjection','stereo','origin',[-90,0],'MapLatLimit',[-90 -20])
        framem
        load coast
        plotm(lat,long,'k')
        gridm %If you want gridlines
    end
end

load(['DataFiles/model_output.mat']);

lon(1) = 0.5; lon(144) = 360;
levels = 0:0.005:0.25;

figure(1)
for j = 1:3
    axes(precip(j));
    contourfm(lat,lon,squeeze(std_precip(j,:,:)),'linestyle','none','levellist',levels)
end
        
figure(2)
for j = 1:3
    axes(sat(j));
    contourfm(lat,lon,squeeze(std_sat(j,:,:)),'linestyle','none','levellist',levels)
end

std_precip = permute(std_precip,[2 3 1]);
std_sat = permute(std_sat,[2 3 1]);

load(['DataFiles/model_output.mat']);

nccreate('corr_std.nc','std_precip_full','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144,'time',3},'Format','classic')
ncwrite('corr_std.nc','std_precip_full',std_precip(1:45,1:144,1:3))
nccreate('corr_std.nc','std_sat_full','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144,'time',3},'Format','classic')
ncwrite('corr_std.nc','std_sat_full',std_sat(1:45,1:144,1:3))

nccreate('corr_std.nc','latitude','Dimensions',...
    {'latitude',45},'Format','classic')
nccreate('corr_std.nc','longitude','Dimensions',...
    {'longitude',144},'Format','classic')
nccreate('corr_std.nc','time','Dimensions',...
    {'time',3},'Format','classic')

ncwrite('corr_std.nc','latitude',lat(1:45,1))
ncwrite('corr_std.nc','longitude',lon(1:144,1))
ncwriteatt('corr_std.nc','latitude','axis','Y')
ncwriteatt('corr_std.nc','longitude','axis','X')
ncwriteatt('corr_std.nc','time','axis','Z')

ncwriteatt('corr_std.nc','latitude','axis','Y')
ncwriteatt('corr_std.nc','latitude','standard_name','latitude')
ncwriteatt('corr_std.nc','latitude','units','degrees_north')

ncwriteatt('corr_std.nc','longitude','axis','Y')
ncwriteatt('corr_std.nc','longitude','standard_name','longitude')
ncwriteatt('corr_std.nc','longitude','units','degrees_east')




