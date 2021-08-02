% Calculate the range in size of the proxy pool that you are sampling from

clear

% load land mask and setup new variables
land = ncread('sftlf_A1.static.nc','sftlf');
load('DataFiles/model_output','precip_detr','sat_detr','SAM')
load('DataFiles/ENSO_corrs.mat','SAT_corr_sig','pre_corr_sig')

% Next we need to create 10 evenly spaced calibration windows in the data,
% of the length of our correlation windows.

NUM_YRS = 500; NUM_TRIALS = 1;
NUM_CAL_WDW = 10; 
STN_MAX = 1;
siterange_precip = NaN(3,10,6); % When running this part, set NUM_TRIALS to 1 to save time
siterange_SAT = NaN(3,10,6);
precip_land = nan(size(precip_detr));
sat_land = nan(size(sat_detr));
irange = 144; jrange = 45;

for windowsize = [31 61 91]
  
  if windowsize > 31
      clear CAL_WDW
  end
  
  for c = 0:9
       overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
       CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
  end
     
% Load topography and eliminate all data that isn't land. This will save time with the correlations.

  for i = 1:irange
    for j = 1:jrange
        if land(i,j) > 0 && isnan(pre_corr_sig(floor(windowsize/30),j,i))
            precip_land(:,j,i) = precip_detr(:,j,i);
        else precip_land(:,j,i) = NaN;
        end
    end
  end
   
  
  for i = 1:irange
    for j = 1:jrange
        if land(i,j) > 0 && isnan(SAT_corr_sig(floor(windowsize/30),j,i))
            sat_land(:,j,i) = sat_detr(:,j,i);
        else sat_land(:,j,i) = NaN;
        end
    end
  end
   
  
% Calculate correlations 

  for c = 1:10
    for i = 1:45
        for j = 1:144
            corr_precip(i,j) = corr(SAM(CAL_WDW(c,:)), precip_land(CAL_WDW(c,:),i,j));
            corr_sat(i,j) = corr(SAM(CAL_WDW(c,:)), sat_land(CAL_WDW(c,:),i,j));
        end
    end
    %mkdir(['Proxies/',num2str(windowsize),'yrWindow/No_Resample/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end))]);


    % Mask correlations so we only have land values and correlations greater
    % than or equal to +-0.3

    for i = 1:irange
        for j = 1:jrange
            if abs(corr_precip(j,i)) < 0.3
                corr_precip(j,i) = NaN;
            end
        end
    end
    
    for i = 1:irange
        for j = 1:jrange
            if abs(corr_sat(j,i)) < 0.3
                corr_sat(j,i) = NaN;
            end
        end
    end
        for region = 1:6
    
            Gl_lon = 1:144; Gl_lat = 1:45;
            SA_lon = 110:130; SA_lat = 15:45;
            Au_lon = 40:75; Au_lat = 15:41;
			AA_lon = 1:144; AA_lat = 15:45;
			AAo_lon = 1:144; AAo_lat = 1:15;
			SoA_lon = 1:25; SoA_lat = 25:45;
        
            if region == 1
                lat = Gl_lat; lon = Gl_lon;
                [I,J] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                siterange_precip(floor(windowsize/30),c,region) = length(I);
                [I_sat,J_sat] = find(~isnan(corr_sat(lat,lon)));
                siterange_SAT(floor(windowsize/30),c,region) = length(I_sat);
            elseif region == 2
                lat = SA_lat; lon = SA_lon;
                [I_SA,J_SA] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                siterange_precip(floor(windowsize/30),c,region) = length(I_SA);
                [I_sat_SA,J_sat_SA] = find(~isnan(corr_sat(lat,lon)));
				siterange_SAT(floor(windowsize/30),c,region) = length(I_sat_SA);
			elseif region == 3
                lat = Au_lat; lon = Au_lon;
                [I_Au,J_Au] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                siterange_precip(floor(windowsize/30),c,region) = length(I_Au);
                [I_sat_Au,J_sat_Au] = find(~isnan(corr_sat(lat,lon)));
				siterange_SAT(floor(windowsize/30),c,region) = length(I_sat_Au);
			elseif region == 4
                lat = AA_lat; lon = AA_lon;
                [I_AA,J_AA] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                siterange_precip(floor(windowsize/30),c,region) = length(I_AA);
                [I_sat_AA,J_sat_AA] = find(~isnan(corr_sat(lat,lon)));
                siterange_SAT(floor(windowsize/30),c,region) = length(I_sat_AA);
            elseif region == 5
                lat = AAo_lat; lon = AAo_lon;
                [I_AAo,J_AAo] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                siterange_precip(floor(windowsize/30),c,region) = length(I_AAo);
				[I_sat_AAo,J_sat_AAo] = find(~isnan(corr_sat(lat,lon)));
				siterange_SAT(floor(windowsize/30),c,region) = length(I_sat_AAo);
			elseif region == 6
                lat = SoA_lat; lon = SoA_lon;
                [I_SoA,J_SoA] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                siterange_precip(floor(windowsize/30),c,region) = length(I_SoA);
        		[I_sat_SoA,J_sat_SoA] = find(~isnan(corr_sat(lat,lon)));
				siterange_SAT(floor(windowsize/30),c,region) = length(I_sat_SoA);				
            end
        end
     c
  end
windowsize
end

for i = 1:6
    for j = 1:3
        precip_min(j,i) = min(siterange_precip(j,:,i));
        sat_min(j,i) = min(siterange_SAT(j,:,i));
        precip_max(j,i) = max(siterange_precip(j,:,i));
        sat_max(j,i) = max(siterange_SAT(j,:,i));
    end
end

save('site_range_ENSO.mat','precip_min','sat_min','precip_max','sat_max')


% sites(1,:,:) = stn_lat_SoA(:,:);
% sites(2,:,:) = stn_lon_SoA(:,:);
% 
% sites_test = squeeze(sites(:,1,:))';
% sites_2 = reshape(sites,[2 70000])';
% 
% u = unique(sites_test,'rows');
% n=histc(test2,u)
