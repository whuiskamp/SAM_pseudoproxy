% Script to generate pseudoproxy groups that exclude sites with a 
% significant correlation with ENSO
% Willem Huiskamp, 2015
%
% Method roughly follows that of Batehup et al. 2015.
% Proxies are only selected from land tiles and must have an absolute
% correlation with SAM greater than 0.3 over a calibration window
% (or over the entire time series). We then generate a network of proxies
% between one and 70 sites.

% This version ensures no site is selected twice for the same network

clear

% load land mask and setup new variables
land = ncread('sftlf_A1.static.nc','sftlf');
load('DataFiles/model_output','precip_detr','sat_detr','SAM')
load('site_range_ENSO','precip_min','sat_min')
load('DataFiles/ENSO_corrs.mat','SAT_corr_sig','pre_corr_sig')
land = land';
% Next we need to create 10 evenly spaced calibration windows in the data,
% of the length of our correlation windows.

NUM_YRS = 500; NUM_TRIALS = 1000;
NUM_CAL_WDW = 10; 
STN_MAX = 70;

for windowsize = [31 61 91]
  
  if windowsize > 31
      clear CAL_WDW
  end
  
  for c = 0:9
       overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
       CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
  end
     
% Load topography and eliminate all data that isn't land. This will save time with the correlations.

  
    for i = 1:45
      for j = 1:144
        if land(i,j) > 0 
               precip_land(:,i,j) = precip_detr(:,i,j);
          else precip_land(:,i,j) = NaN;
        end
      end
    end
  
  
  
    for i = 1:45
      for j = 1:144
        if land(i,j) > 0 
               sat_land(:,i,j) = sat_detr(:,i,j);
          else sat_land(:,i,j) = NaN;
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
    mkdir(['Proxies/NoENSO/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end))]);


    % Mask correlations so we only have land values and correlations greater
    % than or equal to +-0.3

    for i = 1:45
        for j = 1:144
            if abs(corr_precip(i,j)) < 0.3 && isnan(pre_corr_sig(floor(windowsize/30),i,j))
                corr_precip(i,j) = NaN;
            end
        end
    end
    
    for i = 1:45
        for j = 1:144
            if abs(corr_sat(i,j)) < 0.3 && isnan(SAT_corr_sig(floor(windowsize/30),i,j))
                corr_sat(i,j) = NaN;
            end
        end
    end
    % Time to generate our random proxy group - Each region has its own max
    % #STNS determined by the siterange script. We use the minimum # found
    % in any one of the cal. windows (up to 70, of course).
    
    % Make Global (SH)=1, SA=2, Au/NZ=3, SH without AA = 4, Antarctica = 5, Africa = 6 
        for region = 1:6
            Gl_lon = 1:144; Gl_lat = 1:45;
            SA_lon = 110:130; SA_lat = 15:45;
            Au_lon = 40:75; Au_lat = 15:41;
			AA_lon = 1:144; AA_lat = 15:45;
			AAo_lon = 1:144; AAo_lat = 1:15;
			SoA_lon = 1:25; SoA_lat = 25:45;
        
            if region == 1
                for NUM_STNS = 1:STN_MAX
                    stn_lat = zeros(NUM_TRIALS,NUM_STNS); sat_lat = zeros(NUM_TRIALS,NUM_STNS);
                    stn_lon = zeros(NUM_TRIALS,NUM_STNS); sat_lon = zeros(NUM_TRIALS,NUM_STNS);
                    lat = Gl_lat; lon = Gl_lon;
        
                    [I,J] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                                                    
                    for m=1:NUM_TRIALS
                        [stn_lat(m,:),idx] = datasample(I,NUM_STNS,'Replace',false);
                        stn_lon(m,:) = J(idx);
                    end
				
                    [I_sat,J_sat] = find(~isnan(corr_sat(lat,lon)));
				
                    for m=1:NUM_TRIALS
                        [sat_lat(m,:),idx] = datasample(I_sat,NUM_STNS,'Replace',false);
                        sat_lon(m,:) = J_sat(idx);
                    end
				
                    save(['Proxies/NoENSO/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                        'stn_lat','stn_lon','I','J','corr_precip','sat_lat','sat_lon','I_sat','J_sat','corr_sat','windowsize');
                end
            elseif region == 2
                lat = SA_lat; lon = SA_lon;
                [I_SA,J_SA] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                [I_sat_SA,J_sat_SA] = find(~isnan(corr_sat(lat,lon)));
                for NUM_STNS = 1:precip_min(floor(windowsize/30),2)
                    stn_lat_SA = zeros(NUM_TRIALS,NUM_STNS);
                    stn_lon_SA = zeros(NUM_TRIALS,NUM_STNS);
                    for m=1:NUM_TRIALS
                        [stn_lat_SA(m,:),idx] = datasample(I_SA,NUM_STNS,'Replace',false);
                        stn_lon_SA(m,:) = J_SA(idx);
                    end
                    save(['Proxies/NoENSO/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                       'stn_lat_SA','stn_lon_SA','I_SA','J_SA','-append');
                end
                for NUM_STNS = 1:sat_min(floor(windowsize/30),2)
                    sat_lat_SA = zeros(NUM_TRIALS,NUM_STNS);
                    sat_lon_SA = zeros(NUM_TRIALS,NUM_STNS);
                    for m=1:NUM_TRIALS
                        [sat_lat_SA(m,:),idx] = datasample(I_sat_SA,NUM_STNS,'Replace',false);
                        sat_lon_SA(m,:) = J_sat_SA(idx);
                    end
                    save(['Proxies/NoENSO/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                       'sat_lat_SA','sat_lon_SA','I_sat_SA','J_sat_SA','-append');
                end                                
			elseif region == 3
                lat = Au_lat; lon = Au_lon;
                [I_Au,J_Au] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon               
                [I_sat_Au,J_sat_Au] = find(~isnan(corr_sat(lat,lon)));
                for NUM_STNS = 1:precip_min(floor(windowsize/30),3)
                    stn_lat_AuNz = zeros(NUM_TRIALS,NUM_STNS);stn_lon_AuNz = zeros(NUM_TRIALS,NUM_STNS);
                    for m=1:NUM_TRIALS
                        [stn_lat_AuNz(m,:),idx] = datasample(I_Au,NUM_STNS,'Replace',false);
                        stn_lon_AuNz(m,:) = J_Au(idx);
                    end
                    save(['Proxies/NoENSO/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                        'stn_lat_AuNz','stn_lon_AuNz','I_Au','J_Au','-append');
                end
                
                for NUM_STNS = 1:sat_min(floor(windowsize/30),3)
                    sat_lat_AuNz = zeros(NUM_TRIALS,NUM_STNS);sat_lon_AuNz = zeros(NUM_TRIALS,NUM_STNS);
                    for m=1:NUM_TRIALS
                        [sat_lat_AuNz(m,:),idx] = datasample(I_sat_Au,NUM_STNS,'Replace',false);
                        sat_lon_AuNz(m,:) = J_sat_Au(idx);
                    end
                    save(['Proxies/NoENSO/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                        'sat_lat_AuNz','sat_lon_AuNz','I_sat_Au','J_sat_Au','-append');
                end
			elseif region == 4
                lat = AA_lat; lon = AA_lon;
                [I_AA,J_AA] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                [I_sat_AA,J_sat_AA] = find(~isnan(corr_sat(lat,lon)));
                for NUM_STNS = 1:STN_MAX
                    stn_lat_noAA = zeros(NUM_TRIALS,NUM_STNS); sat_lat_noAA = zeros(NUM_TRIALS,NUM_STNS);
                    stn_lon_noAA = zeros(NUM_TRIALS,NUM_STNS); sat_lon_noAA = zeros(NUM_TRIALS,NUM_STNS);
                    for m=1:NUM_TRIALS
                        [stn_lat_noAA(m,:),idx] = datasample(I_AA,NUM_STNS,'Replace',false);
                        stn_lon_noAA(m,:) = J_AA(idx);
                    end
                    for m=1:NUM_TRIALS
                        [sat_lat_noAA(m,:),idx] = datasample(I_sat_AA,NUM_STNS,'Replace',false);
                        sat_lon_noAA(m,:) = J_sat_AA(idx);
                    end
	
                    save(['Proxies/NoENSO/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                        'stn_lat_noAA','stn_lon_noAA','I_AA','J_AA','sat_lat_noAA','sat_lon_noAA','I_sat_AA','J_sat_AA','-append');
                end
            elseif region == 5
                lat = AAo_lat; lon = AAo_lon;
                [I_AAo,J_AAo] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                [I_sat_AAo,J_sat_AAo] = find(~isnan(corr_sat(lat,lon)));
                for NUM_STNS = 1:STN_MAX
                    stn_lat_AAo = zeros(NUM_TRIALS,NUM_STNS);stn_lon_AAo = zeros(NUM_TRIALS,NUM_STNS); 
                    for m=1:NUM_TRIALS
                        [stn_lat_AAo(m,:),idx] = datasample(I_AAo,NUM_STNS,'Replace',false);
                        stn_lon_AAo(m,:) = J_AAo(idx);
                    end
                                   
				    sat_lat_AAo = zeros(NUM_TRIALS,NUM_STNS);sat_lon_AAo = zeros(NUM_TRIALS,NUM_STNS);
                    for m=1:NUM_TRIALS
                        [sat_lat_AAo(m,:),idx] = datasample(I_sat_AAo,NUM_STNS,'Replace',false);
                        sat_lon_AAo(m,:) = J_sat_AAo(idx);
                    end
                    save(['Proxies/NoENSO/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                    'stn_lat_AAo','stn_lon_AAo','I_AAo','J_AAo','sat_lat_AAo','sat_lon_AAo','I_sat_AAo','J_sat_AAo','-append');
                end
			elseif region == 6
                lat = SoA_lat; lon = SoA_lon;
                [I_SoA,J_SoA] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
                [I_sat_SoA,J_sat_SoA] = find(~isnan(corr_sat(lat,lon)));
                for NUM_STNS = 1:precip_min(floor(windowsize/30),6)
                    stn_lat_SoA = zeros(NUM_TRIALS,NUM_STNS); stn_lon_SoA = zeros(NUM_TRIALS,NUM_STNS); 
                    for m=1:NUM_TRIALS
                        [stn_lat_SoA(m,:),idx] = datasample(I_SoA,NUM_STNS,'Replace',false);
                        stn_lon_SoA(m,:) = J_SoA(idx);
                    end
                    save(['Proxies/NoENSO/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                    'stn_lat_SoA','stn_lon_SoA','I_SoA','J_SoA','-append');	
                end
				for NUM_STNS = 1:sat_min(floor(windowsize/30),6)
                    sat_lat_SoA = zeros(NUM_TRIALS,NUM_STNS); sat_lon_SoA = zeros(NUM_TRIALS,NUM_STNS);
                    for m=1:NUM_TRIALS
                        [sat_lat_SoA(m,:),idx] = datasample(I_sat_SoA,NUM_STNS,'Replace',false);
                        sat_lon_SoA(m,:) = J_sat_SoA(idx);
                    end
                    save(['Proxies/NoENSO/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                    'sat_lat_SoA','sat_lon_SoA','I_sat_SoA','J_sat_SoA','-append');	
                end
            end
        end
     c
  end
windowsize
end
 
 
 
 
 
 
