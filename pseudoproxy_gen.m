% Script to generate pseudoproxy groups
% Willem Huiskamp, 2015
%
% Method roughly follows that of Ryan Batehup (Honours Thesis).
% Proxies are only selected from land tiles and must have an absolute
% correlation with SAM greater than 0.3 over a 20 year calibration window
% (or over the entire time series). We then generate a network of proxies
% between one and 20 sites (or more).

clear

% load land mask and setup new variables
land = nc_varget('sftlf_A1.static.nc','sftlf');
load('DataFiles/model_output','precip_detr','sat_detr','SAM')

% Next we need to create 10 evenly spaced calibration windows in the data,
% of the length of our correlation windows.

NUM_SYNRUNS = 1000; NUM_YRS = 500; NUM_TRIALS = 1000;
NUM_CAL_WDW = 10;
STN_MAX = 70;

for windowsize = [31 61 91]
  
  if windowsize > 31
      clear CAL_WDW
  end
  
  for c=0:9
       overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
       CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
  end
     
% Load topography and eliminate all data that isn't land. This will save time with the correlations.

  for t = 1:500
    for i = 1:size(land,1)
      for j = 1:size(land,2)
        if land(i,j) > 0 
               precip_land(:,i,j) = precip_detr(:,i,j);
          else precip_land(:,i,j) = NaN;
        end
      end
    end
  end 
  
  for t = 1:500
    for i = 1:size(land,1)
      for j = 1:size(land,2)
        if land(i,j) > 0 
               sat_land(:,i,j) = sat_detr(:,i,j);
          else sat_land(:,i,j) = NaN;
        end
      end
    end
  end 
  
% Calculate correlations 

  for c = 1:10
    for i = 1:90
        for j = 1:144
            corr_precip(i,j) = corr(SAM(CAL_WDW(c,:))', precip_land(CAL_WDW(c,:),i,j));
            corr_sat(i,j) = corr(SAM(CAL_WDW(c,:))', sat_land(CAL_WDW(c,:),i,j));
        end
    end
    mkdir(['Proxies/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end))]);


    % Mask correlations so we only have land values and correlations greater
    % than or equal to +-0.3

    for i = 1:size(land,1)
        for j = 1:size(land,2)
            if abs(corr_precip(i,j)) < 0.3
            corr_precip(i,j) = NaN;
            end
        end
    end
    
    for i = 1:size(land,1)
        for j = 1:size(land,2)
            if abs(corr_sat(i,j)) < 0.3
            corr_sat(i,j) = NaN;
            end
        end
    end
    % Time to generate or random proxy group

    for NUM_STNS = 1:STN_MAX
    
    % Make Global (SH)=1, SA=2, Au/NZ=3, SH without AA = 4, Antarctica (coast) = 5, Antarctica (interior) = 6 
        for region = [4 5 6]
    
            stn_lat = zeros(NUM_TRIALS,NUM_STNS); sat_lat = zeros(NUM_TRIALS,NUM_STNS);
            stn_lon = zeros(NUM_TRIALS,NUM_STNS); sat_lon = zeros(NUM_TRIALS,NUM_STNS);
            stn_lat_SA = zeros(NUM_TRIALS,NUM_STNS); sat_lat_SA = zeros(NUM_TRIALS,NUM_STNS);
            stn_lon_SA = zeros(NUM_TRIALS,NUM_STNS); sat_lon_SA = zeros(NUM_TRIALS,NUM_STNS);
            stn_lat_AuNz = zeros(NUM_TRIALS,NUM_STNS); sat_lat_AuNz = zeros(NUM_TRIALS,NUM_STNS);
            stn_lon_AuNz = zeros(NUM_TRIALS,NUM_STNS); sat_lon_AuNz = zeros(NUM_TRIALS,NUM_STNS);
			stn_lat_noAA = zeros(NUM_TRIALS,NUM_STNS); sat_lat_noAA = zeros(NUM_TRIALS,NUM_STNS);
			stn_lon_noAA = zeros(NUM_TRIALS,NUM_STNS); sat_lon_noAA = zeros(NUM_TRIALS,NUM_STNS);
			stn_lat_AAc = zeros(NUM_TRIALS,NUM_STNS); sat_lat_AAc = zeros(NUM_TRIALS,NUM_STNS);
			stn_lon_AAc = zeros(NUM_TRIALS,NUM_STNS); sat_lon_AAc = zeros(NUM_TRIALS,NUM_STNS);
			stn_lat_AAi = zeros(NUM_TRIALS,NUM_STNS); sat_lat_AAi = zeros(NUM_TRIALS,NUM_STNS);
			stn_lon_AAi = zeros(NUM_TRIALS,NUM_STNS); sat_lon_AAi = zeros(NUM_TRIALS,NUM_STNS);
    
            Gl_lon = 1:144; Gl_lat = 1:45;
            SA_lon = 110:130; SA_lat = 15:45;
            Au_lon = 40:75; Au_lat = 15:41;
			AA_lon = 1:144; AA_lat = 15:45;
			Ac_lon = 1:144; Ac_lat = 9:15;
			Ai_lon = 1:144; Ai_lat = 1:8;
        
            if region == 1
                lat = Gl_lat; lon = Gl_lon;
        
                [I,J] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
				
                for m=1:NUM_TRIALS
                    [stn_lat(m,:),idx] = datasample(I,NUM_STNS);
                    stn_lon(m,:) = J(idx);
                end
				
				[I_sat,J_sat] = find(~isnan(corr_sat(lat,lon)));
				
				for m=1:NUM_TRIALS
                    [sat_lat(m,:),idx] = datasample(I_sat,NUM_STNS);
                    sat_lon(m,:) = J_sat(idx);
                end
				
                save(['Proxies/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                    'stn_lat','stn_lon','I','J','corr_precip','sat_lat','sat_lon','I_sat','J_sat','corr_sat','windowsize');       
            elseif region == 2
                lat = SA_lat; lon = SA_lon;
                [I_SA,J_SA] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
    
                for m=1:NUM_TRIALS
                    [stn_lat_SA(m,:),idx] = datasample(I_SA,NUM_STNS);
                    stn_lon_SA(m,:) = J_SA(idx);
                end
				
				[I_sat_SA,J_sat_SA] = find(~isnan(corr_sat(lat,lon)));
				
				for m=1:NUM_TRIALS
                    [sat_lat_SA(m,:),idx] = datasample(I_sat_SA,NUM_STNS);
                    sat_lon_SA(m,:) = J_sat_SA(idx);
                end
				
                save(['Proxies/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                    'stn_lat_SA','stn_lon_SA','I_SA','J_SA','sat_lat_SA','sat_lon_SA','I_sat_SA','J_sat_SA','-append');
			elseif region == 3
                lat = Au_lat; lon = Au_lon;
                [I_Au,J_Au] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
    
                for m=1:NUM_TRIALS
                    [stn_lat_AuNz(m,:),idx] = datasample(I_Au,NUM_STNS);
                    stn_lon_AuNz(m,:) = J_Au(idx);
                end
    
				[I_sat_Au,J_sat_Au] = find(~isnan(corr_sat(lat,lon)));
				
				for m=1:NUM_TRIALS
                    [sat_lat_Au(m,:),idx] = datasample(I_sat_Au,NUM_STNS);
                    sat_lon_Au(m,:) = J_sat_Au(idx);
                end
	
                save(['Proxies/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                    'stn_lat_AuNz','stn_lon_AuNz','I_Au','J_Au','sat_lat_AuNz','sat_lon_AuNz','I_sat_Au','J_sat_Au','-append');
			elseif region == 4
                lat = AA_lat; lon = AA_lon;
                [I_AA,J_AA] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
    
                for m=1:NUM_TRIALS
                    [stn_lat_noAA(m,:),idx] = datasample(I_AA,NUM_STNS);
                    stn_lon_noAA(m,:) = J_AA(idx);
                end
    
				[I_sat_AA,J_sat_AA] = find(~isnan(corr_sat(lat,lon)));
				
				for m=1:NUM_TRIALS
                    [sat_lat_AA(m,:),idx] = datasample(I_sat_AA,NUM_STNS);
                    sat_lon_AA(m,:) = J_sat_AA(idx);
                end
	
                save(['Proxies/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                    'stn_lat_noAA','stn_lon_noAA','I_AA','J_AA','sat_lat_noAA','sat_lon_noAA','I_sat_AA','J_sat_AA','-append');
            elseif region == 5
                lat = Ac_lat; lon = Ac_lon;
                [I_Ac,J_Ac] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
    
                for m=1:NUM_TRIALS
                    [stn_lat_AAc(m,:),idx] = datasample(I_Ac,NUM_STNS);
                    stn_lon_AAc(m,:) = J_Ac(idx);
                end
    
				[I_sat_Ac,J_sat_Ac] = find(~isnan(corr_sat(lat,lon)));
				
				for m=1:NUM_TRIALS
                    [sat_lat_AAc(m,:),idx] = datasample(I_sat_Ac,NUM_STNS);
                    sat_lon_AAc(m,:) = J_sat_Ac(idx);
                end
	
                save(['Proxies/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                    'stn_lat_AAc','stn_lon_AAc','I_Ac','J_Ac','sat_lat_AAc','sat_lon_AAc','I_sat_Ac','J_sat_Ac','-append');
			elseif region == 6
                lat = Ai_lat; lon = Ai_lon;
                [I_Ai,J_Ai] = find(~isnan(corr_precip(lat,lon))); % I is lat, J in lon
    
                for m=1:NUM_TRIALS
                    [stn_lat_AAi(m,:),idx] = datasample(I_Ai,NUM_STNS);
                    stn_lon_AAi(m,:) = J_Ai(idx);
                end
    
				[I_sat_Ai,J_sat_Ai] = find(~isnan(corr_sat(lat,lon)));
				
				for m=1:NUM_TRIALS
                    [sat_lat_AAi(m,:),idx] = datasample(I_sat_Ai,NUM_STNS);
                    sat_lon_AAi(m,:) = J_sat_Ai(idx);
                end
				
                save(['Proxies/',num2str(windowsize),'yrWindow/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat'], ...
                    'stn_lat_AAi','stn_lon_AAi','I_Ai','J_Ai','sat_lat_AAi','sat_lon_AAi','I_sat_Ai','J_sat_Ai','-append');		
		
            end 
		end
    end
  c
  end
windowsize
end
 


 
 
 
 
 
 
 