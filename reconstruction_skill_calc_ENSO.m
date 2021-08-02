% This script will calculate the range of the error in the calibration windows 
% for varying station numbers, using staggered windows. It also outputs a
% mat file containing alot of statistics produced for num_stns and calibration window.
% 
% Step 1: Run it with default setup (number of stations = 2:70)
% Step 2: Change NUM_STNS to 1 below in the loop, you also need to remove
% the transpose command here: squeeze(all_stn_precip(n,:,CAL_WDW(c,:)))');
% 
% This file calculates skills for reconsturctions without sites correlated with ENSO
%
% This should take roughly 18 hours to run in its entirety (all window sizes, all regions - run on a desktop i5.
%   
% Ryan Batehup, 2014
% Edited for calculation of SAM reconstruction using CPS only 
% by Willem Huiskamp, 2015

clear

%% Setup
tic;
load DataFiles/model_output.mat
load site_range_ENSO.mat

for windowsize = [31 61 91] % The running window in years
    for region = [1 2 3 4 5] %For the 91yr window, we must remove SA only
%% Loading proxies
    DIR_NAME = ['Proxies/NoENSO/',num2str(windowsize),'yrWindow'];
     % script fails at 1, it is noted in the script what to change if you only want one site
    NUM_YRS = 500; NUM_TRIALS = 1000;
    % Calibration windows set to being 10 overlapping windows over 500 years
    NUM_CAL_WDW = 10; clear CAL_WDW;
    overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
		for c=0:9
			CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
		end

%% Beginning the Loop - extracting the pseudoproxies from the precip file
		for c= 1:size(CAL_WDW,1)
            numstnstocompare = 2:precip_min(floor(windowsize/30),region);
            if max(numstnstocompare) > 70
                numstnstocompare = 2:70;
            end
			for NUM_STNS = numstnstocompare

				all_stn_precip=zeros(NUM_TRIALS,NUM_STNS,NUM_YRS);
                load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat']);
            
% Select region for analysis: 1 = Global (SH really), 2 = South America only, 3 = Australia + New Zealand, 4 = SH without Antarctica
% 5 = Antarctica only, 6 = South Africa only
				if region == 1
					lat_prox = stn_lat;
					lon_prox = stn_lon;
				elseif region == 2
					lat_prox = stn_lat_SA;
					lon_prox = stn_lon_SA;
				elseif region == 3
					lat_prox = stn_lat_AuNz;
					lon_prox = stn_lon_AuNz;
				elseif region == 4
					lat_prox = stn_lat_noAA;
					lon_prox = stn_lon_noAA;
				elseif region == 5
					lat_prox = stn_lat_AAo;
					lon_prox = stn_lon_AAo;
				elseif region == 6
					lat_prox = stn_lat_SoA;
					lon_prox = stn_lon_SoA;
				end
            			
				stn_precip = nan(NUM_TRIALS,NUM_STNS,NUM_YRS,'single');
                for m=1:NUM_TRIALS
					for n=1:NUM_STNS
						stn_precip(m,n,:) = single(precip_detr(:,lat_prox(m,n),lon_prox(m,n)));
                    end
				end

% Standardise the proxies

				stn_precip_mn=mean(stn_precip,3); 
				stn_precip_std=std(stn_precip,0,3);
                for n=1:NUM_TRIALS
					for m=1:NUM_STNS
						all_stn_precip(n,m,:) = single((stn_precip(n,m,:)-stn_precip_mn(n,m))./(stn_precip_std(n,m)));
                    end
				end

				clear lat_prox lon_prox stn_precip stn_precip_mn stn_precip_std

%% Esper et al 2005 CPS Method

                for n=1:NUM_TRIALS
					corr_matrix = corr(SAM(CAL_WDW(c,:))*ones(1,NUM_STNS), squeeze(all_stn_precip(n,:,CAL_WDW(c,:)))'); % If you want to calculate for one site, change NUM_STNS manually and run this last bit once more.
					stn_CPS(n,:) = corr_matrix(1,:)*squeeze(all_stn_precip(n,:,:));
                end

% Normalising (it already has mean ~0)
				for n=1:NUM_TRIALS
					stn_CPS(n,:) = single(squeeze(stn_CPS(n,:))./std(squeeze(stn_CPS(n,:))'));
                end

% Skill Evaluation
				for n=1:NUM_TRIALS
					stn_corr_CPS(n) = single(corr(squeeze(stn_CPS(n,:)'),SAM));
					stn_rmse_CPS(n) = single(sqrt(mean((SAM'-squeeze(stn_CPS(n,:))).^2)));
                end

				all_stn_CPS(NUM_STNS,:,:) = stn_CPS;
				all_stn_corr_CPS(NUM_STNS,:) = stn_corr_CPS;
				all_stn_rmse_CPS(NUM_STNS,:) = stn_rmse_CPS;
            end
            
            numstnstocompare = 2:sat_min(floor(windowsize/30),region);
            if max(numstnstocompare) > 70
                numstnstocompare = 2:70;
            end
			for NUM_STNS = numstnstocompare
                all_stn_sat=zeros(NUM_TRIALS,NUM_STNS,NUM_YRS);
				load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat']);
            
% Select region for analysis: 1 = Global (SH really), 2 = South America only, 3 = Australia + New Zealand, 4 = SH without Antarctica
% 5 = Antarctica only, 6 = South Africa only
				if region == 1
					lat_prox = sat_lat;
					lon_prox = sat_lon;
				elseif region == 2
					lat_prox = sat_lat_SA;
					lon_prox = sat_lon_SA;
				elseif region == 3
					lat_prox = sat_lat_AuNz;
					lon_prox = sat_lon_AuNz;
				elseif region == 4
					lat_prox = sat_lat_noAA;
					lon_prox = sat_lon_noAA;
				elseif region == 5
					lat_prox = sat_lat_AAo;
					lon_prox = sat_lon_AAo;
				elseif region == 6
					lat_prox = sat_lat_SoA;
					lon_prox = sat_lon_SoA;
				end
            			
                stn_sat = nan(NUM_TRIALS,NUM_STNS,NUM_YRS,'single');
				for m=1:NUM_TRIALS
					for n=1:NUM_STNS
                        stn_sat(m,n,:) = single(sat_detr(:,lat_prox(m,n),lon_prox(m,n)));
					end
				end

% Standardise the proxies

                stn_sat_mn=mean(stn_sat,3); 
				stn_sat_std=std(stn_sat,0,3);
                
				for n=1:NUM_TRIALS
					for m=1:NUM_STNS
                        all_stn_sat(n,m,:) = single((stn_sat(n,m,:)-stn_sat_mn(n,m))./(stn_sat_std(n,m)));
					end
				end

				clear lat_prox lon_prox stn_sat stn_sat_mn stn_sat_std

%% Esper et al 2005 CPS Method

                for n=1:NUM_TRIALS
                    corr_matrix_sat = corr(SAM(CAL_WDW(c,:))*ones(1,NUM_STNS), squeeze(all_stn_sat(n,:,CAL_WDW(c,:)))');
					sat_CPS(n,:) = corr_matrix_sat(1,:)*squeeze(all_stn_sat(n,:,:));
				end

% Normalising (it already has mean ~0)
				for n=1:NUM_TRIALS
                    sat_CPS(n,:) = single(squeeze(sat_CPS(n,:))./std(squeeze(sat_CPS(n,:))'));
				end

% Skill Evaluation
				for n=1:NUM_TRIALS
                    sat_corr_CPS(n) = single(corr(squeeze(sat_CPS(n,:)'),SAM));
					sat_rmse_CPS(n) = single(sqrt(mean((SAM'-squeeze(sat_CPS(n,:))).^2)));
				end
                
                all_sat_CPS(NUM_STNS,:,:) = sat_CPS;
				all_sat_corr_CPS(NUM_STNS,:) = sat_corr_CPS;
				all_sat_rmse_CPS(NUM_STNS,:) = sat_rmse_CPS;
            end           
            
% Save each region to file
			if region == 1
				save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
				'all_stn_CPS','all_stn_corr_CPS','all_stn_rmse_CPS','all_sat_CPS','all_sat_corr_CPS','all_sat_rmse_CPS');
			elseif region == 2
				all_stn_CPS_SA = all_stn_CPS;
				all_stn_corr_CPS_SA = all_stn_corr_CPS;
				all_stn_rmse_CPS_SA = all_stn_rmse_CPS;
                all_sat_CPS_SA = all_sat_CPS;
				all_sat_corr_CPS_SA = all_sat_corr_CPS;
				all_sat_rmse_CPS_SA = all_sat_rmse_CPS;
				save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
				'all_stn_CPS_SA','all_stn_corr_CPS_SA','all_stn_rmse_CPS_SA','all_sat_CPS_SA','all_sat_corr_CPS_SA','all_sat_rmse_CPS_SA','-append');
			elseif region == 3
				all_stn_CPS_AuNz = all_stn_CPS;
				all_stn_corr_CPS_AuNz = all_stn_corr_CPS;
				all_stn_rmse_CPS_AuNz = all_stn_rmse_CPS;
                all_sat_CPS_AuNz = all_sat_CPS;
				all_sat_corr_CPS_AuNz = all_sat_corr_CPS;
				all_sat_rmse_CPS_AuNz = all_sat_rmse_CPS;
				save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
				'all_stn_CPS_AuNz','all_stn_corr_CPS_AuNz','all_stn_rmse_CPS_AuNz','all_sat_CPS_AuNz','all_sat_corr_CPS_AuNz','all_sat_rmse_CPS_AuNz','-append');
			elseif region == 4
				all_stn_CPS_AA = all_stn_CPS;
				all_stn_corr_CPS_AA = all_stn_corr_CPS;
				all_stn_rmse_CPS_AA = all_stn_rmse_CPS;
                all_sat_CPS_AA = all_sat_CPS;
				all_sat_corr_CPS_AA = all_sat_corr_CPS;
				all_sat_rmse_CPS_AA = all_sat_rmse_CPS;
				save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
				'all_stn_CPS_AA','all_stn_corr_CPS_AA','all_stn_rmse_CPS_AA','all_sat_CPS_AA','all_sat_corr_CPS_AA','all_sat_rmse_CPS_AA','-append');
			elseif region == 5
				all_stn_CPS_AAo = all_stn_CPS;
				all_stn_corr_CPS_AAo = all_stn_corr_CPS;
				all_stn_rmse_CPS_AAo = all_stn_rmse_CPS;
                all_sat_CPS_AAo = all_sat_CPS;
				all_sat_corr_CPS_AAo = all_sat_corr_CPS;
				all_sat_rmse_CPS_AAo = all_sat_rmse_CPS;
				save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
				'all_stn_CPS_AAo','all_stn_corr_CPS_AAo','all_stn_rmse_CPS_AAo','all_sat_CPS_AAo','all_sat_corr_CPS_AAo','all_sat_rmse_CPS_AAo','-append');
			elseif region == 6
				all_stn_CPS_SoA = all_stn_CPS;
				all_stn_corr_CPS_SoA = all_stn_corr_CPS;
				all_stn_rmse_CPS_SoA = all_stn_rmse_CPS;
                all_sat_CPS_SoA = all_sat_CPS;
				all_sat_corr_CPS_SoA = all_sat_corr_CPS;
				all_sat_rmse_CPS_SoA = all_sat_rmse_CPS;
				save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
				'all_stn_CPS_SoA','all_stn_corr_CPS_SoA','all_stn_rmse_CPS_SoA','all_sat_CPS_SoA','all_sat_corr_CPS_SoA','all_sat_rmse_CPS_SoA','-append');	
            end
            clear all_stn_CPS all_stn_corr_CPS all_stn_rmse_CPS all_sat_CPS all_sat_corr_CPS all_sat_rmse_CPS
		%c
		end
	region
    end
windowsize
end
