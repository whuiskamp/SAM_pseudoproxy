% This script will calculate the range of the error in the calibration windows 
% for varying station numbers, using staggered windows. It also outputs a
% mat file containing alot of statistics produced for num_stns and calibration window.
% 
% Step 1: Run it with default setup (number of stations = 2:70)
% Step 2: Change NUM_STNS to 1 below in the loop, you also need to remove
% the transpose command here: squeeze(all_stn_precip(n,:,CAL_WDW(c,:)))');
%
% Ryan Batehup, 2014
% Edited for calculation of SAM reconstruction using CPS only 
% by Willem Huiskamp, 2015

clear

%% Setup

load DataFiles/model_output.mat

windowsize = 500; % The running window in years

%% Loading proxies
DIR_NAME = ['Proxies/NoSkill/',num2str(windowsize),'yrWindow']; 
numstnstocompare = 2:70; % script fails at 1, it is noted in the script what to change if you only want one site
NUM_YRS = 500; NUM_TRIALS = 1000;
% Calibration windows set to being 10 overlapping windows over 500 years
NUM_CAL_WDW = 1; clear CAL_WDW;

CAL_WDW(1,:) = 1:500;


%% Beginning the Loop
for c=1:size(CAL_WDW,1)
tic;
    %all_stn_CPS = nan(max(numstnstocompare),NUM_TRIALS, NUM_YRS,'single');
    %all_stn_corr_CPS = nan(max(numstnstocompare),NUM_TRIALS,'single');
    %all_stn_rmse_CPS = nan(max(numstnstocompare),NUM_TRIALS,'single');
     
    for NUM_STNS = numstnstocompare %1

        all_stn_precip=zeros(NUM_TRIALS,NUM_STNS,NUM_YRS);
        all_stn_sat=zeros(NUM_TRIALS,NUM_STNS,NUM_YRS);
        
        load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(NUM_STNS),'stns_1000prox.mat']);
        all_stn_lat=stn_lat;
        all_stn_lon=stn_lon;
        
        all_sat_lat=sat_lat;
        all_sat_lon=sat_lon;

        stn_precip = nan(NUM_TRIALS,NUM_STNS,NUM_YRS,'single');
        stn_sat = nan(NUM_TRIALS,NUM_STNS,NUM_YRS,'single');
        for m=1:NUM_TRIALS
            for n=1:NUM_STNS
                stn_precip(m,n,:) = single(precip_detr(:,stn_lat(m,n),stn_lon(m,n)));
                stn_sat(m,n,:) = single(sat_detr(:,sat_lat(m,n),sat_lon(m,n)));
            end
        end

        stn_precip_mn=mean(stn_precip,3); 
        stn_precip_std=std(stn_precip,0,3);
        
        stn_sat_mn=mean(stn_sat,3); 
        stn_sat_std=std(stn_sat,0,3);
        
        for n=1:NUM_TRIALS
            for m=1:NUM_STNS
                all_stn_precip(n,m,:) = single((stn_precip(n,m,:)-stn_precip_mn(n,m))./(stn_precip_std(n,m)));
                all_stn_sat(n,m,:) = single((stn_sat(n,m,:)-stn_sat_mn(n,m))./(stn_sat_std(n,m)));
            end
        end

    clear stn_lat stn_lon stn_precip stn_precip_mn stn_precip_std sat_lat sat_lon stn_sat stn_sat_mn stn_sat_std

%% Esper et al 2005 CPS Method

        %stn_CPS = nan(NUM_TRIALS, NUM_YRS);
        %stn_corr_CPS = nan(NUM_TRIALS,1);
        %stn_rmse_CPS = nan(NUM_TRIALS,1);
        
        for n=1:NUM_TRIALS
            corr_matrix = corr(SAM(CAL_WDW(c,:))'*ones(1,NUM_STNS), squeeze(all_stn_precip(n,:,CAL_WDW(c,:)))'); % If you want to calculate for one site, change NUM_STNS manually and run this last bit once more.
            stn_CPS(n,:) = corr_matrix(1,:)*squeeze(all_stn_precip(n,:,:));
            
            corr_matrix_sat = corr(SAM(CAL_WDW(c,:))'*ones(1,NUM_STNS), squeeze(all_stn_sat(n,:,CAL_WDW(c,:)))'); % If you want to calculate for one site, change NUM_STNS manually and run this last bit once more.
            sat_CPS(n,:) = corr_matrix_sat(1,:)*squeeze(all_stn_sat(n,:,:));
        end

% Normalising (it already has mean ~0)
        for n=1:NUM_TRIALS
            stn_CPS(n,:) = single(squeeze(stn_CPS(n,:))./std(squeeze(stn_CPS(n,:))'));
            sat_CPS(n,:) = single(squeeze(sat_CPS(n,:))./std(squeeze(sat_CPS(n,:))'));
        end

% Skill Evaluation
        for n=1:NUM_TRIALS
            stn_corr_CPS(n) = single(corr(squeeze(stn_CPS(n,:)'),SAM'));
            stn_rmse_CPS(n) = single(sqrt(mean((SAM-squeeze(stn_CPS(n,:))).^2)));
            sat_corr_CPS(n) = single(corr(squeeze(sat_CPS(n,:)'),SAM'));
            sat_rmse_CPS(n) = single(sqrt(mean((SAM-squeeze(sat_CPS(n,:))).^2)));
        end

        all_stn_CPS(NUM_STNS,:,:) = stn_CPS;
        all_stn_corr_CPS(NUM_STNS,:) = stn_corr_CPS;
        all_stn_rmse_CPS(NUM_STNS,:) = stn_rmse_CPS;
        
        all_sat_CPS(NUM_STNS,:,:) = sat_CPS;
        all_sat_corr_CPS(NUM_STNS,:) = sat_corr_CPS;
        all_sat_rmse_CPS(NUM_STNS,:) = sat_rmse_CPS;

    end
% That took about an hour per calibration window
save([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
     'all_stn_CPS','all_stn_corr_CPS','all_stn_rmse_CPS','all_sat_CPS','all_sat_corr_CPS','all_sat_rmse_CPS');
 toc;
 
end
