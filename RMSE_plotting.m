%% Plot the median RMSE for our reconstructions

clear all
NUM_YRS = 500; NUM_TRIALS = 1000;
NUM_CAL_WDW = 10; 
STN_MAX = 70;

for windowsize = [31 61 91]
    if windowsize > 31
      clear CAL_WDW
    end
    rmse_sat_SH = nan(NUM_CAL_WDW,STN_MAX,NUM_TRIALS);
    rmse_sat_AA = nan(NUM_CAL_WDW,STN_MAX,NUM_TRIALS);
    rmse_sat_AAo = nan(NUM_CAL_WDW,STN_MAX,NUM_TRIALS);
    rmse_sat_AuNz = nan(NUM_CAL_WDW,STN_MAX,NUM_TRIALS);
    rmse_sat_SA = nan(NUM_CAL_WDW,STN_MAX,NUM_TRIALS);
    rmse_stn_SH = nan(NUM_CAL_WDW,STN_MAX,NUM_TRIALS);
    rmse_stn_AA = nan(NUM_CAL_WDW,STN_MAX,NUM_TRIALS);
    rmse_stn_AAo = nan(NUM_CAL_WDW,STN_MAX,NUM_TRIALS);
    rmse_stn_AuNz = nan(NUM_CAL_WDW,STN_MAX,NUM_TRIALS);
    rmse_stn_SA = nan(NUM_CAL_WDW,STN_MAX,NUM_TRIALS);
    
    for c = 0:9
       overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
       CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
    end
    
    for c = 1:NUM_CAL_WDW
        DIR = ['Proxies/NoResample/',num2str(windowsize),'yrWindow/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'];
        load(DIR,'all_stn_rmse_CPS','all_stn_rmse_CPS_AA','all_stn_rmse_CPS_AAo','all_stn_rmse_CPS_AuNz','all_stn_rmse_CPS_SA',...
            'all_sat_rmse_CPS','all_sat_rmse_CPS_AA','all_sat_rmse_CPS_AAo','all_sat_rmse_CPS_AuNz','all_sat_rmse_CPS_SA')
        % SAT
        rmse_sat_SH(c,:,:) = all_sat_rmse_CPS;
        rmse_sat_AA(c,:,:) = all_sat_rmse_CPS_AA;
        rmse_sat_AAo(c,:,:) = all_sat_rmse_CPS_AAo;
        rmse_sat_AuNz(c,1:size(all_sat_rmse_CPS_AuNz,1),:) = all_sat_rmse_CPS_AuNz;
        rmse_sat_SA(c,1:size(all_sat_rmse_CPS_SA,1),:) = all_sat_rmse_CPS_SA;
        % Precip
        rmse_stn_SH(c,:,:) = all_stn_rmse_CPS;
        rmse_stn_AA(c,:,:) = all_stn_rmse_CPS_AA;
        rmse_stn_AAo(c,:,:) = all_stn_rmse_CPS_AAo;
        rmse_stn_AuNz(c,1:size(all_stn_rmse_CPS_AuNz,1),:) = all_stn_rmse_CPS_AuNz;
        rmse_stn_SA(c,1:size(all_stn_rmse_CPS_SA,1),:) = all_stn_rmse_CPS_SA;
        clear all_stn_rmse_CPS all_stn_rmse_CPS_AA all_stn_rmse_CPS_AAo all_stn_rmse_CPS_AuNz all_stn_rmse_CPS_SA ...
              all_sat_rmse_CPS all_sat_rmse_CPS_AA all_sat_rmse_CPS_AAo all_sat_rmse_CPS_AuNz all_sat_rmse_CPS_SA 
    end
    
    rmse_sat_SH_M = median(median(rmse_sat_SH,3),1);
    rmse_sat_AA_M = median(median(rmse_sat_AA,3),1);
    rmse_sat_AAo_M = median(median(rmse_sat_AAo,3),1);
    rmse_sat_AuNz_M = median(median(rmse_sat_AuNz,3),1);
    rmse_sat_SA_M = median(median(rmse_sat_SA,3),1);
   
    rmse_stn_SH_M = median(median(rmse_stn_SH,3),1);
    rmse_stn_AA_M = median(median(rmse_stn_AA,3),1);
    rmse_stn_AAo_M = median(median(rmse_stn_AAo,3),1);
    rmse_stn_AuNz_M = median(median(rmse_stn_AuNz,3),1);
    rmse_stn_SA_M = median(median(rmse_stn_SA,3),1);
    
    save(['RMSE_',num2str(windowsize),'yrWdw_recons.mat'],'rmse_sat_SH_M','rmse_sat_AA_M','rmse_sat_AAo_M','rmse_sat_AuNz_M','rmse_sat_SA_M',...
        'rmse_stn_SH_M','rmse_stn_AA_M','rmse_stn_AAo_M','rmse_stn_AuNz_M','rmse_stn_SA_M')
        
end


%% Making the plots

clear all

figure(1)
for windowsize = [31 61 91]
    load(['RMSE_',num2str(windowsize),'yrWdw_recons.mat'])
    subplot(3,3,floor(windowsize/30))
    grid on
    hold on
    h1 = plot(rmse_sat_SH_M(2:70),'linewidth',1.5);
    %h2 = plot(rmse_sat_AA_M(2:70),'linewidth',1.5);
    h3 = plot(rmse_sat_AAo_M(2:70),'linewidth',1.5);
    h4 = plot(rmse_sat_AuNz_M(2:70),'linewidth',1.5);
    h5 = plot(rmse_sat_SA_M(2:70),'linewidth',1.5);
    xlim([0 70]); ylim([0.88 1])
    leg = legend([h1 h3 h4 h5],...
        'SH','AA only','AuNz','SA');
end
for windowsize = [31 61 91]
    load(['RMSE_',num2str(windowsize),'yrWdw_recons.mat'])
    subplot(3,3,floor(windowsize/30)+3)
    grid on
    hold on
    h6 = plot(rmse_stn_SH_M(2:70),'linewidth',1.5);
    %h7 = plot(rmse_stn_AA_M(2:70),'linewidth',1.5);
    h8 = plot(rmse_stn_AAo_M(2:70),'linewidth',1.5);
    h9 = plot(rmse_stn_AuNz_M(2:70),'linewidth',1.5);
    h10 = plot(rmse_stn_SA_M(2:70),'linewidth',1.5);
    xlim([0 70]); ylim([0.88 1])
    leg = legend([h6 h8 h9 h10],...
        'SH','AA only','AuNz','SA');
end
subplot(3,3,1); ylabel('RMSE')
subplot(3,3,4); ylabel('RMSE')
    
subplot(3,3,5); xlabel('Number of proxies');     
    
subplot(3,3,2); title('SAT'); subplot(3,3,5); title('Precip.');  

print('RMSE.pdf','-dpdf','-painters')    
    
    
    
    
    
    
    
    
    
    

