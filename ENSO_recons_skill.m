%% Figures for recon skill comparison between all sites and ENSO sites removed.

clear
numstnstocompare = 2:70; skilful_threshold = sqrt(0.5); % Using a linear statistical model, we want r^2 (Coefficient of determination) to be greater than 0.5, so r must be greater than sqrt(0.5).
NUM_CAL_WDW = 10; clear CAL_WDW; NUM_YRS = 500; NUM_TRIALS = 1000; 

pcnt_skilful_CPS = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that do better than the skill threshold
pcnt_skilful_CPS_SH = nan(max(numstnstocompare),length([31 61 91]));

corr_CPS_qn = nan(70,10,3);
corr_CPS_qn_SH = nan(70,10,3);

corr_CPS_qn_rng = nan(3,size(corr_CPS_qn,1),2,size(corr_CPS_qn,3)); % creates the range for the 3 windows, 70 sites, 2 bounds of the range, 3 percentiles
corr_CPS_qn_rng_SH = nan(3,size(corr_CPS_qn_SH,1),2,size(corr_CPS_qn_SH,3)); 

SAT_corr_CPS = nan(3,70,10,1000);
SAT_corr_CPS_SH = nan(3,70,10,1000);

% For these figures, instead of comparing to the 500 year recons, we compare to the recons with ENSO sites

for windowsize = [31 61 91]
    DIR_NAME = ['Proxies/NoResample/',num2str(windowsize),'yrWindow'];
    DIR_NAME_E = ['Proxies/NoENSO/',num2str(windowsize),'yrWindow'];
    overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
        clear CAL_WDW;
        for c=0:9
            CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
        end
        
        for c=1:size(CAL_WDW,1)
            load([DIR_NAME_E,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
            'all_sat_corr_CPS')
            SAT_corr_CPS_SH(floor(windowsize/30),:,c,:) = all_sat_corr_CPS; % this is (3,70,10,1000)
			load([DIR_NAME,'/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
            'all_sat_corr_CPS')
            SAT_corr_CPS(floor(windowsize/30),:,c,:) = all_sat_corr_CPS;
        end

               
        for n = numstnstocompare
            pcnt_skilful_CPS_SH(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(SAT_corr_CPS_SH(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
            
            pcnt_skilful_CPS(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(SAT_corr_CPS(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end

        % Southern Hemisphere (normal) recon.
        corr_CPS_qn = quantile(squeeze(SAT_corr_CPS(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng(floor(windowsize/30),:,1,:) = min(corr_CPS_qn,[],2);
        corr_CPS_qn_rng(floor(windowsize/30),:,2,:) = max(corr_CPS_qn,[],2);
        
        % Southern Hemisphere (ENSO) recon.
        corr_CPS_qn_SH = quantile(squeeze(SAT_corr_CPS_SH(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_SH(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_SH,[],2);
        corr_CPS_qn_rng_SH(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_SH,[],2);
end


%% Plotting
figure(1)
clf
s_Hnd = tight_subplot(2,3,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

for windowsize = 31
    axes(s_Hnd(1)); hold on; xlim([0,70]); ylim([0,1]); grid on
    jbfill(2:70,squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    plot(pcnt_skilful_CPS(:,floor(windowsize/30)),'Color','k','LineWidth',3); 
    
    axes(s_Hnd(1));
    jbfill(2:70,squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'r','k',[],0.5);
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
end

for windowsize = 61
    axes(s_Hnd(2)); hold on; xlim([0,70]); ylim([0,1]); grid on
    jbfill(2:70,squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    plot(pcnt_skilful_CPS(:,floor(windowsize/30)),'Color','k','LineWidth',3); 
    
    axes(s_Hnd(2));
    jbfill(2:70,squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'r','k',[],0.5);
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
end

for windowsize = 91
    axes(s_Hnd(3)); hold on; xlim([0,70]); ylim([0,1]); grid on
    jbfill(2:70,squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    plot(pcnt_skilful_CPS(:,floor(windowsize/30)),'Color','k','LineWidth',3); 
    
    axes(s_Hnd(3));
    jbfill(2:70,squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'r','k',[],0.5);
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
end
  
numstnstocompare = 2:70; skilful_threshold = sqrt(0.5); % Using a linear statistical model, we want r^2 (Coefficient of determination) to be greater than 0.5, so r must be greater than sqrt(0.5).
NUM_CAL_WDW = 10; clear CAL_WDW; NUM_YRS = 500; NUM_TRIALS = 1000; 

pcnt_skilful_CPS = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that do better than the skill threshold
pcnt_skilful_CPS_SH = nan(max(numstnstocompare),length([31 61 91]));

corr_CPS_qn = nan(70,10,3);
corr_CPS_qn_SH = nan(70,10,3);

corr_CPS_qn_rng = nan(3,size(corr_CPS_qn,1),2,size(corr_CPS_qn,3)); % creates the range for the 3 windows, 70 sites, 2 bounds of the range, 3 percentiles
corr_CPS_qn_rng_SH = nan(3,size(corr_CPS_qn_SH,1),2,size(corr_CPS_qn_SH,3)); 

precip_corr_CPS = nan(3,70,10,1000);
precip_corr_CPS_SH = nan(3,70,10,1000);
        
for windowsize = [31 61 91]
    DIR_NAME = ['Proxies/NoResample/',num2str(windowsize),'yrWindow'];
    DIR_NAME_E = ['Proxies/NoENSO/',num2str(windowsize),'yrWindow'];
    overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
        clear CAL_WDW;
        for c=0:9
            CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
        end
        
        for c=1:size(CAL_WDW,1)
            load([DIR_NAME_E,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
            'all_stn_corr_CPS')
            precip_corr_CPS_SH(floor(windowsize/30),:,c,:) = all_stn_corr_CPS; % this is (3,70,10,1000)
			load([DIR_NAME,'/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
            'all_stn_corr_CPS')
            precip_corr_CPS(floor(windowsize/30),:,c,:) = all_stn_corr_CPS;
        end

        for n = numstnstocompare
            pcnt_skilful_CPS_SH(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS_SH(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);   
        end
        for n = numstnstocompare
            pcnt_skilful_CPS(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end

        % Southern Hemisphere (normal) recon.
        corr_CPS_qn = quantile(squeeze(precip_corr_CPS(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng(floor(windowsize/30),:,1,:) = min(corr_CPS_qn,[],2);
        corr_CPS_qn_rng(floor(windowsize/30),:,2,:) = max(corr_CPS_qn,[],2);
        
        % Southern Hemisphere (ENSO) recon.
        corr_CPS_qn_SH = quantile(squeeze(precip_corr_CPS_SH(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_SH(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_SH,[],2);
        corr_CPS_qn_rng_SH(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_SH,[],2);
        clear corr_CPS_qn_SH all_stn_corr_CPS_SH
end
        
% Plotting
for windowsize = 31
    axes(s_Hnd(4)); hold on; xlim([0,70]); ylim([0,1]); grid on
    jbfill(2:70,squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    plot(pcnt_skilful_CPS(:,floor(windowsize/30)),'Color','k','LineWidth',3); 
    
    axes(s_Hnd(4));
    jbfill(2:70,squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'r','k',[],0.5);
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
end

for windowsize = 61
    axes(s_Hnd(5)); hold on; xlim([0,70]); ylim([0,1]); grid on
    jbfill(2:70,squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    plot(pcnt_skilful_CPS(:,floor(windowsize/30)),'Color','k','LineWidth',3); 
    
    axes(s_Hnd(5));
    jbfill(2:70,squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'r','k',[],0.5);
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
end

for windowsize = 91
    axes(s_Hnd(6)); hold on; xlim([0,70]); ylim([0,1]); grid on
    jbfill(2:70,squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    plot(pcnt_skilful_CPS(:,floor(windowsize/30)),'Color','k','LineWidth',3); 
    
    axes(s_Hnd(6));
    jbfill(2:70,squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'r','k',[],0.5);
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
end
  
for i = [1 4]
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[0:0.2:1]);
    %ylabel(['Proportion'])
end

for i = 1:6
    axes(s_Hnd(i));
    set(gca,'XTickLabel',[0:20:70]);
    line([0 70],[skilful_threshold skilful_threshold],'linestyle','--','color','k','linewidth',0.75)
end

print('ENSO_recon_skill.pdf','-dpdf','-painters','-bestfit')
