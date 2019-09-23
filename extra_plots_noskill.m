%% Skilfullness of different reconstructions - Precip
% First for precip
% Do this for CPS_AAo; AuNZ; SA; SH

clear

numstnstocompare = 2:70; skilful_threshold = sqrt(0.5); % Using a linear statistical model, we want r^2 (Coefficient of determination) to be greater than 0.5, so r must be greater than sqrt(0.5).
NUM_CAL_WDW = 10; clear CAL_WDW; NUM_YRS = 500; NUM_TRIALS = 1000; 

pcnt_skilful_CPS_SH = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_SA = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_AuNZ = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_AAo = nan(max(numstnstocompare),length([31 61 91]));

corr_CPS_qn = nan(70,2);
corr_CPS_qn_SH = nan(70,10,3);
corr_CPS_qn_SA = nan(70,10,3);
corr_CPS_qn_AuNZ = nan(70,10,3);
corr_CPS_qn_AAo = nan(70,10,3);

corr_CPS_qn_rng_SH = nan(3,size(corr_CPS_qn_SH,1),2,size(corr_CPS_qn_SH,3));
corr_CPS_qn_rng_SA = nan(3,size(corr_CPS_qn_SA,1),2,size(corr_CPS_qn_SA,3));
corr_CPS_qn_rng_AuNZ = nan(3,size(corr_CPS_qn_AuNZ,1),2,size(corr_CPS_qn_AuNZ,3));
corr_CPS_qn_rng_AAo = nan(3,size(corr_CPS_qn_AAo,1),2,size(corr_CPS_qn_AAo,3));

precip_corr_CPS_SH = nan(3,70,10,1000);
precip_corr_CPS_SA = nan(3,70,10,1000);
precip_corr_CPS_AuNZ = nan(3,70,10,1000);
precip_corr_CPS_AAo = nan(3,70,10,1000);

windowsize = 500;
DIR_NAME = ['Proxies/NoSkill/',num2str(windowsize),'yrWindow'];
load([DIR_NAME,'/CalWdw:1-500/tonsofstats.mat'], ...
            'all_stn_corr*')
precip_corr_CPS(:,:) = all_stn_corr_CPS;
        
for n=numstnstocompare
    pcnt_skilful_CPS(n) = sum(sum(sum(...
            squeeze(precip_corr_CPS(n,:))>skilful_threshold )))/(1*NUM_TRIALS);
end

for windowsize = [31 61 91]
    DIR_NAME = ['/srv/ccrc/data40/z3215716/SAM/Proxies/NoSkill/',num2str(windowsize),'yrWindow'];
    overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
        clear CAL_WDW;
        for c=0:9
            CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize); 
        end
        
        for c=1:size(CAL_WDW,1)
            load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
            'all_stn_corr*')
            precip_corr_CPS_SH(floor(windowsize/30),:,c,:) = all_stn_corr_CPS; % this is (3,70,10,1000)
			precip_corr_CPS_SA(floor(windowsize/30),:,c,:) = all_stn_corr_CPS_SA;
			precip_corr_CPS_AuNZ(floor(windowsize/30),:,c,:) = all_stn_corr_CPS_AuNz;
			precip_corr_CPS_AAo(floor(windowsize/30),:,c,:) = all_stn_corr_CPS_AAo;
        end

        fprintf('Part 1 successful: %s',num2str(windowsize))
        
        for n = numstnstocompare
            pcnt_skilful_CPS_SH(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS_SH(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
            
            pcnt_skilful_CPS_AAo(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS_AAo(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        
            pcnt_skilful_CPS_SA(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS_SA(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
       
            pcnt_skilful_CPS_AuNZ(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS_AuNZ(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end
        
        fprintf('Part 2 successful: %s',num2str(windowsize))
        
        % Southern Hemisphere recon.
        corr_CPS_qn_SH = quantile(squeeze(precip_corr_CPS_SH(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_SH(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_SH,[],2);
        corr_CPS_qn_rng_SH(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_SH,[],2);
        
        % South America recon.
        corr_CPS_qn_SA = quantile(squeeze(precip_corr_CPS_SA(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_SA(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_SA,[],2);
        corr_CPS_qn_rng_SA(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_SA,[],2);
        
        % Aus/ NZ recon.
        corr_CPS_qn_AuNZ = quantile(squeeze(precip_corr_CPS_AuNZ(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_AuNZ(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_AuNZ,[],2);
        corr_CPS_qn_rng_AuNZ(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_AuNZ,[],2);
        
        % Antarctica recon.
        corr_CPS_qn_AAo = quantile(squeeze(precip_corr_CPS_AAo(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_AAo(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_AAo,[],2);
        corr_CPS_qn_rng_AAo(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_AAo,[],2);
        
        clear corr_CPS_qn_SA corr_CPS_qn_AuNZ all_stn_corr_CPS_AuNz all_stn_corr_CPS_SA % these change size, so need to be cleared after each windowsize is processed.
        fprintf('Part 3 successful: %s',num2str(windowsize))
end
        
        % 500yr recon.
        corr_CPS_qn = quantile(squeeze(precip_corr_CPS(:,:)),[.05 .5 .95], 2);  
        corr_CPS_qn_rng(:,1) = min(corr_CPS_qn,[],2);   
        corr_CPS_qn_rng(:,2) = max(corr_CPS_qn,[],2);
        
        
        
% Plotting
figure(1)
clf
s_Hnd = tight_subplot(3,4,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

for windowsize = 31
    for i = 1:12
        axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
        plot(pcnt_skilful_CPS,'Color',[0 0 0.8],'LineWidth',3); 
        jbfill([2:70],corr_CPS_qn_rng(2:70,2)',corr_CPS_qn_rng(2:70,1)','b','k',[],0.5);
    end
    axes(s_Hnd(1));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(2));
    plot(pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(3));
    plot(pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(4));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(2:70,floor(windowsize/30))-pcnt_skilful_CPS_AAo(2:70,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end

for windowsize = 61
    axes(s_Hnd(5));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(6));
    plot(pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(7));
    plot(pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(8));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end

for windowsize = 91
    axes(s_Hnd(9));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(10));
    plot(pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(11));
    plot(pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(12));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end
  
   
% Formatting
% axes(s_Hnd(1)); title(['South Africa'],'FontSize',16);
% axes(s_Hnd(2)); title(['South America'],'FontSize',16);
% axes(s_Hnd(3)); title(['Aus/ Nz'],'FontSize',16);
% axes(s_Hnd(4)); title(['Antarctica'],'FontSize',16);
% leg_h = legend(s_Hnd(3),[prop_31_Hnd(1,1) prop_61_Hnd(1,1) prop_91_Hnd(1,1)], ...
%                '31yr','61yr','91yr','orientation','horizontal');
% set(leg_h,'FontSize',14,'orientation','horizontal');
% leg_h = legend(s_Hnd(4),[prop_31_Hnd(1,2) prop_61_Hnd(1,2) prop_91_Hnd(1,2)], ...
%                'STAT 31yr','STAT 61yr','STAT 91yr','orientation','horizontal');
% set(leg_h,'FontSize',14,'orientation','horizontal');
% 


for i = [1 5 9]
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[0:0.2:1]);
    %ylabel(['Proportion'])
end

for i = 1:12
    axes(s_Hnd(i));
    set(gca,'XTickLabel',[0:20:70]);
    line([0 70],[skilful_threshold skilful_threshold],'linestyle','--','color','k','linewidth',0.75)
end

h=gcf;

set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

print(gcf, '-dpdf', 'recon_skill_multi_precip_noskill.pdf')


plot2svg('recon_skill_multi_precip_noskill.svg')

%% Skilfullness of different reconstructions - SAT
% First for precip
% Do this for CPS_AAo; AuNZ; SA; SH

clear

numstnstocompare = 2:70; skilful_threshold = sqrt(0.5); % Using a linear statistical model, we want r^2 (Coefficient of determination) to be greater than 0.5, so r must be greater than sqrt(0.5).
NUM_CAL_WDW = 10; clear CAL_WDW; NUM_YRS = 500; NUM_TRIALS = 1000; 

pcnt_skilful_CPS = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that do better than the skill threshold
pcnt_skilful_CPS_SH = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_SA = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_AuNZ = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_AAo = nan(max(numstnstocompare),length([31 61 91]));

corr_CPS_qn = nan(70,3);
corr_CPS_qn_SH = nan(70,10,3);
corr_CPS_qn_SA = nan(70,10,3);
corr_CPS_qn_AuNZ = nan(70,10,3);
corr_CPS_qn_AAo = nan(70,10,3);

corr_CPS_qn_rng = nan(70,2); 
corr_CPS_qn_rng_SH = nan(3,size(corr_CPS_qn_SH,1),2,size(corr_CPS_qn_SH,3)); % creates the range for the 3 windows, 70 sites, 2 bounds of the range, 3 percentiles
corr_CPS_qn_rng_SA = nan(3,size(corr_CPS_qn_SA,1),2,size(corr_CPS_qn_SA,3));
corr_CPS_qn_rng_AuNZ = nan(3,size(corr_CPS_qn_AuNZ,1),2,size(corr_CPS_qn_AuNZ,3));
corr_CPS_qn_rng_AAo = nan(3,size(corr_CPS_qn_AAo,1),2,size(corr_CPS_qn_AAo,3));

sat_corr_CPS_SH = nan(3,70,10,1000);
sat_corr_CPS_SA = nan(3,70,10,1000);
sat_corr_CPS_AuNZ = nan(3,70,10,1000);
sat_corr_CPS_AAo = nan(3,70,10,1000);

windowsize = 500;
DIR_NAME = ['Proxies/NoSkill/',num2str(windowsize),'yrWindow'];
load([DIR_NAME,'/CalWdw:1-500/tonsofstats.mat'], ...
            'all_sat_corr_CPS')
sat_corr_CPS(:,:) = all_sat_corr_CPS;
        
for n=numstnstocompare
    pcnt_skilful_CPS(n) = sum(sum(sum(...
            squeeze(sat_corr_CPS(n,:))>skilful_threshold )))/(1*NUM_TRIALS);
end
for windowsize = [31 61 91]
    DIR_NAME = ['/srv/ccrc/data40/z3215716/SAM/Proxies/NoSkill/',num2str(windowsize),'yrWindow'];
    overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
        clear CAL_WDW;
        for c=0:9
            CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize); %#ok<SAGROW>
        end
        
      
        for c=1:size(CAL_WDW,1)
            load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
            'all_sat_corr*')
            sat_corr_CPS_SH(floor(windowsize/30),:,c,:) = all_sat_corr_CPS; % this is (3,70,10,1000)
			sat_corr_CPS_SA(floor(windowsize/30),:,c,:) = all_sat_corr_CPS_SA;
			sat_corr_CPS_AuNZ(floor(windowsize/30),:,c,:) = all_sat_corr_CPS_AuNz;
			sat_corr_CPS_AAo(floor(windowsize/30),:,c,:) = all_sat_corr_CPS_AAo;
        end

        fprintf('Part 1 successful: %s',num2str(windowsize))
        
        for n = numstnstocompare
            pcnt_skilful_CPS_SH(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(sat_corr_CPS_SH(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
            
            pcnt_skilful_CPS_AAo(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(sat_corr_CPS_AAo(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
       
            pcnt_skilful_CPS_SA(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(sat_corr_CPS_SA(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
      
            pcnt_skilful_CPS_AuNZ(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(sat_corr_CPS_AuNZ(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end
        
        fprintf('Part 2 successful: %s',num2str(windowsize))
        
        % Southern Hemisphere recon.
        corr_CPS_qn_SH = quantile(squeeze(sat_corr_CPS_SH(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_SH(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_SH,[],2);
        corr_CPS_qn_rng_SH(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_SH,[],2);
        
        % South America recon.
        corr_CPS_qn_SA = quantile(squeeze(sat_corr_CPS_SA(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_SA(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_SA,[],2);
        corr_CPS_qn_rng_SA(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_SA,[],2);
        
        % Aus/ NZ recon.
        corr_CPS_qn_AuNZ = quantile(squeeze(sat_corr_CPS_AuNZ(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_AuNZ(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_AuNZ,[],2);
        corr_CPS_qn_rng_AuNZ(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_AuNZ,[],2);
        
        % Antarctica recon.
        corr_CPS_qn_AAo = quantile(squeeze(sat_corr_CPS_AAo(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_AAo(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_AAo,[],2);
        corr_CPS_qn_rng_AAo(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_AAo,[],2);
        
        clear corr_CPS_qn_SA corr_CPS_qn_AuNZ all_sat_corr_CPS_AuNz all_sat_corr_CPS_SA % these change size, so need to be cleared after each windowsize is processed.
        fprintf('Part 3 successful: %s',num2str(windowsize))
end
        

        % 500yr recon.
        corr_CPS_qn = quantile(squeeze(sat_corr_CPS(:,:)),[.05 .5 .95], 2);  
        corr_CPS_qn_rng(:,1) = min(corr_CPS_qn,[],2);   
        corr_CPS_qn_rng(:,2) = max(corr_CPS_qn,[],2);
        
        
% Plotting
figure(2)
clf
s_Hnd = tight_subplot(3,4,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

for windowsize = 31
    for i = 1:4
        axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
        plot(pcnt_skilful_CPS(:,floor(windowsize/30)),'Color',[0 0 0.8],'LineWidth',3); 
        jbfill([2:70],corr_CPS_qn_rng(2:70,2)',corr_CPS_qn_rng(2:70,1)','b','k',[],0.5);
    end
    axes(s_Hnd(1));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(2));
    plot(pcnt_skilful_CPS_SA(2:70,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(3));
    plot(pcnt_skilful_CPS_AuNZ(2:70,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(4));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(2:70,floor(windowsize/30))-pcnt_skilful_CPS_AAo(2:70,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end

for windowsize = 61
    for i = 5:8
        axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
        plot(pcnt_skilful_CPS(:,1),'Color',[0 0 0.8],'LineWidth',3); 
        jbfill([2:70],corr_CPS_qn_rng(2:70,2)',corr_CPS_qn_rng(2:70,1)','b','k',[],0.5);
    end
    axes(s_Hnd(5));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(6));
    plot(pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(7));
    plot(pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(8));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end

for windowsize = 91
    for i = 9:12
        axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
        plot(pcnt_skilful_CPS(:,1),'Color',[0 0 0.8],'LineWidth',3); 
        jbfill([2:70],corr_CPS_qn_rng(2:70,2)',corr_CPS_qn_rng(2:70,1)','b','k',[],0.5);
    end
    axes(s_Hnd(9));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(10));
    plot(pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_SA(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(11));
    plot(pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AuNZ(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(12));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end
  
   
% Formatting
% axes(s_Hnd(1)); title(['South Africa'],'FontSize',16);
% axes(s_Hnd(2)); title(['South America'],'FontSize',16);
% axes(s_Hnd(3)); title(['Aus/ Nz'],'FontSize',16);
% axes(s_Hnd(4)); title(['Antarctica'],'FontSize',16);
% leg_h = legend(s_Hnd(3),[prop_31_Hnd(1,1) prop_61_Hnd(1,1) prop_91_Hnd(1,1)], ...
%                '31yr','61yr','91yr','orientation','horizontal');
% set(leg_h,'FontSize',14,'orientation','horizontal');
% leg_h = legend(s_Hnd(4),[prop_31_Hnd(1,2) prop_61_Hnd(1,2) prop_91_Hnd(1,2)], ...
%                'STAT 31yr','STAT 61yr','STAT 91yr','orientation','horizontal');
% set(leg_h,'FontSize',14,'orientation','horizontal');
% 


for i = [1 5 9]
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[0:0.2:1]);
    %ylabel(['Proportion'])
end

for i = 1:12
    axes(s_Hnd(i));
    set(gca,'XTickLabel',[0:20:70]);
    line([0 70],[skilful_threshold skilful_threshold],'linestyle','--','color','k','linewidth',0.75)
end

%h=gcf;

%set(h,'PaperOrientation','landscape');
%set(h,'PaperUnits','normalized');
%set(h,'PaperPosition', [0 0 1 1]);

print(gcf, '-dpdf', 'recon_skill_multi_SAT_noskill.pdf')

plot2svg('recon_skill_multi_SAT_noskill.svg')

