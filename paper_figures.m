%% 1 Model and Reanalysis SAM/climatologies correlations

clear
load('DataFiles/model_output.mat','lat','lon')
load('DataFiles/SAM_corrs.mat','model_SAT_corr','model_pre_corr')
model_SAT_corr = double(model_SAT_corr); model_pre_corr = double(model_pre_corr);
load('DataFiles/GISS_corrs.mat')
load('DataFiles/GISS_model_output.mat')
load('model_fit.mat')
irange = 144;
jrange = 45;
% land = ncread('DataFiles/sftlf_A1.static.nc','sftlf');
for i = 1:2
    figure(i)
    axesm('MapProjection','pcarree','MapLatLimit',[-90 -0],'maplonlim',[0 360])
    framem;
    gridm;
    mlabel;
    plabel;
    tightmap
end

levels = -1:0.1:1;
levels_n = -1:0.1:-0.3;
levels_p = 0.3:0.1:1;

% model SAT
figure(1)
contourfm(lat(1:45),lon,model_SAT_corr(1:45,:),levels,'linestyle','none')
hold on;
colormap(b2r(-1,1))
for j = 1:jrange
    for i = 1:irange
        if CM2_fit_sat(i,j) ~= 1
            scatterm(lat(j),lon(i),6,CM2_fit_sat(i,j),'k','filled')
        %else
        %    scatterm(lat(j),lon(i),6,CM2_fit_sat(i,j),'g','filled')
        end
    end
end

% model Precip
figure(2)
contourfm(lat(1:45),lon,model_pre_corr(1:45,:),levels,'linestyle','none')
hold on;
colormap(b2r(-1,1))
for i = 1:irange
    for j = 1:jrange
        if CM2_fit_pre(i,j) ~= 1
            scatterm(lat(j),lon(i),6,CM2_fit_pre(i,j),'k','filled')
        %else
        %    scatterm(lat(j),lon(i),6,CM2_fit_pre(i,j),'g','filled')
        end
    end
end

for i = 1:2
    figure(i)
    load coast
    plotm(lat,long,'k','linewidth',2)
end

figure(1)
print('CM2p1_sat_fit.pdf','-dpdf','-painters')%,'-bestfit')
figure(2)
print('CM2p1_precip_fit.pdf','-dpdf','-painters')%,'-bestfit')
% What % of sites fit within the model range?
% Total # of SH land cells  = 2241 
num_land = sum(sum(~isnan(land(:,1:45))));

%% Figure 2 - Apparent vs. true correlation percentiles

% Setup
clear
load DataFiles/model_output.mat

windowsize = 31; % The running window in years. Change for each set of figures

corr_sat = zeros(45,size(sat_detr,3));
P_sat = zeros(45,size(sat_detr,3));
corr_precip = zeros(45,size(precip_detr,3));
P_precip = zeros(45,size(precip_detr,3));

for i=1:45
    for j=1:size(sat_detr,3)
        corr_sat(i,j) = corr(SAM',sat_detr(:,i,j));
    end
end

for i=1:45
    for j=1:size(precip_detr,3)
        corr_precip(i,j) = corr(SAM',precip_detr(:,i,j));
    end
end

load(['DataFiles/runcorr',num2str(windowsize),'yrwdw.mat']);
% Temp
bin = -1.0:0.03:1.0;
[sorted_corr, sorted_corr_ind] = sort(corr_sat(:));
sat_runcorr = sat_runcorr(:,1:45,:);
sorted_sat_runcorr = sat_runcorr_31(:,sorted_corr_ind); 
bin_sizes = histc(squeeze(sorted_corr),bin);

current_index = 1; sat_runcorr_quan = nan(7,length(bin_sizes));
for m=1:length(bin_sizes)
    sat_runcorr_quan(:,m) = quantile(reshape(sorted_sat_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end

% Precip

[sorted_corr_p, sorted_corr_ind_p] = sort(corr_precip(:));
precip_runcorr = precip_runcorr(:,1:45,:);
sorted_precip_runcorr = precip_runcorr(:,sorted_corr_ind_p);
bin_sizes_p = histc(squeeze(sorted_corr_p),bin);

current_index = 1; precip_runcorr_quan = nan(7,length(bin_sizes_p));
for m=1:length(bin_sizes_p)
    precip_runcorr_quan(:,m) = quantile(reshape(sorted_precip_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes_p(m);
end

% Plotting

figure(1)
hold on; Hnd = nan(1,7);
for n=1:7
    HA(n) = plot(bin,sat_runcorr_quan(n,:));  % repeat for other two windows
    set(HA(n),'Color','k','LineWidth',3);
end
plot([-0.3,-0.3],[-1 1],'k','LineWidth',2); plot([0.3,0.3],[-1 1],'k','LineWidth',2)
plot([-1 1],[-0.3,-0.3],'k','LineWidth',2); plot([-1 1],[0.3,0.3],'k','LineWidth',2)
hold off; grid on; axis equal; axis([-1 1 -1 1]);
set(HA([3,5]),'Visible','off')
set(HA([1,7]),'LineStyle','-','LineWidth',2);
set(HA([4]),'LineStyle','-','LineWidth',4);
set(text(0,0.30,'0.95'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.1,0.25,'0.75'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,0,'0.50'),'FontSize',15,'Rotation',40,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.15,0.07,'0.25'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,-0.28,'0.05'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.26,0.76,'0.99'),'FontSize',15,'Rotation',35,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.31,-0.33,'0.01'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
ylabel(['r(',num2str(windowsize),'yr)'],'FontSize',16); % Change windowsize when doing each figure
xlabel('r(500yr)','FontSize',16);
title(['Percentiles'],'FontSize',16,'FontWeight','bold')
set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.05 .01] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , -1:0.5:1, ...
    'XTick'       , -1:0.5:1, ...
    'LineWidth'   , 2,  ...    
    'YLim'        , [-1 1], ...
    'XLim'        , [-1 1], ...
    'FontSize'    , 16   ...
    );
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [19 19]);
set(gcf, 'PaperPosition', [0 0 19 19]); %x_width=19cm y_width=28cm
text(-1+0.08,1-0.18,'(a)','FontSize',22,'FontWeight','bold');

figure(2)

hold on; Hnd = nan(1,7);
for n=1:7
    HA(n) = plot(bin,precip_runcorr_quan(n,:)); % Repeat for other two windows
    set(HA(n),'Color','k','LineWidth',3);
end
plot([-0.3,-0.3],[-1 1],'k','LineWidth',2); plot([0.3,0.3],[-1 1],'k','LineWidth',2)
plot([-1 1],[-0.3,-0.3],'k','LineWidth',2); plot([-1 1],[0.3,0.3],'k','LineWidth',2)
hold off; grid on; axis equal; axis([-1 1 -1 1]);
set(HA([3,5]),'Visible','off')
set(HA([1,7]),'LineStyle','-','LineWidth',2);
set(HA([4]),'LineStyle','-','LineWidth',4);
set(text(0,0.30,'0.95'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.1,0.25,'0.75'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,0,'0.50'),'FontSize',15,'Rotation',40,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.15,0.07,'0.25'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center','Visible','off')
set(text(0,-0.28,'0.05'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.26,0.76,'0.99'),'FontSize',15,'Rotation',35,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
set(text(0.31,-0.33,'0.01'),'FontSize',15,'Rotation',45,'Margin',0.1,'BackgroundColor',[1 1 1],'HorizontalAlignment','center')
ylabel(['r(',num2str(windowsize),'yr)'],'FontSize',16);
xlabel('r(500yr)','FontSize',16);
title(['Percentiles'],'FontSize',16,'FontWeight','bold')
set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.05 .01] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , -1:0.5:1, ...
    'XTick'       , -1:0.5:1, ...
    'LineWidth'   , 2,  ...    
    'YLim'        , [-1 1], ...
    'XLim'        , [-1 1], ...
    'FontSize'    , 16   ...
    );
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [19 19]);
set(gcf, 'PaperPosition', [0 0 19 19]); %x_width=19cm y_width=28cm
text(-1+0.08,1-0.18,'(b)','FontSize',22,'FontWeight','bold');


%% Figure 3 - SH wide SAT recons with and without AA

clear

NUM_YRS = 500;
numstnstocompare = 2:70;
NUM_TRIALS = 1000;
load site_range.mat

for region = [1 4] % Paper only requires these two
	for windowsize = [31, 61, 91] % 500
      DIR_NAME = ['Proxies/NoResample/',num2str(windowsize),'yrWindow/'];
      NUM_CAL_WDW = 10; clear CAL_WDW;
      overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
      for c=0:9  % Change to 0 for 500yr window
        CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
      end
    
      SAT_corr_CPS = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
    
		for c=1:size(CAL_WDW,1)
			load([DIR_NAME,'/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
			'all_sat_corr_CPS','all_sat_corr_CPS_SA','all_sat_corr_CPS_AuNz','all_sat_corr_CPS_AA','all_sat_corr_CPS_AAo','all_sat_corr_CPS_SoA')
            if region == 1
				SAT_corr_CPS(:,c,:) = all_sat_corr_CPS;     % Whole SH 
			elseif region == 2
				SAT_corr_CPS(:,c,:) = all_sat_corr_CPS_SA;  % S. America only
			elseif region == 3
				SAT_corr_CPS(:,c,:) = all_sat_corr_CPS_AuNz;% Aus/NZ only
			elseif region == 4
				SAT_corr_CPS(:,c,:) = all_sat_corr_CPS_AA;  % SH minus Antarctica
			elseif region == 5
				SAT_corr_CPS(:,c,:) = all_sat_corr_CPS_AAo; % Antarctica only
			elseif region == 6
				SAT_corr_CPS(:,c,:) = all_sat_corr_CPS_SoA; % S. Africa only
			end
		end

% Plotting CPS
      if region == 1
        figure(1)
      elseif region == 2
        figure(2)
      elseif region == 3
        figure(3)
	  elseif region == 4
        figure(4)
	  elseif region == 5
        figure(5)
	  elseif region == 6
        figure(6)
	  end
      if windowsize == 31
        subplot(1,3,1)
      elseif windowsize == 61
        subplot(1,3,2)
      elseif windowsize == 91
        subplot(1,3,3)
      end
      
      corr_RV_qn = quantile(SAT_corr_CPS,[.05 .5 .95], 3);
      corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
      corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
      corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
      
      if windowsize == 500
         range = 2:70;
      elseif windowsize < 500
      range = 2:sat_min(floor(windowsize/30),region);
        if max(range) > 70
         range = 2:70;
        end
      end
      
      jbfill(range,squeeze(corr_RV_qn_rng(range,2,1))',squeeze(corr_RV_qn_rng(range,1,1))','b','k',[],0.5);
      jbfill(range,squeeze(corr_RV_qn_rng(range,2,3))',squeeze(corr_RV_qn_rng(range,1,3))','r','k','add',0.5);
      jbfill(range,squeeze(corr_RV_qn_rng(range,2,2))',squeeze(corr_RV_qn_rng(range,1,2))','y','k','add',0.5);
       
      xlim([0,70]); ylim([0,1]); grid on
      
      ylabel(['r(',num2str(windowsize),'yr Window)']);
      %title(['CPS'])
      set(gca, 'FontSize',18, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
      set(gca,'XTickLabel',0:20:70);
      if windowsize == 61  
        legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','southeast');
        set(legendH, 'FontSize',10);
        xlabel('Number of proxies in reconstruction');
      end
      %set(gcf, 'PaperPosition', [0 0 19 28]);
      set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on','YTick', [0:0.1:1]);
   end
end
figure(1)
h = gcf;
set(h,'Units','Inches');
%pos = get(h,'Position'); Below is the value we choose
pos = [2.20833333333333,0.697916666666667,12.1041666666667,9.38541666666667];
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'Figures/SAM_reconskill_glob_SAT','-dpdf','-bestfit')

figure(4)
h = gcf;
set(h,'Units','Inches');
pos = [2.20833333333333,0.697916666666667,12.1041666666667,9.38541666666667];
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'Figures/SAM_reconskill_noAA_SAT','-dpdf','-bestfit')

%% Figure 3 - SH wide precip recons with and without AA

clear

NUM_YRS = 500;
numstnstocompare = 2:70;
NUM_TRIALS = 1000;
load site_range.mat

for region = [1 4] % Paper uses only these two
	for windowsize = [31, 61, 91] % 500
      DIR_NAME = ['Proxies/NoResample/',num2str(windowsize),'yrWindow/'];
      NUM_CAL_WDW = 10; clear CAL_WDW;
      overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
      for c=0:9  % Change to 0 for 500yr window
        CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
      end
    
      precip_corr_CPS = nan(max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');
    
		for c=1:size(CAL_WDW,1)
			load([DIR_NAME,'/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
			'all_stn_corr_CPS','all_stn_corr_CPS_SA','all_stn_corr_CPS_AuNz','all_stn_corr_CPS_AA','all_stn_corr_CPS_AAo','all_stn_corr_CPS_SoA')
            if region == 1
				precip_corr_CPS(:,c,:) = all_stn_corr_CPS;     % Whole SH
			elseif region == 2
				precip_corr_CPS(:,c,:) = all_stn_corr_CPS_SA;  % S. America only
			elseif region == 3
				precip_corr_CPS(:,c,:) = all_stn_corr_CPS_AuNz;% Aus/NZ only
			elseif region == 4
				precip_corr_CPS(:,c,:) = all_stn_corr_CPS_AA;  % SH minus Antarctica
			elseif region == 5
				precip_corr_CPS(:,c,:) = all_stn_corr_CPS_AAo; % Antarctica only
			elseif region == 6
				precip_corr_CPS(:,c,:) = all_stn_corr_CPS_SoA; % S. Africa only
			end
		end

% Plotting CPS
      if region == 1
        figure(1)
      elseif region == 2
        figure(2)
      elseif region == 3
        figure(3)
	  elseif region == 4
        figure(4)
	  elseif region == 5
        figure(5)
	  elseif region == 6
        figure(6)
	  end
      if windowsize == 31
        subplot(1,3,1)
      elseif windowsize == 61
        subplot(1,3,2)
      elseif windowsize == 91
        subplot(1,3,3)
      end
      
      corr_RV_qn = quantile(precip_corr_CPS,[.05 .5 .95], 3);
      corr_RV_qn_rng = nan(size(corr_RV_qn,1),2,size(corr_RV_qn,3));
      corr_RV_qn_rng(:,1,:) = min(corr_RV_qn,[],2);
      corr_RV_qn_rng(:,2,:) = max(corr_RV_qn,[],2);
      
      if windowsize == 500
         range = 2:70;
      elseif windowsize < 500
      range = 2:precip_min(floor(windowsize/30),region);
        if max(range) > 70
         range = 2:70;
        end
      end
      
      jbfill(range,squeeze(corr_RV_qn_rng(range,2,1))',squeeze(corr_RV_qn_rng(range,1,1))','b','k',[],0.5);
      jbfill(range,squeeze(corr_RV_qn_rng(range,2,3))',squeeze(corr_RV_qn_rng(range,1,3))','r','k','add',0.5);
      jbfill(range,squeeze(corr_RV_qn_rng(range,2,2))',squeeze(corr_RV_qn_rng(range,1,2))','y','k','add',0.5);
       
      xlim([0,70]); ylim([0,1]); grid on
      
      ylabel(['r(',num2str(windowsize),'yr Window)']);
      %title(['CPS'])
      set(gca, 'FontSize',18, 'LineWidth', 1.0, 'Box', 'on', 'YTick', [0:0.1:1]); 
      set(gca,'XTickLabel',[0:20:70]);
      if windowsize == 61
        legendH = legend('5^t^h Percentile Range','95^t^h Percentile Range','Median Range','location','southeast');
        set(legendH, 'FontSize',10);
        xlabel('Number of proxies in reconstruction');
      end
      set(gcf, 'PaperPosition', [0 0 19 28]);
      set(gca, 'FontSize',14, 'LineWidth', 1.0, 'Box', 'on','YTick', [0:0.1:1]);
   end
end

figure(1)
h = gcf;
set(h,'Units','Inches');
%pos = get(h,'Position'); Below is the value we choose
pos = [2.20833333333333,0.697916666666667,12.1041666666667,9.38541666666667];
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'Figures/SAM_reconskill_glob_precip','-dpdf','-bestfit')

figure(4)
h = gcf;
set(h,'Units','Inches');
pos = [2.20833333333333,0.697916666666667,12.1041666666667,9.38541666666667];
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'Figures/SAM_reconskill_noAA_precip','-dpdf','-bestfit')

%% Figure 4 - SH and regional reconstruction skill for SAT

clear
load('site_range.mat','sat_min')
numstnstocompare = 2:70; skilful_threshold = sqrt(0.5); % Using a linear statistical model, we want r^2 (Coefficient of determination) to be greater than 0.5, so r must be greater than sqrt(0.5).
NUM_CAL_WDW = 10; clear CAL_WDW; NUM_YRS = 500; NUM_TRIALS = 1000; 

pcnt_skilful_CPS = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that do better than the skill threshold
pcnt_skilful_CPS_SH = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_SA = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_AuNZ = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_AAo = nan(max(numstnstocompare),length([31 61 91]));

corr_CPS_qn = nan(70,10,3);
corr_CPS_qn_SH = nan(70,10,3);
corr_CPS_qn_SA = nan(70,10,3);
corr_CPS_qn_AuNZ = nan(70,10,3);
corr_CPS_qn_AAo = nan(70,10,3);

corr_CPS_qn_rng = nan(size(corr_CPS_qn,1),2); % Creates the range for 70 sites, max and min (only one window)
corr_CPS_qn_rng_SH = nan(3,size(corr_CPS_qn_SH,1),2,size(corr_CPS_qn_SH,3)); % creates the range for the 3 windows, 70 sites, 2 bounds of the range, 3 percentiles
corr_CPS_qn_rng_SA = nan(3,size(corr_CPS_qn_SA,1),2,size(corr_CPS_qn_SA,3));
corr_CPS_qn_rng_AuNZ = nan(3,size(corr_CPS_qn_AuNZ,1),2,size(corr_CPS_qn_AuNZ,3));
corr_CPS_qn_rng_AAo = nan(3,size(corr_CPS_qn_AAo,1),2,size(corr_CPS_qn_AAo,3));

SAT_corr_CPS_SH = nan(3,70,10,1000);
SAT_corr_CPS_SA = nan(3,70,10,1000);
SAT_corr_CPS_AuNZ = nan(3,70,10,1000);
SAT_corr_CPS_AAo = nan(3,70,10,1000);

windowsize = 500;
DIR_NAME = ['/media/huiskamp/My Book/CM2.1/SAM/Proxies/NoResample/',num2str(windowsize),'yrWindow'];
load([DIR_NAME,'/CalWdw1_500/tonsofstats.mat'], ...
            'all_sat_corr*')
SAT_corr_CPS(:,:) = all_sat_corr_CPS;
        
for n=numstnstocompare
    pcnt_skilful_CPS(n) = sum(sum(sum(...
            squeeze(SAT_corr_CPS(n,:))>skilful_threshold )))/(1*NUM_TRIALS);
end
for windowsize = [31 61 91]
    DIR_NAME = ['/media/huiskamp/My Book/CM2.1/SAM/Proxies/NoResample/',num2str(windowsize),'yrWindow'];
    overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
        clear CAL_WDW;
        for c=0:9
            CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
        end
        
        SA = 1:sat_min(floor(windowsize/30),2);
        AUNZ = 1:sat_min(floor(windowsize/30),3);
        for c=1:size(CAL_WDW,1)
            load([DIR_NAME,'/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
            'all_sat_corr*')
            SAT_corr_CPS_SH(floor(windowsize/30),:,c,:) = all_sat_corr_CPS; % this is (3,70,10,1000)
			SAT_corr_CPS_SA(floor(windowsize/30),SA,c,:) = all_sat_corr_CPS_SA;
			SAT_corr_CPS_AuNZ(floor(windowsize/30),AUNZ,c,:) = all_sat_corr_CPS_AuNz;
			SAT_corr_CPS_AAo(floor(windowsize/30),:,c,:) = all_sat_corr_CPS_AAo;
        end

        fprintf('Part 1 successful: %s',num2str(windowsize))
        
        for n = numstnstocompare
            pcnt_skilful_CPS_SH(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(SAT_corr_CPS_SH(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
            
            pcnt_skilful_CPS_AAo(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(SAT_corr_CPS_AAo(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end
        for n = 2:sat_min(floor(windowsize/30),2)
            pcnt_skilful_CPS_SA(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(SAT_corr_CPS_SA(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end
        for n = 2:sat_min(floor(windowsize/30),3)
            pcnt_skilful_CPS_AuNZ(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(SAT_corr_CPS_AuNZ(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end
        
        fprintf('Part 2 successful: %s',num2str(windowsize))
        
        % Southern Hemisphere recon.
        corr_CPS_qn_SH = quantile(squeeze(SAT_corr_CPS_SH(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_SH(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_SH,[],2);
        corr_CPS_qn_rng_SH(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_SH,[],2);
        
        % South America recon.
        corr_CPS_qn_SA = quantile(squeeze(SAT_corr_CPS_SA(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_SA(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_SA,[],2);
        corr_CPS_qn_rng_SA(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_SA,[],2);
        
        % Aus/ NZ recon.
        corr_CPS_qn_AuNZ = quantile(squeeze(SAT_corr_CPS_AuNZ(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_AuNZ(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_AuNZ,[],2);
        corr_CPS_qn_rng_AuNZ(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_AuNZ,[],2);
        
        % Antarctica recon.
        corr_CPS_qn_AAo = quantile(squeeze(SAT_corr_CPS_AAo(floor(windowsize/30),:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng_AAo(floor(windowsize/30),:,1,:) = min(corr_CPS_qn_AAo,[],2);
        corr_CPS_qn_rng_AAo(floor(windowsize/30),:,2,:) = max(corr_CPS_qn_AAo,[],2);
        
        clear corr_CPS_qn_SA corr_CPS_qn_AuNZ all_sat_corr_CPS_AuNz all_sat_corr_CPS_SA % these change size, so need to be cleared after each windowsize is processed.
        fprintf('Part 3 successful: %s',num2str(windowsize))
end
        
% 500yr recon.
corr_CPS_qn = quantile(squeeze(SAT_corr_CPS(:,:)),[.05 .5 .95], 2);  % This creates a (70,3) file ie: 70 sites and 3 percentiles
corr_CPS_qn_rng(:,1) = min(corr_CPS_qn,[],2); % This creates a (3,70,2,3) variable: 3 windows, 70 sites, max and min, 3 percentiles
corr_CPS_qn_rng(:,2) = max(corr_CPS_qn,[],2);


% Plotting
figure(2)
clf
s_Hnd = tight_subplot(3,4,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

for windowsize = 31
    range_SA = 2:sat_min(floor(windowsize/30),2);
    range_AuNz = 2:sat_min(floor(windowsize/30),3);
    for i = 1:4
        axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
        plot(pcnt_skilful_CPS(:,1),'Color',[0 0 0.8],'LineWidth',3); 
        jbfill(2:70,corr_CPS_qn_rng(2:70,2)',corr_CPS_qn_rng(2:70,1)','b','k',[],0.5);
    end
    axes(s_Hnd(1));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(2));
    plot(pcnt_skilful_CPS_SA(range_SA,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_SA,squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,1,1)),'y','k',[],0.5);
    axes(s_Hnd(3));
    plot(pcnt_skilful_CPS_AuNZ(range_AuNz,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_AuNz,squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,1,1)),'y','k',[],0.5);
    axes(s_Hnd(4));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end

for windowsize = 61
    range_SA = 2:sat_min(floor(windowsize/30),2);
    range_AuNz = 2:sat_min(floor(windowsize/30),3);
    for i = 5:8
        axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
        plot(pcnt_skilful_CPS(:,1),'Color',[0 0 0.8],'LineWidth',3); 
        jbfill(2:70,corr_CPS_qn_rng(2:70,2)',corr_CPS_qn_rng(2:70,1)','b','k',[],0.5);
    end
    axes(s_Hnd(5));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(6));
    plot(pcnt_skilful_CPS_SA(range_SA,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_SA,squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,1,1)),'y','k',[],0.5);
    axes(s_Hnd(7));
    plot(pcnt_skilful_CPS_AuNZ(range_AuNz,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_AuNz,squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,1,1)),'y','k',[],0.5);
    axes(s_Hnd(8));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end

for windowsize = 91
    range_SA = 2:sat_min(floor(windowsize/30),2);
    range_AuNz = 2:sat_min(floor(windowsize/30),3);
    for i = 9:12
        axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
        plot(pcnt_skilful_CPS(:,1),'Color',[0 0 0.8],'LineWidth',3); 
        jbfill(2:70,corr_CPS_qn_rng(2:70,2)',corr_CPS_qn_rng(2:70,1)','b','k',[],0.5);
    end
    axes(s_Hnd(9));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(10));
    plot(pcnt_skilful_CPS_SA(range_SA,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_SA,squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,1,1)),'y','k',[],0.5);
    axes(s_Hnd(11));
    plot(pcnt_skilful_CPS_AuNZ(range_AuNz,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_AuNz,squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,1,1)),'y','k',[],0.5);
    axes(s_Hnd(12));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end
  
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

%% Figure 5 - SH and regional reconstruction skill for precip

clear
load('site_range.mat','precip_min')
numstnstocompare = 2:70; skilful_threshold = sqrt(0.5); % Using a linear statistical model, we want r^2 (Coefficient of determination) to be greater than 0.5, so r must be greater than sqrt(0.5).
NUM_CAL_WDW = 10; clear CAL_WDW; NUM_YRS = 500; NUM_TRIALS = 1000; 

pcnt_skilful_CPS = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that do better than the skill threshold
pcnt_skilful_CPS_SH = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_SA = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_AuNZ = nan(max(numstnstocompare),length([31 61 91]));
pcnt_skilful_CPS_AAo = nan(max(numstnstocompare),length([31 61 91]));

corr_CPS_qn = nan(70,10,3);
corr_CPS_qn_SH = nan(70,10,3);
corr_CPS_qn_SA = nan(70,10,3);
corr_CPS_qn_AuNZ = nan(70,10,3);
corr_CPS_qn_AAo = nan(70,10,3);

corr_CPS_qn_rng = nan(3,size(corr_CPS_qn,1),2,size(corr_CPS_qn,3)); % creates the range for the 3 windows, 70 sites, 2 bounds of the range, 3 percentiles
corr_CPS_qn_rng_SH = nan(3,size(corr_CPS_qn_SH,1),2,size(corr_CPS_qn_SH,3));
corr_CPS_qn_rng_SA = nan(3,size(corr_CPS_qn_SA,1),2,size(corr_CPS_qn_SA,3));
corr_CPS_qn_rng_AuNZ = nan(3,size(corr_CPS_qn_AuNZ,1),2,size(corr_CPS_qn_AuNZ,3));
corr_CPS_qn_rng_AAo = nan(3,size(corr_CPS_qn_AAo,1),2,size(corr_CPS_qn_AAo,3));

precip_corr_CPS_SH = nan(3,70,10,1000);
precip_corr_CPS_SA = nan(3,70,10,1000);
precip_corr_CPS_AuNZ = nan(3,70,10,1000);
precip_corr_CPS_AAo = nan(3,70,10,1000);

windowsize = 500;
DIR_NAME = ['/media/huiskamp/My Book/CM2.1/SAM/Proxies/NoResample/',num2str(windowsize),'yrWindow'];
load([DIR_NAME,'/CalWdw1_500/tonsofstats.mat'], ...
            'all_stn_corr*')
precip_corr_CPS(:,:) = all_stn_corr_CPS;
        
for n=numstnstocompare
    pcnt_skilful_CPS(n) = sum(sum(sum(...
            squeeze(precip_corr_CPS(n,:))>skilful_threshold )))/(1*NUM_TRIALS);
end
for windowsize = [31 61 91]
    DIR_NAME = ['/media/huiskamp/My Book/CM2.1/SAM/Proxies/NoResample/',num2str(windowsize),'yrWindow'];
    overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
        clear CAL_WDW;
        for c=0:9
            CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize); %#ok<SAGROW>
        end
        
        SA = 1:precip_min(floor(windowsize/30),2);
        AUNZ = 1:precip_min(floor(windowsize/30),3);
        for c=1:size(CAL_WDW,1)
            load([DIR_NAME,'/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
            'all_stn_corr*')
            precip_corr_CPS_SH(floor(windowsize/30),:,c,:) = all_stn_corr_CPS; % this is (3,70,10,1000)
			precip_corr_CPS_SA(floor(windowsize/30),SA,c,:) = all_stn_corr_CPS_SA;
			precip_corr_CPS_AuNZ(floor(windowsize/30),AUNZ,c,:) = all_stn_corr_CPS_AuNz;
			precip_corr_CPS_AAo(floor(windowsize/30),:,c,:) = all_stn_corr_CPS_AAo;
        end

        fprintf('Part 1 successful: %s',num2str(windowsize))
        
        for n = numstnstocompare
            pcnt_skilful_CPS_SH(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS_SH(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
            
            pcnt_skilful_CPS_AAo(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS_AAo(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end
        for n = 2:precip_min(floor(windowsize/30),2)
            pcnt_skilful_CPS_SA(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS_SA(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end
        for n = 2:precip_min(floor(windowsize/30),3)
            pcnt_skilful_CPS_AuNZ(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS_AuNZ(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end
        
        fprintf('Part 2 successful: %s',num2str(windowsize))
        
        % 500yr recon.
        corr_CPS_qn = quantile(squeeze(precip_corr_CPS(:,:,:)),[.05 .5 .95], 3);  % This creates a (70,10,3) file ie: 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng(floor(windowsize/30),:,1,:) = min(corr_CPS_qn,[],2);      % This creates a (3,70,2,3) variable: 3 windows, 70 sites, max and min, 3 percentiles
        corr_CPS_qn_rng(floor(windowsize/30),:,2,:) = max(corr_CPS_qn,[],2);
        
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
        
% Plotting
figure(1)
clf
s_Hnd = tight_subplot(3,4,[0.05 0.01],[0.10 0.01],[0.1 0.01]);

for windowsize = 31
    range_SA = 2:precip_min(floor(windowsize/30),2);
    range_AuNz = 2:precip_min(floor(windowsize/30),3);
    for i = 1:4
        axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
        plot(pcnt_skilful_CPS(:,1),'Color',[0 0 0.8],'LineWidth',3); 
        jbfill([2:70],squeeze(corr_CPS_qn_rng(1,2:70,2,3)),squeeze(corr_CPS_qn_rng(1,2:70,1,1)),'b','k',[],0.5);
    end
    axes(s_Hnd(1));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(2));
    plot(pcnt_skilful_CPS_SA(range_SA,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_SA,squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,1,1)),'y','k',[],0.5);
    axes(s_Hnd(3));
    plot(pcnt_skilful_CPS_AuNZ(range_AuNz,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_AuNz,squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,1,1)),'y','k',[],0.5);
    axes(s_Hnd(4));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end

for windowsize = 61
    range_SA = 2:precip_min(floor(windowsize/30),2);
    range_AuNz = 2:precip_min(floor(windowsize/30),3);
    for i = 5:8
        axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
        plot(pcnt_skilful_CPS(:,1),'Color',[0 0 0.8],'LineWidth',3); 
        jbfill([2:70],squeeze(corr_CPS_qn_rng(1,2:70,2,3)),squeeze(corr_CPS_qn_rng(1,2:70,1,1)),'b','k',[],0.5);
    end
    axes(s_Hnd(5));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(6));
    plot(pcnt_skilful_CPS_SA(range_SA,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_SA,squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,1,1)),'y','k',[],0.5);
    axes(s_Hnd(7));
    plot(pcnt_skilful_CPS_AuNZ(range_AuNz,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_AuNz,squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,1,1)),'y','k',[],0.5);
    axes(s_Hnd(8));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    %plot(pcnt_skilful_CPS(:,floor(windowsize/30))-pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','k','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end

for windowsize = 91
    range_SA = 2:precip_min(floor(windowsize/30),2);
    range_AuNz = 2:precip_min(floor(windowsize/30),3);
    for i = 9:12
        axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
        plot(pcnt_skilful_CPS(:,1),'Color',[0 0 0.8],'LineWidth',3); 
        jbfill([2:70],squeeze(corr_CPS_qn_rng(1,2:70,2,3)),squeeze(corr_CPS_qn_rng(1,2:70,1,1)),'b','k',[],0.5);
    end
    axes(s_Hnd(9));
    plot(pcnt_skilful_CPS_SH(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_SH(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
    axes(s_Hnd(10));
    plot(pcnt_skilful_CPS_SA(range_SA,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_SA,squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,2,3)),squeeze(corr_CPS_qn_rng_SA(floor(windowsize/30),range_SA,1,1)),'y','k',[],0.5);
    axes(s_Hnd(11));
    plot(pcnt_skilful_CPS_AuNZ(range_AuNz,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill(range_AuNz,squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,2,3)),squeeze(corr_CPS_qn_rng_AuNZ(floor(windowsize/30),range_AuNz,1,1)),'y','k',[],0.5);
    axes(s_Hnd(12));
    plot(pcnt_skilful_CPS_AAo(:,floor(windowsize/30)),'Color','red','LineWidth',3);
    jbfill([2:70],squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,2,3)),squeeze(corr_CPS_qn_rng_AAo(floor(windowsize/30),2:70,1,1)),'y','k',[],0.5);
end
  
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

%% Figure 6  Nonstationary cells 

clear
wdwsize = 61;
load(['DataFiles/nonstat_map_',num2str(wdwsize),'yrwdw.mat']) % Change for each window size
load('DataFiles/model_output.mat')
levels = 10;
levels_sig = [47 60 70 80 90 100; 44 60 70 80 90 100; 41 60 70 80 90 100];

for i = 1:2
    figure(i)
    axesm('MapProjection','stereo','origin',[-90,0],'MapLatLimit',[-90 0])
    framem
    gridm
    %mlabel  % commented out - we add our own labels later
    tightmap
end

figure(1)
contourfm(lat,lon,nonstat_precipmap,'linestyle','none');
hold on
contourm(lat,lon,nonstat_precipmap,levels_sig(floor(wdwsize/30),:),'color','k')
colorbar;
caxis([0, 100]);
colormap(flipud(hot(levels)));

figure(2)
contourfm(lat,lon,nonstat_satmap,'linestyle','none');
hold on
contourm(lat,lon,nonstat_satmap,levels_sig(floor(wdwsize/30),:),'color','k')
colorbar;
caxis([0, 100]);
colormap(flipud(hot(levels)));

load coast
for i=1:2
    figure(i)
    plotm(lat,long,'k','linewidth',2)
end

figure(1);print(['nonstat_precipmap',num2str(wdwsize),'yrwdw.pdf'],'-dpdf','-painters')
figure(2);print(['nonstat_SATmap',num2str(wdwsize),'yrwdw.pdf'],'-dpdf','-painters')


%% Figure 7 Impact of non-stationary sites

clear all
NUM_YRS = 500;
num_prox = 70;
NUM_CAL_WDW = 10;
nstat_corrs = nan(3,2,69); % windows, vars, network size
regr_slopes = nan(3,2,69);
for windowsize = [31 61 91]
    for group_size = 2:num_prox
        clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
        for c = 0:9
            CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize);
        end
        for c = 1:size(CAL_WDW,1)
            load(['Proxies/NoResample/',num2str(windowsize),'yrWindow/CalWdw',num2str(CAL_WDW(c,1)),'_',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'],...
                'all_sat_corr_CPS','all_stn_corr_CPS')
            sat_rvals(c,:) = all_sat_corr_CPS(group_size,:);
            stn_rvals(c,:) = all_stn_corr_CPS(group_size,:);
        end
        load(['DataFiles/num_nstat',num2str(windowsize),'yrWdw_',num2str(group_size),'prox.mat'])
        sat_r = reshape(sat_rvals,1,10*1000);
        pre_r = reshape(stn_rvals,1,10*1000);
        nstat_sat = reshape(nstat_nstns_sat,1,10*1000)./group_size;
        nstat_pre = reshape(nstat_nstns_pre,1,10*1000)./group_size;
        [Rs1 Ps1] = corrcoef(sat_r,nstat_sat);
        [Rs2 Ps2] = corrcoef(pre_r,nstat_pre);
        nstat_corrs(floor(windowsize/30),1,group_size) = Rs1(1,2); nstat_Ps(floor(windowsize/30),1,group_size) = Ps1(1,2);
        nstat_corrs(floor(windowsize/30),2,group_size) = Rs2(1,2); nstat_Ps(floor(windowsize/30),2,group_size) = Ps2(1,2);
        p_sat = polyfit(sat_r,nstat_sat,1); p_pre = polyfit(pre_r,nstat_pre,1);
        regr_slopes(floor(windowsize/30),1,group_size) = p_sat(1);
        regr_slopes(floor(windowsize/30),2,group_size) = p_pre(1);
    end
end
        
figure(1)
subplot(1,2,1)
hold on 
grid on
% SAT
plot(squeeze(nstat_corrs(1,1,:)),'o','color',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(nstat_corrs(2,1,:)),'o','color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(nstat_corrs(3,1,:)),'o','color',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'markersize',3,'linestyle','-','linewidth',2)     
% Precip
plot(squeeze(nstat_corrs(1,2,:)),'o','color',[0.4940, 0.1840, 0.5560],'MarkerFaceColor',[0.4940, 0.1840, 0.5560],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(nstat_corrs(2,2,:)),'o','color',[0.4660, 0.6740, 0.1880],'MarkerFaceColor',[0.4660, 0.6740, 0.1880],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(nstat_corrs(3,2,:)),'o','color',[0.3010, 0.7450, 0.9330],'MarkerFaceColor',[0.3010, 0.7450, 0.9330],'markersize',3,'linestyle','-','linewidth',2)
line([0 70],[0 0],'linestyle','--','color','k') % 0 line
line([0 70],[-0.02 -0.02],'color','r'); line([0 70],[0.027 0.027],'color','r') % region where corrs are not significant at p < 0.05
legend('SAT-31yrWdw','SAT-61yrWdw','SAT-91yrWdw','Pre-31yrWdw','Pre-61yrWdw','Pre-91yrWdw')
xlabel('Network Size'); ylabel('r');

% Plot regression slopes...
subplot(1,2,2)
hold on 
grid on
% SAT
plot(squeeze(regr_slopes(1,1,:)),'o','color',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(regr_slopes(2,1,:)),'o','color',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(regr_slopes(3,1,:)),'o','color',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'markersize',3,'linestyle','-','linewidth',2)     
% Precip
plot(squeeze(regr_slopes(1,2,:)),'o','color',[0.4940, 0.1840, 0.5560],'MarkerFaceColor',[0.4940, 0.1840, 0.5560],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(regr_slopes(2,2,:)),'o','color',[0.4660, 0.6740, 0.1880],'MarkerFaceColor',[0.4660, 0.6740, 0.1880],'markersize',3,'linestyle','-','linewidth',2)
plot(squeeze(regr_slopes(3,2,:)),'o','color',[0.3010, 0.7450, 0.9330],'MarkerFaceColor',[0.3010, 0.7450, 0.9330],'markersize',3,'linestyle','-','linewidth',2)
line([0 70],[0 0],'linestyle','--','color','k') % 0 line
line([0 70],[-0.02 -0.02],'color','r'); line([0 70],[0.027 0.027],'color','r') % region where corrs are not significant at p < 0.05
ylabel('regression slope');


print('nonstat_impact.pdf','-dpdf','-painters')%,'-bestfit')


%% Figure 8
% This figure is made in python using the script std_corr.py

%% Figures 9 and 10
% These figures are made in python, using the script ENSO_plotting.py

%% Figure 11

clear
load('DataFiles/ENSO_corrs.mat','SAT_corr_sig','pre_corr_sig')
load('DataFiles/ENSO_regression.mat','n34_SAT_regr','n34_pre_regr')
load('DataFiles/SAM_corrs.mat','model_pre_corr_land','model_SAT_corr_land')
land = ncread('DataFiles/sftlf_A1.static.nc','sftlf')';

figure(1)
% SAT
nonsig_count = zeros(3);
for w = [1 2 3]
    subplot(2,3,w)
    hold on
    xlabel('SAM-SAT r')
    ylabel('n34-SAM/SAT running corr regr. coeff')
    for i = 1:size(SAT_corr_sig,2)
        for j = 1:size(SAT_corr_sig,3)
            if isnan(SAT_corr_sig(w,i,j)) && ~isnan(land(i,j))  % Plot non-significant correlations first
                scatter(abs(model_SAT_corr_land(i,j)),abs(n34_SAT_regr(w,i,j)),'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
                nonsig_count(w,1) = nonsig_count(w,1) + 1;
            end           
        end
    end
end
% To make sure significant corrs are not covered, we need a separate
% loop...
sig_count = zeros(3);
for w = [1 2 3]
    subplot(2,3,w)
    for i = 1:size(SAT_corr_sig,2)
        for j = 1:size(SAT_corr_sig,3)
            if ~isnan(SAT_corr_sig(w,i,j)) && ~isnan(land(i,j))  % Plot non-significant correlations first
                scatter(abs(model_SAT_corr_land(i,j)),abs(n34_SAT_regr(w,i,j)),'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
                sig_count(w,1) = sig_count(w,1) + 1;
            end           
        end
    end
end
pcnt_sig_SAT = sig_count./(nonsig_count+sig_count);

% Precip
nonsig_count = zeros(3);
for w = [1 2 3]
    subplot(2,3,w+3)
    hold on
    xlabel('SAM-precip r')
    ylabel('n34-SAM/precip running corr regr. coeff')
    for i = 1:size(pre_corr_sig,2)
        for j = 1:size(pre_corr_sig,3)
            if isnan(pre_corr_sig(w,i,j)) && ~isnan(land(i,j))
                scatter(abs(model_pre_corr_land(i,j)),abs(n34_pre_regr(w,i,j)),'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
                nonsig_count(w,1) = nonsig_count(w,1) + 1;
            end           
        end
    end
end
sig_count = zeros(3);
for w = [1 2 3]
    subplot(2,3,w+3)
    hold on
    for i = 1:size(pre_corr_sig,2)
        for j = 1:size(pre_corr_sig,3)
            if ~isnan(pre_corr_sig(w,i,j)) && ~isnan(land(i,j))
                scatter(abs(model_pre_corr_land(i,j)),abs(n34_pre_regr(w,i,j)),'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
                sig_count(w,1) = sig_count(w,1) + 1;
            end           
        end
    end
end
pcnt_sig_precip = sig_count./(nonsig_count+sig_count);

h=gcf;
set(h,'PaperOrientation','landscape');
print(gcf,'figures/zscoreENSO_regr_SAM_r_scatter_2','-dpdf','-painters','-bestfit')

%% Figure A1 Probability of inclusion of non-stationary stations
% SAT
NUM_YRS = 500;
NUM_TRIALS = 1000;
for window = [31 61 91]
    if window == 31
        figure(1)
    elseif window == 61 
        figure(2)
    elseif window == 91
        figure(3)
    end
    load(['DataFiles/nonstat_map_',num2str(window),'yrwdw.mat'])
    NUM_GROUPS=6; i=1; NUM_CAL_WDW = 10; cmap = hsv(NUM_GROUPS); cmap(2,2)=0.8;
    clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
    for c=0:9
        CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window);
    end
    for group_size=[2,5,10,30,50,70]
        corrs = nan([1,NUM_TRIALS]);
        nstat_nstns = nan(NUM_CAL_WDW,NUM_TRIALS);
        for c=1:size(CAL_WDW,1)

            DIR_NAME = ['Proxies/',num2str(window),'yrWindow/'];
            load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(group_size),'stns_1000prox.mat'], ...
             'sat_lat','sat_lon');

            nstat_yrs = nan(size(sat_lat));
            for m=1:group_size
                for tr=1:NUM_TRIALS
                    nstat_yrs(tr,m) = nonstat_satmap(sat_lat(tr,m),sat_lon(tr,m));
                end
            end

            nstat_avyrs = mean(nstat_yrs,2);
            nstat_nstns(c,:) = sum(nstat_yrs > ceil(0.1*(NUM_YRS-window)),2)/group_size;
            
        end
        xbins = (0:1.0/group_size:1)-0.0000001; xbins(1)=0; xbins(end)=1;
        h=histc(nstat_nstns(:),xbins)/100;
        hold on; Hnd(i) = plot(xbins,h,'Color',cmap(i,:),'LineWidth',2);
        h(~h) = NaN;
        plot(xbins,h,'o','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:),'MarkerSize',7);
        i=i+1;
    end
        
    hold off;
    grid on
    xlim([0 1]); ylim([0 100]); xlabel(' Proportion of non-stationary proxies ','FontSize',16);
    legend([Hnd],'Network Size of 2','Network Size of 5','Network Size of 10','Network Size of 30','Network Size of 50','Network Size of 70');
    ylabel('Percentage of reconstructions','FontSize',16)


    set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'YLim'        , [0 100], ...
    'XLim'        , [0 1], ...
    'FontSize'    , 16 ...
    );

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 12]); %x_width=19cm y_width=28cm
set(gcf, 'PaperSize', [18 14]);
end

% Precip
for window = [31 61 91]
    if window == 31
        figure(4)
    elseif window == 61 
        figure(5)
    elseif window == 91
        figure(6)
    end
    load(['DataFiles/nonstat_map_',num2str(window),'yrwdw.mat'])
    NUM_GROUPS=6; i=1; NUM_CAL_WDW = 10; cmap = hsv(NUM_GROUPS); cmap(2,2)=0.8;
    clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
    for c=0:9
        CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window);
    end
    for group_size=[2,5,10,30,50,70]
        corrs = nan([1,NUM_TRIALS]);
        nstat_nstns = nan(NUM_CAL_WDW,NUM_TRIALS);
        for c=1:size(CAL_WDW,1)

            DIR_NAME = ['Proxies/',num2str(window),'yrWindow/'];
            load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/',num2str(group_size),'stns_1000prox.mat'], ...
             'stn_lat','stn_lon');

            nstat_yrs = nan(size(stn_lat));
            for m=1:group_size
                for tr=1:NUM_TRIALS
                    nstat_yrs(tr,m) = nonstat_precipmap(stn_lat(tr,m),stn_lon(tr,m));
                end
            end

            nstat_avyrs = mean(nstat_yrs,2);
            nstat_nstns(c,:) = sum(nstat_yrs > ceil(0.1*(NUM_YRS-window)),2)/group_size;
            
        end
        xbins = (0:1.0/group_size:1)-0.0000001; xbins(1)=0; xbins(end)=1;
        h=histc(nstat_nstns(:),xbins)/100;
        hold on; Hnd(i) = plot(xbins,h,'Color',cmap(i,:),'LineWidth',2);
        h(~h) = NaN;
        plot(xbins,h,'o','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:),'MarkerSize',7);
        i=i+1;
    end
        
    hold off;
    grid on
    xlim([0 1]); ylim([0 100]); xlabel(' Proportion of non-stationary proxies ','FontSize',16);
    legend([Hnd],'Network Size of 2','Network Size of 5','Network Size of 10','Network Size of 30','Network Size of 50','Network Size of 70');
    ylabel('Percentage of reconstructions','FontSize',16)

    set(gca, ...
    'TickDir','in', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 2,  ...    
    'YLim'        , [0 100], ...
    'XLim'        , [0 1], ...
    'FontSize'    , 16 ...
    );

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 12]); %x_width=19cm y_width=28cm
set(gcf, 'PaperSize', [18 14]);
end

%% Figure A2

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
