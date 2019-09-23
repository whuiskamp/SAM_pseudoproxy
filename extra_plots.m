%% Setup

load DataFiles/model_output.mat

windowsize = 31; % The running window in years

corr_sat = zeros(45,size(sat_detr,3));
P_sat = zeros(45,size(sat_detr,3));
corr_precip = zeros(45,size(precip_detr,3));
P_precip = zeros(45,size(precip_detr,3));

for i=1:45
    for j=1:size(sat_detr,3)
        [corr_sat(i,j) P_sat(i,j)]= corr(SAM',sat_detr(:,i,j));
    end
end

corr_precip = zeros(size(45),size(precip_detr,3));
for i=1:45
    for j=1:size(precip_detr,3)
        corr_precip(i,j) = corr(SAM',precip_detr(:,i,j));
    end
end

NUM_YRS=500; NUM_TRIALS=1000; numstnstocompare=2:70;


%% Figure 1-1
load(['DataFiles/runcorr31yrwdw.mat'])
load(['DataFiles/nonstat_map_31yrwdw.mat'])
precip_runcorr_91 = load('DataFiles/runcorr91yrwdw.mat','precip_runcorr');
precip_runcorr_91 = precip_runcorr_91.('precip_runcorr');
sat_runcorr_91 = load('DataFiles/runcorr91yrwdw.mat','sat_runcorr');
sat_runcorr_91 = sat_runcorr_91.('sat_runcorr');
load(['DataFiles/nonstat_map_91yrwdw.mat'])

% Temp
bin = -1.0:0.03:1.0;
[sorted_corr sorted_corr_ind] = sort(corr_sat(:));
sat_runcorr = sat_runcorr(:,1:45,:);
sat_runcorr_91 = sat_runcorr_91(:,1:45,:);
sorted_sat_runcorr = sat_runcorr(:,sorted_corr_ind); % This will need to be checked
sorted_sat_runcorr_91 = sat_runcorr_91(:,sorted_corr_ind);
bin_sizes = histc(squeeze(sorted_corr),bin);
current_index = 1; sat_runcorr_quan = nan(7,length(bin_sizes)); 
for m=1:length(bin_sizes)
    sat_runcorr_quan(:,m) = quantile(reshape(sorted_sat_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end
current_index = 1; sat_runcorr_91_quan = nan(7,length(bin_sizes));
for m=1:length(bin_sizes)
    sat_runcorr_91_quan(:,m) = quantile(reshape(sorted_sat_runcorr_91(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes(m);
end


% Precip

[sorted_corr_p sorted_corr_ind_p] = sort(corr_precip(:));
precip_runcorr = precip_runcorr(:,1:45,:);
precip_runcorr_91 = precip_runcorr_91(:,1:45,:);
sorted_precip_runcorr = precip_runcorr(:,sorted_corr_ind_p); % This will need to be checked
sorted_precip_runcorr_91 = precip_runcorr_91(:,sorted_corr_ind_p);
bin_sizes_p = histc(squeeze(sorted_corr_p),bin);
current_index = 1; precip_runcorr_quan = nan(7,length(bin_sizes_p));
for m=1:length(bin_sizes_p)
    precip_runcorr_quan(:,m) = quantile(reshape(sorted_precip_runcorr(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes_p(m);
end
current_index = 1; precip_runcorr_quan_91 = nan(7,length(bin_sizes_p));
for m=1:length(bin_sizes_p)
    precip_runcorr_quan_91(:,m) = quantile(reshape(sorted_precip_runcorr_91(:,current_index:(current_index-1+bin_sizes(m))),[],1), ...
                 [0.01,0.05,0.25,0.5,0.75,0.95,0.99]);
    current_index = current_index + bin_sizes_p(m);
end

%% The plotting part

figure(1)
hold on; Hnd = nan(1,7);
for n=1:7
    HA(n) = plot(bin,sat_runcorr_quan(n,:));
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
text(-1+0.08,1-0.18,'a)','FontSize',22,'FontWeight','bold');

figure(2)

hold on; Hnd = nan(1,7);
for n=1:7
    HA(n) = plot(bin,precip_runcorr_quan(n,:));
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
text(-1+0.08,1-0.18,'b)','FontSize',22,'FontWeight','bold');

figure(3)

hold on; Hnd = nan(1,7);
for n=1:7
    HA(n) = plot(bin,sat_runcorr_91_quan(n,:));
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
ylabel(['r(91yr)'],'FontSize',16);
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
text(-1+0.08,1-0.18,'c)','FontSize',22,'FontWeight','bold');

figure(4)

hold on; Hnd = nan(1,7);
for n=1:7
    HA(n) = plot(bin,precip_runcorr_quan_91(n,:));
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
ylabel(['r(91yr)'],'FontSize',16);
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
text(-1+0.08,1-0.18,'d)','FontSize',22,'FontWeight','bold');


% Figure 1-2
subplot(2,2,2);

contourf(lon,lat,corr_sat); axis equal; plotworld

% Observations overlay
load DataFiles/GISTEMP_ersst_n34Corr.mat
hold on; [c_m, c_h]=contour(lon1,lat1,corr_ssta_an,[0, 0.2, 0.4, 0.6, 0.8 1]); axis equal; hold off;
hold on; [c_m, c_h]=contour(lon1,lat1,corr_ssta_an,[-0.2, -0.4, -0.6, -0.8 -1]); axis equal; hold off;

shading flat
title('r(ts,Nino3.4)','FontSize',14,'FontWeight','bold');
plotworld;
colormap(b2r(-1,1));
set(gca, ...
    'TickDir','out', ...
    'Box', 'on',    ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , -90:30:90, ...
    'YTickLabel'  , {'90S','60S','30S','EQ','30N','60N','90N'} ,...
    'LineWidth'   , 2,  ...    
    'YLim'        , [-90 90], ...
    'XLim'        , [0 360],  ...
    'FontSize'    , 16   ...
    );
set(gcf,'renderer','zbuffer');
colorbar;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 19 19]);
set(gcf, 'PaperSize', [35 13]);
text(0+1,90-1,'b)','FontSize',22,'FontWeight','bold');

set(gcf, 'PaperSize', [35 13]);

%% Probability of inclusion of non-stationary stations
% SAT
NUM_YRS = 500;
NUM_TRIALS = 1000;
for window = [31 91]
    if window == 31
        figure(1)
    elseif window == 91 figure (2)
    end
    load(['DataFiles/nonstat_map_',num2str(window),'yrwdw.mat'])
    NUM_GROUPS=6; i=1; NUM_CAL_WDW = 10; cmap = hsv(NUM_GROUPS); cmap(2,2)=0.8;
    clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
    for c=0:9
        CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
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
    % h = findobj(gca,'Type','Patch');
    % set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
    ylabel('Percentage of reconstructions','FontSize',16)
%     title('Proportion of nstat stns per reconstruction','FontSize',14);

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
for window = [31 91]
    if window == 31
        figure(3)
    elseif window == 91 
        figure(4)
    end
    load(['DataFiles/nonstat_map_',num2str(window),'yrwdw.mat'])
    NUM_GROUPS=6; i=1; NUM_CAL_WDW = 10; cmap = hsv(NUM_GROUPS); cmap(2,2)=0.8;
    clear CAL_WDW; overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*window)/9.0);
    for c=0:9
        CAL_WDW(c+1,:) = (1+c*(window-overlap)):((c*(window-overlap))+window); %#ok<SAGROW>
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
    % h = findobj(gca,'Type','Patch');
    % set(h,'FaceColor',[1 1 1], 'EdgeColor','black');
    ylabel('Percentage of reconstructions','FontSize',16)
%     title('Proportion of nstat stns per reconstruction','FontSize',14);

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

%% Save stuff to netcdf

nccreate('corr.nc','corr_precip','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat',corr_sat(1:45,1:144))

load('Proxies/31yrWindow/CalWdw:1-31/1stns_1000prox')

nccreate('corr.nc','corr_precip_wdw1','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip_wdw1',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat_wdw1','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat_wdw1',corr_sat(1:45,1:144))

load('Proxies/31yrWindow/CalWdw:53-83/1stns_1000prox')

nccreate('corr.nc','corr_precip_wdw2','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip_wdw2',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat_wdw2','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat_wdw2',corr_sat(1:45,1:144))

load('Proxies/31yrWindow/CalWdw:105-135/1stns_1000prox')

nccreate('corr.nc','corr_precip_wdw3','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip_wdw3',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat_wdw3','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat_wdw3',corr_sat(1:45,1:144))

load('Proxies/31yrWindow/CalWdw:157-187/1stns_1000prox')

nccreate('corr.nc','corr_precip_wdw4','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip_wdw4',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat_wdw4','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat_wdw4',corr_sat(1:45,1:144))

load('Proxies/31yrWindow/CalWdw:209-239/1stns_1000prox')

nccreate('corr.nc','corr_precip_wdw5','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip_wdw5',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat_wdw5','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat_wdw5',corr_sat(1:45,1:144))

load('Proxies/31yrWindow/CalWdw:261-291/1stns_1000prox')

nccreate('corr.nc','corr_precip_wdw6','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip_wdw6',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat_wdw6','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat_wdw6',corr_sat(1:45,1:144))

load('Proxies/31yrWindow/CalWdw:313-343/1stns_1000prox')

nccreate('corr.nc','corr_precip_wdw7','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip_wdw7',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat_wdw7','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat_wdw7',corr_sat(1:45,1:144))

load('Proxies/31yrWindow/CalWdw:365-395/1stns_1000prox')

nccreate('corr.nc','corr_precip_wdw8','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip_wdw8',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat_wdw8','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat_wdw8',corr_sat(1:45,1:144))

load('Proxies/31yrWindow/CalWdw:417-447/1stns_1000prox')

nccreate('corr.nc','corr_precip_wdw9','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip_wdw9',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat_wdw9','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat_wdw9',corr_sat(1:45,1:144))

load('Proxies/31yrWindow/CalWdw:469-499/1stns_1000prox')

nccreate('corr.nc','corr_precip_wdw10','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_precip_wdw10',corr_precip(1:45,1:144))
nccreate('corr.nc','corr_sat_wdw10','Datatype','single','Dimensions',...
{'latitude',45,'longitude',144},'Format','classic')
ncwrite('corr.nc','corr_sat_wdw10',corr_sat(1:45,1:144))

nccreate('corr.nc','latitude','Dimensions',...
    {'latitude',45},'Format','classic')
nccreate('corr.nc','longitude','Dimensions',...
{'longitude',144},'Format','classic')

ncwrite('corr.nc','latitude',lat(1:45,1))
ncwrite('corr.nc','longitude',lon(1:144,1))
ncwriteatt('corr.nc','latitude','axis','Y')
ncwriteatt('corr.nc','longitude','axis','X')

ncwriteatt('corr.nc','latitude','axis','Y')
ncwriteatt('corr.nc','latitude','standard_name','latitude')
ncwriteatt('corr.nc','latitude','units','degrees_north')

ncwriteatt('corr.nc','longitude','axis','Y')
ncwriteatt('corr.nc','longitude','standard_name','longitude')
ncwriteatt('corr.nc','longitude','units','degrees_east')

%% Skilfullness of different reconstructions
% First for precip
% Do this for CPS_AAo; AuNZ; SA; SoA
clf
clear

s_Hnd = tight_subplot(3,4,[0.05 0.01],[0.10 0.01],[0.1 0.01]);
numstnstocompare = 2:70; skilful_threshold = 0.6;
windowsize = [31 61 91];
pcnt31_skilful_CPS = nan(max(numstnstocompare),length([31 61 91])); % Percentage of reconstructions that get the same as the 500yr recons
pcnt61_skilful_CPS = nan(max(numstnstocompare),length([31 61 91]));
pcnt91_skilful_CPS = nan(max(numstnstocompare),length([31 61 91]));

GROUP_NAME = 'pcnt_ts'; % Change group name to get other figs
NUM_CAL_WDW = 10; clear CAL_WDW; NUM_YRS = 500; NUM_TRIALS = 1000; 

for windowsize = [31 61 91]
    for region = 1:5
        DIR_NAME = ['/srv/ccrc/data37/z3215716/SAM/Proxies/',num2str(windowsize),'yrWindow/'];
        overlap = ceil(-(NUM_YRS-NUM_CAL_WDW*windowsize)/9.0);
        for c=0:9
            CAL_WDW(c+1,:) = (1+c*(windowsize-overlap)):((c*(windowsize-overlap))+windowsize); %#ok<SAGROW>
        end

        precip_corr_CPS = nan(3,max(numstnstocompare),size(CAL_WDW,1),NUM_TRIALS,'single');

        for c=1:size(CAL_WDW,1)
            load([DIR_NAME,'/CalWdw:',num2str(CAL_WDW(c,1)),'-',num2str(CAL_WDW(c,end)),'/tonsofstats.mat'], ...
            'all_stn_corr*')
            if region == 1
				precip_corr_CPS(floor(windowsize/30),:,c,:) = all_stn_corr_CPS_SoA;
			elseif region == 2
				precip_corr_CPS(floor(windowsize/30),:,c,:) = all_stn_corr_CPS_SA;
			elseif region == 3
				precip_corr_CPS(floor(windowsize/30),:,c,:) = all_stn_corr_CPS_AuNz;
			elseif region == 4
				precip_corr_CPS(floor(windowsize/30),:,c,:) = all_stn_corr_CPS_AAo;
			elseif region == 5
				precip_corr_CPS(floor(windowsize/30),:,c,:) = all_stn_corr_CPS;
            end
        end


        for n=numstnstocompare
            pcnt31_skilful_CPS(n,floor(windowsize/30)) = sum(sum(sum(...
            squeeze(precip_corr_CPS(floor(windowsize/30),n,:,:))>skilful_threshold )))/(NUM_CAL_WDW*NUM_TRIALS);
        end
        
        corr_CPS_qn = quantile(precip_corr_CPS,[.05 .5 .95], 4);  % This creates a (3,70,10,3) file ie: 3 windows, 70 sites, 10 cal windows and 3 percentiles
        corr_CPS_qn_rng = nan(3,size(corr_CPS_qn,2),2,size(corr_CPS_qn,4)); % creates the range for the 3 windows, 70 sites, 2 bounds of the range, 3 percentiles
        corr_CPS_qn_rng(:,:,1,:) = min(corr_CPS_qn,[],3);
        corr_CPS_qn_rng(:,:,2,:) = max(corr_CPS_qn,[],3);
                              
        % Plotting
        
        if windowsize == 31
            skip = 0;
        elseif windowsize == 61
            skip = 4;
        elseif windowsize == 91
            skip = 8;
        end
        
        if region == 5
            for i = 1:12
                axes(s_Hnd(i)); hold on; xlim([0,70]); ylim([0,1]); grid on
                plot(pcnt31_skilful_CPS(:,floor(windowsize/30)),'Color',[0 0 0.8],'LineWidth',3); 
                jbfill([2:70],squeeze(corr_CPS_qn_rng(2:70,2,3))',squeeze(corr_CPS_qn_rng(2:70,1,1))','b','k',[],0.5);
            end
        elseif any(region == 1:4) 
               for i = region+skip
                    axes(s_Hnd(i));
                    %subplot(3,4,i);
                    plot(pcnt31_skilful_CPS(:,floor(windowsize/30)),'Color','red','LineWidth',3);
                    jbfill([2:70],squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,2,3))',squeeze(corr_CPS_qn_rng(floor(windowsize/30),2:70,1,1))','y','k',[],0.5);
                end
         end
    end
end

% Formatting
axes(s_Hnd(1)); title(['EPC\_RV'],'FontSize',16);
axes(s_Hnd(2)); title(['CPS\_RV'],'FontSize',16);
axes(s_Hnd(3)); title(['MRV'],'FontSize',16);
axes(s_Hnd(4)); title(['RVM'],'FontSize',16);
leg_h = legend(s_Hnd(3),[prop_31_Hnd(1,1) prop_61_Hnd(1,1) prop_91_Hnd(1,1)], ...
               'PNEOF 31yr','PNEOF 61yr','PNEOF 91yr','orientation','horizontal');
set(leg_h,'FontSize',14,'orientation','horizontal');
leg_h = legend(s_Hnd(4),[prop_31_Hnd(1,2) prop_61_Hnd(1,2) prop_91_Hnd(1,2)], ...
               'STAT 31yr','STAT 61yr','STAT 91yr','orientation','horizontal');
set(leg_h,'FontSize',14,'orientation','horizontal');


for i=1:3
    axes(s_Hnd(i));
    set(gca,'YTickLabel',[],'XTickLabel',[])
    set(gca, 'FontSize',16, 'LineWidth', 2.0, 'Box', 'on', 'YTick', [0:0.2:1],'XTick', [0:20:70]); 
    xlim([0 70]); ylim([0 1]); grid on
    text(0+1,1-0.06,[letters(i),')'],'FontSize',22,'FontWeight','bold');
end

axes(s_Hnd(1)); xlabel('Network Size');
for i=1:1
    axes(s_Hnd(i*2-1));
    set(gca,'YTickLabel',[0:0.2:1]);
    ylabel(['Proportion'])
end

for i=1:4
    axes(s_Hnd(i));
    set(gca,'XTickLabel',[0:20:70]);
end

suptitle('PNEOF1');
