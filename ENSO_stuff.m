% Calculating Nino3.4 and comparing with SAM/SAT running corr. PCs
% Data has been downloaded from:
% ftp://nomads.gfdl.noaa.gov/gfdl_cm2_1/CM2.1U_Control-1860_D4/pp/ocean_tripolar/ts/monthly/
% It was combined with ncarcat then DJF was extraced with cdo select,season=DJF temp_tot.nc temp_DJF.nc
% lat and lon come from these files (geolat, geolon).

% Import file and remove first Jan/Feb and last Dec.
% Import lat and lon grids. We only need values at the equator, so extract
% a vector for these.

%data = ncread('DataFiles/monthly_data/temp_DJF.nc','thetao',[1,1,1,1],[360,200,1,1500]);
%sst = squeeze(data); sst(:,:,[1 2 1500]) = [];   
load('DataFiles/model_output.mat','sst_DJF');
load('DataFiles/model_output.mat','lat_oce','lon_oce'); % For reference

% Create DJF means
sst_DJF = nan(499,360,200);
for i = 1:3:1497
    sst_DJF(ceil(i/3),:,:) = mean(sst(:,:,i:i+2),3);
end

%save('DataFiles/model_output.mat','sst_DJF','-append');

% Nino3.4 Boundaries from Ryan's script.
EN_S = 84; EN_N = 109; EN_W = 191; EN_E = 241;

n34_ind = mean(mean(sst_DJF(:,EN_W:EN_E,EN_S:EN_N),3),2);

% Detrend and standardise

n34 = detrend(n34_ind,'linear');
n34 = zscore(n34);

save('DataFiles/model_output.mat','n34','-append');

%% 1: Low pass filtering of ENSO index (method from Shayne's email)
% This is done on the cluster - data too big to run locally
load('DataFiles/monthly_data/n34_monthly.mat','n34_m')
for wdw = [30 60 90]
    if wdw == 30
        % 30 -year filter would be:
        f1 = 1.0717e-9; freq = 3.8580e-7; Rp = 1; Rs = 5; falloff = .2e-9;
    elseif wdw == 60
        % 60-year filter would be:
        f1 = 5.3584e-10; freq = 3.8580e-7; Rp = 1; Rs = 5 ; falloff = .9e-10;
    elseif wdw == 90
        % 90 year filter would be:
        f1 = 3.5722e-10; freq = 3.8580e-7; Rp = 1; Rs = 5 ; falloff = 0.6e-10;
    end
    nyq = freq/2.0;
    l_Wp = f1 / nyq;
    l_Ws = (f1-falloff)/nyq;  
    [low_n,low_Wn] = buttord(l_Wp,l_Ws,Rp,Rs);
    [low_b,low_a] = butter(low_n,low_Wn);
    if wdw == 30
        filt_n34_30 = filtfilt(low_b,low_a,double(n34_m));
    elseif wdw == 60
        filt_n34_60 = filtfilt(low_b,low_a,double(n34_m));
    elseif wdw == 90
        filt_n34_90 = filtfilt(low_b,low_a,double(n34_m)); 
    end        
end

% Save and do correlations with this
for i = 1:499
    Start = i*12;
    End = Start+2;
    filt_ann_30(i) = mean(filt_n34_30(Start:End));
    filt_ann_60(i) = mean(filt_n34_60(Start:End));
    filt_ann_90(i) = mean(filt_n34_90(Start:End));      
end

save('DataFiles/monthly_data/filtered_n34.mat','filt_ann_30','filt_ann_60','filt_ann_90')
%% 2: Compare Nino3.4 index to Principal components derived from the running correlation
% of SAM and SAT/PRECIP over geographical regions
load('DataFiles/model_output.mat')
load('DataFiles/filtered_n34.mat')
% First, create running mean of Nino3.4
wdw = [31 61 91];
for i = wdw
   nino_smth(floor(i/30),:) = movmean(n34,i,'Endpoints','fill');
end

nino_filt(1,:) = filt_ann_30'; nino_filt(2,:) = filt_ann_60'; nino_filt(3,:) = filt_ann_90';
% Compare to PCs from AA
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_AA.mat'],'pre_PC','sat_PC')
    for i = 1:3
        corrs_pre_n34_AA(floor(wdw/30),i) = corr2(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))',pre_PC(2:end,i)); % starts at 2 because ENSO record is one year shorter than our running corr record
        corrs_sat_n34_AA(floor(wdw/30),i) = corr2(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))',sat_PC(2:end,i));
    end
end

% Compare to PCs from SH
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_SH.mat'],'pre_PC','sat_PC')
    for i = 1:3
        corrs_pre_n34_SH(floor(wdw/30),i) = corr2(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))',pre_PC(2:end,i));
        corrs_sat_n34_SH(floor(wdw/30),i) = corr2(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))',sat_PC(2:end,i));
    end
end

% Compare to PCs from Au/NZ 
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_Au.mat'],'pre_PC','sat_PC')
    for i = 1:3
        corrs_pre_n34_Au(floor(wdw/30),i) = corr2(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))',pre_PC(2:end,i));
        corrs_sat_n34_Au(floor(wdw/30),i) = corr2(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))',sat_PC(2:end,i));
    end
end

% Compare to PCs from South America
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_SA.mat'],'pre_PC','sat_PC')
    for i = 1:3
        corrs_pre_n34_SA(floor(wdw/30),i) = corr2(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))',pre_PC(2:end,i));
        corrs_sat_n34_SA(floor(wdw/30),i) = corr2(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))',sat_PC(2:end,i));
    end
end


% Plotting 
figure(2)
dim = [0.13 0.82 0.1 0.1]; % Top left placement
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_AA.mat'],'pre_PC','sat_PC')
    j = (floor(wdw/30)-1)*3;
    for i = 1:3
        subplot(3,3,j+i)
        plot(zscore(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))))
        hold on
        plot(zscore(sat_PC(2:end,i)))
        xlim([0 499-wdw]);
        legend({'Nino3.4','EOF'})
        b = annotation('textbox',dim,'String',['r = ',num2str(corrs_sat_n34_AA(floor(wdw/30),i))],'FitBoxToText','on');
        b.LineStyle = 'none';
    end
end

% Add overall title, add labels for year wdws and update legend for EOFs

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('ENSO_filtered_PCs','-dpdf','-besfit')





% Test for 90 year band

%f1 = 3.5722e-10; freq = 3.215e-8; Rp = 1; Rs = 5 ; falloff = .1e-9; 



% Test figure comparing methods
figure(2)
plot(zscore(sat_PC(:,2)),'linewidth',2,'color','red')
hold
plot(zscore(nino_smth(3,46:454)))
plot(zscore(filt_n34_90(46:454)))
plot(zscore(n34_90_loess(46:454)))


figure(1)
subplot(3,1,1)
plot(nino_smth(1,:))
hold
plot(filt_n34_30(1:484,1))
legend({'31yr runn. mean.','30yr filter'})
subplot(3,1,2)
plot(nino_smth(2,:))
hold
plot(filt_n34_60(1:484,1))
legend({'61yr runn. mean.','60yr filter'})
subplot(3,1,3)
plot(nino_smth(3,:))
hold
plot(filt_n34_90(1:484,1))
legend({'91yr runn. mean.','90yr filter'})

set(gcf,...
    'PaperPosition',[0 0 15 32],...
    'PaperSize',[15 32]);
print('n34_filtered_method_2','-dpdf','-besfit')

set(gcf,...
    'PaperPosition',[0 0 32 15],...
    'PaperSize',[32 15]);
print('method_compar','-dpdf','-besfit')

%% Lets try the filtering with monthly data...
clear
load('DataFiles/sst_monthly_DJF.mat')
EN_S = 84; EN_N = 109; EN_W = 191; EN_E = 241;

n34_ind = squeeze(mean(mean(sst(EN_W:EN_E,EN_S:EN_N,:),2),1));
n34_m = detrend(n34_ind,'linear');

for wdw = 90
    if wdw == 10
        % Decadal filter
        f1 = 3.215e-9; freq = 3.8580e-7; Rp = 1 ; Rs = 5 ; falloff = .5e-9; % This sampling freq. should actually be multiplied by 3
    elseif wdw == 30
        % 30 -year filter would be:
        f1 = 1.0717e-9; freq = 3.8580e-7; Rp = 1; Rs = 5; falloff = .7e-10;
    elseif wdw == 60
        % 60-year filter would be:
        f1 = 5.3584e-10; freq = 3.8580e-7; Rp = 1; Rs = 5 ; falloff = .5e-10;
    elseif wdw == 90
        % 90 year filter would be:
        f1 = 1.07166e-9; freq = 11574e-6; Rp = 1; Rs = 5 ; falloff = .5e-10;
    end
end

nyq = freq/2.0;
l_Wp = f1 / nyq;
l_Ws = (f1-falloff)/nyq;  
[low_n,low_Wn] = buttord(l_Wp,l_Ws,Rp,Rs);
[low_b,low_a] = butter(low_n,low_Wn);
filt_n34_90 = filtfilt(low_b,low_a,double(n34_m));
%[sos, g] = zp2sos(z,low_b,low_a);
%freqz(sos)
%zplane(z,low_b)
%fvtool(low_b,low_a,'FrequencyScale','log')

% My additions
% figure(wdw)
% [h,w]=freqz(low_b,low_a);
% p=abs(h);
% f=w/(pi*2);
% plot(f,p)


% Create DJF means
n34_DJF_m = nan(499,1);
for i = 1:3:1497
    n34_DJF_m(ceil(i/3),:) = mean(filt_n34_90(i:i+2,1));
end

figure(1)
subplot(3,1,1)
plot(nino_smth(1,:))
hold
plot(filt_n34_30(1:484,1))
legend({'31yr runn. mean.','30yr filter'})
subplot(3,1,2)
plot(nino_smth(2,:))
hold
plot(filt_n34_60(1:484,1))
legend({'61yr runn. mean.','60yr filter'})
subplot(3,1,3)
plot(nino_smth(3,:))
hold
plot(filt_n34_90(1:484,1))
legend({'91yr runn. mean.','90yr filter'})

set(gcf,...
    'PaperPosition',[0 0 15 32],...
    'PaperSize',[15 32]);
print('proxies_Marshall_regressions','-dpdf','-besfit')
print('n34_filtered_compare','-dpdf','-besfit')

% Figures comparing methods:

figure(1)
subplot(3,1,1)
plot(n34)
hold on
plot(n34_31_loess,'linewidth',2,'color','r')
plot(filt_n34_30,'linewidth',2,'color','k')
plot(nino_smth(1,:),'linewidth',2,'color','magenta')
legend({'n34','30yr rloess filt','30yr bandpass filt','31yr run mean'})

subplot(3,1,2)
plot(n34)
hold on
plot(n34_61_loess,'linewidth',2,'color','r')
plot(filt_n34_60,'linewidth',2,'color','k')
plot(nino_smth(2,:),'linewidth',2,'color','magenta')
legend({'n34','60yr rloess filt','60yr bandpass filt','61yr run mean'})

subplot(3,1,3)
plot(n34)
hold on
plot(n34_91_loess,'linewidth',2,'color','r')
plot(filt_n34_90,'linewidth',2,'color','k')
plot(nino_smth(3,:),'linewidth',2,'color','magenta')
legend({'n34','90yr rloess filt','90yr bandpass filt','91yr run mean'})

set(gcf,...
    'PaperPosition',[0 0 15 32],...
    'PaperSize',[15 32]);
print('methods_compare_2','-dpdf','-besfit')



% Annual data methods compare
plot(n34_91_rloess,'linewidth',2,'color',[0.85 0.33 0.1])
hold
plot(filt_n34_90,'linewidth',2,'color','cyan')
plot(n34_91_loess,'linewidth',2,'color',[0 0.45 0.74])
plot(nino_smth(3,:),'color','magenta')
line([45 45],[-0.3 0.4],'linestyle','--','color','k')
line([455 455],[-0.3 0.4],'linestyle','--','color','k')
legend({'rloess','90yr filter','loess','runn. mean.'})

set(gcf,...
    'PaperPosition',[0 0 32 15],...
    'PaperSize',[32 15]);
print('methods_compare_ann','-dpdf','-besfit')










