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
% SAT, Antarctica
figure(1)
dim = [0.13 0.82 0.1 0.1]; % Top left placement
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_AA.mat'],'sat_PC')
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
print('figures/ENSO_filtered_PCs_SAT_AA','-dpdf','-besfit')

% Precip, Antarctica
figure(2)
dim = [0.13 0.82 0.1 0.1]; % Top left placement
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_AA.mat'],'pre_PC')
    j = (floor(wdw/30)-1)*3;
    for i = 1:3
        subplot(3,3,j+i)
        plot(zscore(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))))
        hold on
        plot(zscore(pre_PC(2:end,i)))
        xlim([0 499-wdw]);
        legend({'Nino3.4','EOF'})
        b = annotation('textbox',dim,'String',['r = ',num2str(corrs_pre_n34_AA(floor(wdw/30),i))],'FitBoxToText','on');
        b.LineStyle = 'none';
    end
end

% Add overall title, add labels for year wdws and update legend for EOFs

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('figures/ENSO_filtered_PCs_PRE_AA','-dpdf','-besfit')

% SAT, SH
figure(3)
dim = [0.13 0.82 0.1 0.1]; % Top left placement
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_SH.mat'],'sat_PC')
    j = (floor(wdw/30)-1)*3;
    for i = 1:3
        subplot(3,3,j+i)
        plot(zscore(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))))
        hold on
        plot(zscore(sat_PC(2:end,i)))
        xlim([0 499-wdw]);
        legend({'Nino3.4','EOF'})
        b = annotation('textbox',dim,'String',['r = ',num2str(corrs_sat_n34_SH(floor(wdw/30),i))],'FitBoxToText','on');
        b.LineStyle = 'none';
    end
end

% Add overall title, add labels for year wdws and update legend for EOFs

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('figures/ENSO_filtered_PCs_SAT_SH','-dpdf','-besfit')

% Precip, SH
figure(4)
dim = [0.13 0.82 0.1 0.1]; % Top left placement
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_SH.mat'],'pre_PC')
    j = (floor(wdw/30)-1)*3;
    for i = 1:3
        subplot(3,3,j+i)
        plot(zscore(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))))
        hold on
        plot(zscore(pre_PC(2:end,i)))
        xlim([0 499-wdw]);
        legend({'Nino3.4','EOF'})
        b = annotation('textbox',dim,'String',['r = ',num2str(corrs_pre_n34_SH(floor(wdw/30),i))],'FitBoxToText','on');
        b.LineStyle = 'none';
    end
end

% Add overall title, add labels for year wdws and update legend for EOFs

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('figures/ENSO_filtered_PCs_PRE_SH','-dpdf','-besfit')

% SAT, Au
figure(5)
dim = [0.13 0.82 0.1 0.1]; % Top left placement
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_Au.mat'],'sat_PC')
    j = (floor(wdw/30)-1)*3;
    for i = 1:3
        subplot(3,3,j+i)
        plot(zscore(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))))
        hold on
        plot(zscore(sat_PC(2:end,i)))
        xlim([0 499-wdw]);
        legend({'Nino3.4','EOF'})
        b = annotation('textbox',dim,'String',['r = ',num2str(corrs_sat_n34_Au(floor(wdw/30),i))],'FitBoxToText','on');
        b.LineStyle = 'none';
    end
end

% Add overall title, add labels for year wdws and update legend for EOFs

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('figures/ENSO_filtered_PCs_SAT_Au','-dpdf','-besfit')

% Precip, Au
figure(6)
dim = [0.13 0.82 0.1 0.1]; % Top left placement
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_Au.mat'],'pre_PC')
    j = (floor(wdw/30)-1)*3;
    for i = 1:3
        subplot(3,3,j+i)
        plot(zscore(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))))
        hold on
        plot(zscore(pre_PC(2:end,i)))
        xlim([0 499-wdw]);
        legend({'Nino3.4','EOF'})
        b = annotation('textbox',dim,'String',['r = ',num2str(corrs_pre_n34_Au(floor(wdw/30),i))],'FitBoxToText','on');
        b.LineStyle = 'none';
    end
end

% Add overall title, add labels for year wdws and update legend for EOFs

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('figures/ENSO_filtered_PCs_PRE_Au','-dpdf','-besfit')

% SAT, SA
figure(7)
dim = [0.13 0.82 0.1 0.1]; % Top left placement
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_SA.mat'],'sat_PC')
    j = (floor(wdw/30)-1)*3;
    for i = 1:3
        subplot(3,3,j+i)
        plot(zscore(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))))
        hold on
        plot(zscore(sat_PC(2:end,i)))
        xlim([0 499-wdw]);
        legend({'Nino3.4','EOF'})
        b = annotation('textbox',dim,'String',['r = ',num2str(corrs_sat_n34_SA(floor(wdw/30),i))],'FitBoxToText','on');
        b.LineStyle = 'none';
    end
end

% Add overall title, add labels for year wdws and update legend for EOFs

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('figures/ENSO_filtered_PCs_SAT_SA','-dpdf','-besfit')

% Precip, Au
figure(8)
dim = [0.13 0.82 0.1 0.1]; % Top left placement
for wdw = [31 61 91]
    load(['DataFiles/reg_EOFs_',num2str(wdw),'yrwdw_SA.mat'],'pre_PC')
    j = (floor(wdw/30)-1)*3;
    for i = 1:3
        subplot(3,3,j+i)
        plot(zscore(nino_filt(floor(wdw/30),ceil(wdw/2):499-floor(wdw/2))))
        hold on
        plot(zscore(pre_PC(2:end,i)))
        xlim([0 499-wdw]);
        legend({'Nino3.4','EOF'})
        b = annotation('textbox',dim,'String',['r = ',num2str(corrs_pre_n34_SA(floor(wdw/30),i))],'FitBoxToText','on');
        b.LineStyle = 'none';
    end
end

% Add overall title, add labels for year wdws and update legend for EOFs

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 50 32],...
    'PaperSize',[50 32]);
print('figures/ENSO_filtered_PCs_PRE_SA','-dpdf','-besfit')





