% Calculating Nino3.4 and comparing with SAM/SAT running corr. PCs
% Data has been downloaded from:
% ftp://nomads.gfdl.noaa.gov/gfdl_cm2_1/CM2.1U_Control-1860_D4/pp/ocean_tripolar/ts/monthly/
% It was combined with ncarcat then DJF was extraced with cdo select,season=DJF temp_tot.nc temp_DJF.nc

% Import files and extract DJF mean. Files are decadal blocks

data = ncread('DataFiles/monthly_data/temp_DJF.nc','thetao')';
sst_DJF = data;        
        
% Nino3.4 Boundaries from Ryan's script.
[~,nN] = min(abs(lat-5));
[~,nS] = min(abs(lat+5));
[~,nE] = min(abs(lon-240));
[~,nW] = min(abs(lon-190));

SOMETHING IN BETWEEN

n34_ind = mean(mean(ats(:,nS:nN,nW:nE),3),2);