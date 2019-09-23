%% Plotting
% Set up map
for i=1:3
    figure(i)
    axesm('MapProjection','stereo','origin',[-90,0],'MapLatLimit',[-90 -20])
    load coast
    framem
    plotm(lat,long,'k','linewidth',3)
    gridm %If you want gridlines
end
clear lat long

%% Spatial correlations
clear
ncid = netcdf.open('pr_A0.0001-0500.nc','NC_NOWRITE');
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i=0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
assignin('base', char(varname), netcdf.getVar(ncid,i));
end

ncid = netcdf.open('ua_A0.0001-0500.nc','NC_NOWRITE');
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i=0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
assignin('base', char(varname), netcdf.getVar(ncid,i));
end

ncid = netcdf.open('uas_A0.0001-0500.nc','NC_NOWRITE');
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
for i=0:numvars-1
[varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
assignin('base', char(varname), netcdf.getVar(ncid,i));
end

lon(1) = 0;
lon(144) = 360;

wind_s = uas;
wind = squeeze(ua(:,:,3,:)); %for 850hPa
precip = pr;

% Correlation with 850hPa winds

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(wind(j,i,1:500),squeeze(precip(j,i,1:500)));
     wind_corr(i,j)=corr(1,2); 
     wind_p(i,j)=pval(1,2);
  end
end

[row,col] = find(wind_p>0.05);
wind_corr2 = wind_corr;

for i = 1:length(col)
    wind_corr2(row(i),col(i)) = NaN;
end

% Correlation with surface wind speed

for i=1:length(lat),
  for j=1:length(lon),
     [corr,pval]=corrcoef(wind_s(j,i,1:500),squeeze(precip(j,i,1:500)));
     wind_s_corr(i,j)=corr(1,2); 
     wind_s_p(i,j)=pval(1,2);
  end
end

[row,col] = find(wind_s_p>0.05);
wind_s_corr2 = wind_s_corr;

for i = 1:length(col)
    wind_s_corr2(row(i),col(i)) = NaN;
end

figure(1)
wind_corr3 = double(wind_corr2);
contourfm(lat,lon,wind_corr3,30,'levellist',[-0.1 -0.2 -0.3 -0.4 -0.5 -0.7 0.1 0.2 0.3 0.4 0.5 0.7])
hcb=colorbar;
set(hcb,'YTick',-0.7:0.1:0.7)

figure(2)
wind_s_corr3 = double(wind_s_corr2);
contourfm(lat,lon,wind_s_corr3,30,'levellist',[-0.1 -0.2 -0.3 -0.4 -0.5 -0.7 0.1 0.2 0.3 0.4 0.5 0.7])
hcb=colorbar;
set(hcb,'YTick',-0.7:0.1:0.7)

% Running correlations

window = [11 21 51 71];
% output is of dimensions (144,45,500,4)
wipe = zeros(144,45);

for i = 1:144
    i
    for j = 1:45 %SH only
       % j
        for k = 1:4
            indexColumn = 1;
            windowSize = window(k); 
            dataMatrix(1:500,1) = squeeze(wind(i,j,:));
            dataMatrix(1:500,2) = squeeze(precip(i,j,:));
            [correlationTS correlationTS2 PvalTS] = movingCorrelation(dataMatrix, windowSize, indexColumn);
            correlation(i,j,1:500,k) = correlationTS;
            pval(i,j,1:500,k) = PvalTS;
            correlation2(i,j,1:500,k) = correlationTS2;
        end
    end
end

figure(1)
% 21 year  window
data = correlation2(:,:,:,2);
for i=11:489
    i
    contourfm(lat(1:45),lon,squeeze(correlation2(:,:,i,1))',30,'levellist',[-0.1 -0.2 -0.3 -0.4 -0.5 -0.7 0.1 0.2 0.3 0.4 0.5 0.7]);
    F(i) = getframe(gcf);
    h = get(gca,'Children');
    delete(h(1,1));
end

movie2avi(F,'correlation21.avi','fps',10,'quality',30);  %if it doesn't work, you've fucked up the plotting (blank plots).

h = get(gca,'Children');
delete(h(1,1));
contourfm(lat(1:45),lon,squeeze(correlation(:,:,21,2)'),30,'levellist',[-0.1 -0.2 -0.3 -0.4 -0.5 -0.7 0.1 0.2 0.3 0.4 0.5 0.7]);
