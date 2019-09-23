% This matlab script will produce synthetic running correlations for
% precipitation and temperature by using the same method as specified in
% {Gallant et al, 2013, 'Nonstationary Australasian Teleconnections and
%   Implications for Paleoclimate Reconstructions', Journal of Climate,
%   vol. 26 pp. 8827-8849} as equation (1):

% nu(t) = a0 + a1*c(t) + sigma_nu*sqrt(1-r^2)*[eta_nu(t) + B*eta_nu(t-1)]

% where: nu(t) is the synthetic precipitation/temperature series
%        a0 is the first regression coefficient (also the mean of the data)
%        a1 is the second regression coefficient
%        c(t) is the Nino3.4 Index
%        sigma_nu is the standard deviation of the actual precip/temp series
%        r is the correlation between Nino3.4 index and precip/temp
%        eta_nu is random Gaussian noise
%        B is the autocorrelation of the climate variable at lag 1
%        [eta_nu(t) + B*eta_nu(t-1)] is red noise

% Note: This script requires the mexcdf package to be installed, and the
% following files: plotworld.m, coast_v2.mat, b2r.m to be in the current
% directory and in DataFiles

%% Setup

load DataFiles/model_output.mat

reg_wind_detr = nan*zeros(size(wind_detr,2),size(wind_detr,3)); % Regress the SAM index onto the wind field
for i=1:size(wind_detr,2)
    for j=1:size(wind_detr,3)
        if ~isnan((wind_detr(:,i,j)))
        reg_wind_detr(i,j) = regress(wind_detr(:,i,j),SAM);
        end
    end
end

reg_precip_detr = nan*zeros(size(precip_detr,2),size(precip_detr,3)); % Regress the SAM index onto the precip field
for i=1:size(precip_detr,2)
    for j=1:size(precip_detr,3)
        if ~isnan((precip_detr(:,i,j)))
        reg_precip_detr(i,j) = regress(precip_detr(:,i,j),SAM);
        end
    end
end

corr_precip = nan*zeros(size(precip_detr,2),size(precip_detr,3)); % Correlate SAM and precip
for i=1:size(precip_detr,2)
    for j=1:size(precip_detr,3)
        corr_precip(i,j) = corr(SAM,precip_detr(:,i,j));
    end
end

corr_wind = zeros(size(wind_detr,2),size(wind_detr,3)); % Correlate SAM and wind
for i=1:size(wind_detr,2)
    for j=1:size(wind_detr,3)
        corr_wind(i,j) = corr(SAM,wind_detr(:,i,j));
    end
end

atcorr_precip = zeros(2,size(precip_detr,2),size(precip_detr,3)); % Autocorrelation of Precip
for i=1:size(precip_detr,2)
    for j=1:size(precip_detr,3)
        atcorr_precip(:,i,j) = autocorr(precip_detr(:,i,j),1);
    end
end

atcorr_wind = zeros(2,size(wind_detr,2),size(wind_detr,3)); % Autocorrelation of winds
for i=1:size(wind_detr,2)
    for j=1:size(wind_detr,3)
        atcorr_wind(:,i,j) = autocorr(wind_detr(:,i,j),1);
    end
end
atcorr_precip = squeeze(atcorr_precip(2,:,:));
atcorr_wind = squeeze(atcorr_wind(2,:,:));

sigma_precip = zeros(size(precip_detr,2),size(precip_detr,3)); % Calculate the standard deviation of the precip time-series
for i=1:size(precip_detr,2)
    for j=1:size(precip_detr,3)
        sigma_precip(i,j) = std(precip_detr(:,i,j));
    end
end

sigma_wind = zeros(size(wind_detr,2),size(wind_detr,3)); % Calculate the standard deviation of the wind time-series
for i=1:size(wind_detr,2)
    for j=1:size(wind_detr,3)
        sigma_wind(i,j) = std(wind_detr(:,i,j));
    end
end

%% Calculating the Synthetic Series
mkdir('Synth_Data')
for n=1:1000
    tic;
    eta_nu = randn(length(SAM),1);
    nu_wind=zeros(length(SAM),size(wind_detr,2),size(wind_detr,3),'single');
    nu_wind(1,:,:)=NaN;

    for i=1:size(wind_detr,2)
        for j=1:size(wind_detr,3)
            nu_wind(2:end,i,j) = wind_detr(i,j) + reg_wind_detr(i,j)*SAM(2:end) + ...
                sigma_wind(i,j)*sqrt(1.0-corr_wind(i,j)^2) * ...
                (eta_nu(2:end) + atcorr_wind(i,j)*eta_nu(1:end-1));
        end
    end

    nu_precip=zeros(length(SAM),size(precip_detr,2),size(precip_detr,3),'single');
    nu_precip(1,:,:)=NaN;

    for i=1:size(precip_detr,2)
        for j=1:size(precip_detr,3)
            nu_precip(2:end,i,j) = precip_detr(i,j) + reg_precip_detr(i,j)*SAM(2:end) + ...
                sigma_precip(i,j)*sqrt(1.0-corr_precip(i,j)^2) * ...
                (eta_nu(2:end) + atcorr_precip(i,j)*eta_nu(1:end-1));
        end
    end

   save(['Synth_Data/run',num2str(n),'syn.mat'],'nu_wind','nu_precip','eta_nu')
   toc;
end

%% Power Spectral Density

% rednoise=(eta_nu(2:end) + atcorr_ts(i,j)*eta_nu(1:end-1));
% Fs=length(rednoise);
% t = 0:1/Fs:1-1/Fs;
% 
% N = length(rednoise);
% xdft = fft(rednoise);
% xdft = xdft(1:N/2+1);
% psdx = (1/(Fs*N)).*abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:Fs/length(rednoise):Fs/2;
% plot(freq,10*log10(psdx)); grid on;
% title('Periodogram Using FFT');
% xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');