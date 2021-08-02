% r val for noise proxies. i = 101:110
% matlab -nosplash -nodesktop -r "r_noise_calc_i101_110; exit"
% Run as a job - this takes a LONG time

clear all
land = ncread('DataFiles/sftlf_A1.static.nc','sftlf')';
load('DataFiles/model_output.mat')
irange = 1:144;
jrange = 45;
nsim = 1000;

tic
for windowsize = 61
    for i= irange
        for j= 1:jrange
            if ~isnan(land(j,i)) && ~exist(['noise_proxies/noise_rvals',num2str(windowsize),'yrWDW_i',num2str(i),'_j',num2str(j),'.mat'],'file')
            load(['noise_proxies/noise_prx_i',num2str(i),'_j',num2str(j),'.mat'])
                for k = 1:nsim
                    r_sat(:,k) = movingCorrelation([squeeze(noise_sat(:,k)),SAM],windowsize,2);
                    r_pre(:,k) = movingCorrelation([squeeze(noise_pre(:,k)),SAM],windowsize,2);
                end
                save(['noise_proxies/noise_rvals',num2str(windowsize),'yrWDW_i',num2str(i),'_j',num2str(j),'.mat'],'r_sat','r_pre')
                clear r_sat r_pre
            end
        end
        i
    end
    windowsize
end
toc
