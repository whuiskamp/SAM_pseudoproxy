%% What's going on with the 31yr precip proxies

clear
load('DataFiles/CM2_31yr_runcorr.mat','pre_CM_runcorr')
load('DataFiles/nonstat_map_31yrwdw.mat','nonstat_precipmap')

pre_CM_runcorr(isnan(pre_CM_runcorr)) = 0;
for i = 1:144
    for j = 1:90
        if nonstat_precipmap(j,i) < 50
            pre_CM_runcorr(:,i,j) = 0;
        end
        if all(abs(pre_CM_runcorr(:,i,j))<0.3)
            pre_CM_runcorr(:,i,j) = 0;
        end
    end
end

figure(1); hold on
count = 0;
for i=71:144
    for j = 1:12
        if mean(pre_CM_runcorr(:,i,j)) ~= 0
            plot(pre_CM_runcorr(:,i,j))
            count = count + 1;
        end
    end
end


test = squeeze(pre_CM_runcorr(:,5,3));
test2 = squeeze(pre_CM_runcorr(300,:,:));





















