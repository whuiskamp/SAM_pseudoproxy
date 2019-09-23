% Putting the data back together again
% make two different directories, as these files have the same name.
% send one loop to winds, the other to precip and output in the final,
% original directory

wind_synthcorr = zeros(500,45,144);
for i=1:1000
   filenames = dir(['run',num2str(i),'i_*']);
   for j = 1:9
       load(filenames(j).name);
       j_start = (j-1)*5+1;
       j_end = j*5;
       wind_synthcorr(:,j_start:j_end,:) = wind_synruncorr(:,j_start:j_end,:);
   end
   save(['Synth_runcorr/31yrWindow/run',num2str(i),'corrfinal.mat'],'wind_synthcorr')
end

precip_synthcorr = zeros(500,45,144);
for i=1:1000
   filenames = dir(['run',num2str(i),'i_*']);
   for j = 1:9
       load(filenames(j).name);
       j_start = (j-1)*5+1;
       j_end = j*5;
       precip_synthcorr(:,j_start:j_end,:) = precip_synruncorr(:,j_start:j_end,:);
   end
   save(['Synth_runcorr/31yrWindow/run',num2str(i),'corrfinal.mat'],'precip_synthcorr','-append')
end

