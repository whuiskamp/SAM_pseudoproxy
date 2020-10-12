% Putting the data back together again
% Willem Huiskamp, 2015
% 
% If you messed up file names for line 15/26, use this in the terminal:
% rename -v 's/_5syncorr.mat/_05syncorr.mat/' *_5syncorr.mat


clear
for k = [61] % Window size
  sat_synthcorr = zeros(500,45,144);
  precip_synthcorr = zeros(500,45,144);
  for i=1:1000 
      i
    filenames = dir(['Synth_runcorr/',num2str(k),'yrWindow/run',num2str(i),'i_*']); % Looks into the directory and pulls a list of all filenames for a given window size and synthetic run
    for j = 1:9 
        load(['Synth_runcorr/',num2str(k),'yrWindow/',filenames(j).name]); % loads in the 9 chunks of data for a synth run
        j_start = (j-1)*5+1;                                          % establishing that we want to loop from 1 to 45 in steps of 5
        j_end = j*5;
        sat_synthcorr(:,j_start:j_end,:) = sat_synruncorr(:,j_start:j_end,:); % writes out data in the correct latitudes
        precip_synthcorr(:,j_start:j_end,:) = precip_synruncorr(:,j_start:j_end,:);
    end
    save(['Synth_runcorr/',num2str(k),'yrWindow/run',num2str(i),'corrfinal.mat'],'sat_synthcorr','precip_synthcorr') % Save reconstructed data to new directory
  end
end
