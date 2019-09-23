%% This script calculates the EOF of running correlations between SAT/Precip
% and SAM over individual countries/ regions rather than the whole SH.
% This code is based on the key_regions script.
%
% August 2019

% Setup
clear
load('DataFiles/model_output.mat','lat');
land = ncread('sftlf_A1.static.nc','sftlf')'; 
lat_2 = 45; % change this from 45 to get some sub-set of the SH
weight = zeros(lat_2,size(land,2));
for i=1:size(land,2)
    weight(:,i) = sqrt(cos(lat(1:lat_2)*3.141/180));
end

%Define regions
Gl_lon = 1:144; Gl_lat = 1:45;     % S. Hemisphere
SA_lon = 110:130; SA_lat = 15:45;  % South America
Au_lon = 40:75; Au_lat = 15:41;    % Australia and NZ
AA_lon = 1:144; AA_lat = 1:15;     % Antarctica
SoA_lon = 1:25; SoA_lat = 25:45;   % South Africa

for windowsize = [31 61 91]
    load(['DataFiles/runcorr',num2str(windowsize),'yrwdw.mat']);
    for i = 1:size(sat_runcorr,1)
        sat_runcorr(i,isnan(land)) = nan;
        precip_runcorr(i,isnan(land)) = nan;
    end
    sat_runcorr = sat_runcorr(ceil(windowsize/2):500-floor(windowsize/2),1:lat_2,:); % we only need the SH and eliminate the NaNs from the running correlation at either end.
    pre_runcorr = precip_runcorr(ceil(windowsize/2):500-floor(windowsize/2),1:lat_2,:);
    
    % Apply latitude weighting here
    sat_runcorr_w = zeros(size(sat_runcorr)); pre_runcorr_w = zeros(size(pre_runcorr));
    for i = 1:size(sat_runcorr,1)
        sat_runcorr_w(i,:,:) = squeeze(sat_runcorr(i,:,:)) .* weight;
        pre_runcorr_w(i,:,:) = squeeze(pre_runcorr(i,:,:)) .* weight;
    end
    
    % test = squeeze(sat_runcorr_w(1,:,:));
    
    for region = 1:5
        if region == 1
            lat = Gl_lat; lon = Gl_lon; reg_name = 'SH';
        elseif region == 2
            lat = SA_lat; lon = SA_lon; reg_name = 'SA';
        elseif region == 3
            lat = Au_lat; lon = Au_lon; reg_name = 'Au';
        elseif region == 4
            lat = AA_lat; lon = AA_lon; reg_name = 'AA';
        elseif region == 5
            lat = SoA_lat; lon = SoA_lon; reg_name = 'SoA';
        end
        
        sat_rc_reg = sat_runcorr_w(:,lat,lon);
        pre_rc_reg = pre_runcorr_w(:,lat,lon);
        % Create index arrays
        I = nan(size(lat,2),size(lon,2)); J = nan(size(lat,2),size(lon,2)); 
        for i = 1:size(lat,2)
            I(i,:) = i;
        end
        for j = 1:size(lon,2)
            J(:,j) = j;
        end
            
        % Reshape so that we have a 2D matrix of time X cells
        sat_rc_2 = reshape(sat_rc_reg,size(sat_rc_reg,1),size(sat_rc_reg,2)*size(sat_rc_reg,3)); 
        pre_rc_2 = reshape(pre_rc_reg,size(pre_rc_reg,1),size(pre_rc_reg,2)*size(pre_rc_reg,3));
        % Do the same with our indices
        I_2 = reshape(I,1,size(sat_rc_reg,2)*size(sat_rc_reg,3)); J_2 = reshape(J,1,size(sat_rc_reg,2)*size(sat_rc_reg,3));
        
        % Create anomalies
        [n,p] = size(sat_rc_2); [q,w] = size(pre_rc_2);
        Xbar = nanmean(sat_rc_2,1); Xbar_p = nanmean(pre_rc_2,1);
        sat_anom = sat_rc_2-ones(n,1)*Xbar; pre_anom = pre_rc_2-ones(q,1)*Xbar_p;
                    
        % Remove non-land cells
        [~, K] = find(isnan(sat_anom(1,:))); % We need to eliminate cells (cols) with NaN values. We only need to look at one year, as mask never changes.
        I_2(K) = []; J_2(K) = [];  % Remove these cells from our index, so we can later put the 2D array back together
        sat_anom(:,K) = []; pre_anom(:,K) = []; % Now remove these cells from the array itself so we can do the PCA
        
        % Calculate EOFs
        [sat_EOFs,sat_PC,~,~,sat_var_exp] = pca(sat_anom);
        [pre_EOFs,pre_PC,~,~,pre_var_exp] = pca(pre_anom);
        
        % To reshape back into lat/lon space, we need our index
        sat_EOF_reg = nan(size(sat_EOFs,2),size(sat_rc_reg,2),size(sat_rc_reg,3));
        pre_EOF_reg = nan(size(pre_EOFs,2),size(pre_rc_reg,2),size(pre_rc_reg,3));
        
        for i = 1:size(sat_EOFs,2)
            for j = 1:size(I_2,2)
                sat_EOF_reg(i,I_2(1,j),J_2(1,j)) = sat_EOFs(j,i);
                pre_EOF_reg(i,I_2(1,j),J_2(1,j)) = pre_EOFs(j,i);
            end
        end
        % test = squeeze(sat_EOF_reg(1,:,:));
        save(['DataFiles/reg_EOFs_',num2str(windowsize),'yrwdw_',reg_name,'.mat'],'sat_EOF_reg','sat_PC','sat_var_exp','pre_EOF_reg','pre_PC','pre_var_exp')
        clear lat lon sat_rc_reg pre_rc_reg sat_rc_2 pre_rc_2 I J K I_2 J_2 n p q w Xbar Xbar_p sat_anom pre_anom sat_EOFs ...
            sat_PC sat_var_exp sat_EOF_reg pre_EOFs pre_PC pre_var_exp pre_EOF_reg
        reg_name
    end
    windowsize    
end
    
%% Post-processing/ plotting 
clear
load(['DataFiles/model_output.mat'],'sat_detr','precip_detr');
%Define regions
Gl_lon = 1:144; Gl_lat = 1:45;     % S. Hemisphere
SA_lon = 110:130; SA_lat = 15:45;  % South America
Au_lon = 40:75; Au_lat = 15:41;    % Australia and NZ
AA_lon = 1:144; AA_lat = 1:15;     % Antarctica
SoA_lon = 1:25; SoA_lat = 25:45;   % South Africa
levels = -0.1:0.01:0.1;
levels_p = 0.005:0.02:0.1;
levels_n = -0.1:0.02:0;

for region = 1:5
        if region == 1
            reg_lat = Gl_lat; reg_lon = Gl_lon; reg_name = 'SH';
        elseif region == 2
            reg_lat = SA_lat; reg_lon = SA_lon; reg_name = 'SA';
        elseif region == 3
            reg_lat = Au_lat; reg_lon = Au_lon; reg_name = 'Au';
        elseif region == 4
            reg_lat = AA_lat; reg_lon = AA_lon; reg_name = 'AA';
        elseif region == 5
            reg_lat = SoA_lat; reg_lon = SoA_lon; reg_name = 'SoA';
        end
        for w = [31 61 91]
            % Standardise the PCs and find their standard deviation
            load(['DataFiles/reg_EOFs_',num2str(w),'yrwdw_',reg_name,'.mat'],'sat_EOF_reg','sat_PC','sat_var_exp','pre_EOF_reg','pre_PC','pre_var_exp')
            for i = 1:(size(sat_PC,2))
                std_sat_PC(floor(w/30),1,i) = std(sat_PC(:,i));
                std_pre_PC(floor(w/30),1,i) = std(pre_PC(:,i));
                norm_sat_PC(:,i) = sat_PC(:,i)/squeeze(std_sat_PC(floor(w/30),1,i)); % Standardise by / by std
                norm_pre_PC(:,i) = pre_PC(:,i)/squeeze(std_pre_PC(floor(w/30),1,i));
                sat_EOFs(i,:,:) = sat_EOF_reg(i,:,:)*squeeze(std_sat_PC(floor(w/30),1,i)); % In order to be able to re-make the original data, we have to multiply the EOFs by the std
                pre_EOFs(i,:,:) = pre_EOF_reg(i,:,:)*squeeze(std_pre_PC(floor(w/30),1,i));
            end
            % Plot the PCs
            figure(floor(w/30))
            for i = 1:3
                subplot(3,1,i)
                plot(norm_sat_PC(:,i),'color','black','linewidth',2)
                title([reg_name,' PC',num2str(i),'for ',num2str(w),'yrWdw: SAT - Black, Precip - Red'])
                hold
                plot(norm_pre_PC(:,i),'color','red','linewidth',2)
            end
            h=gcf;
            set(h,'PaperPositionMode','auto');         
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);
        
            %print(gcf,'-dpdf','-painters',['figures/3PCs_runncorr_',num2str(w),'yrWdW_',reg_name,'.pdf'])
            close(floor(w/30))
            
            % This plots the first 3 EOFs of the running correlation between SAT and SAM
            subplot1(3, 1, 'Gap', [.01 .03], 'XTickL', 'Margin', 'YTickL', 'Margin');
            for i = 1:3
                subplot1(i)
                axesm('MapProjection','eqdcylin','origin',[0,180],'MapLatLimit',[-90 0],'meridianlabel','on','parallellabel','on','PLabelMeridian',0,'Mlabelparallel','south')
                framem
                gridm %If you want gridlinesload coast
                load('DataFiles/model_output.mat','lat','lon')
                contourfm(lat(reg_lat),lon(reg_lon),squeeze(sat_EOFs(i,:,:)),levels,'linestyle','none') 
                contourm(lat(reg_lat),lon(reg_lon),squeeze(sat_EOFs(i,:,:)),levels_p,'k'); contourm(lat(reg_lat),lon(reg_lon),squeeze(sat_EOFs(i,:,:)),levels_n,'k','linestyle','--');
                colormap(b2r(-0.1,0.1))
                load coast
                plotm(lat,long,'k','linewidth',2)
                title(['SAT EOF',num2str(i),' for',num2str(w),'yr Window - ',num2str(round(sat_var_exp(i,:))),'% var. exp.'])
                clear lat long
            end
            h=gcf;
            set(h,'PaperPositionMode','auto');         
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);
            print(gcf,'-dpdf','-painters',['figures/',reg_name,'_SAT_EOFs_',num2str(w),'yrs'])
            close(1)
            
            % This plots the first 3 EOFs of the running correlation between Precip and SAM
            subplot1(3, 1, 'Gap', [.01 .03], 'XTickL', 'Margin', 'YTickL', 'Margin');
            for i = 1:3
                subplot1(i)
                axesm('MapProjection','eqdcylin','origin',[0,180],'MapLatLimit',[-90 0],'meridianlabel','on','parallellabel','on','PLabelMeridian',0,'Mlabelparallel','south')
                framem
                gridm %If you want gridlinesload coast
                load('DataFiles/model_output.mat','lat','lon')
                contourfm(lat(reg_lat),lon(reg_lon),squeeze(pre_EOFs(i,:,:)),levels,'linestyle','none') 
                contourm(lat(reg_lat),lon(reg_lon),squeeze(pre_EOFs(i,:,:)),levels_p,'k'); contourm(lat(reg_lat),lon(reg_lon),squeeze(pre_EOFs(i,:,:)),levels_n,'k','linestyle','--');
                colormap(b2r(-0.1,0.1))
                load coast
                plotm(lat,long,'k','linewidth',2)
                title(['Precip EOF',num2str(i),' for',num2str(w),'yr Window - ',num2str(round(pre_var_exp(i,:))),'% var. exp.'])
                clear lat long
            end
            h=gcf;
            set(h,'PaperPositionMode','auto');         
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);
            print(gcf,'-dpdf','-painters',['figures/',reg_name,'_Pre_EOFs_',num2str(w),'yrs'])
            close(1)
            
            % Regress PCs against SAT/Precip fields
            N = 3; % select number of PCs to ivestigate
            
            sat_detr_2 = sat_detr(ceil(w/2):500-floor(w/2),:,:);
            pre_detr_2 = precip_detr(ceil(w/2):500-floor(w/2),:,:);
            for n = 1:N
                for i = 1:90
                    for j = 1:144
                        SAM_SAT_regr(n,i,j) = regress(norm_sat_PC(:,n),squeeze(sat_detr_2(:,i,j)));
                        SAM_pre_regr(n,i,j) = regress(norm_pre_PC(:,n),squeeze(pre_detr_2(:,i,j)));
                    end
                end
            end
            %save(['DataFiles/reg_EOFs_',num2str(w),'yrwdw_',reg_name,'.mat'],'SAM_SAT_regr','SAM_precip_regr','-append')
            
            
            
        clear norm_sat_PC norm_pre_PC pre_EOFs sat_EOFs std_sat_PC std_pre_PC sat_detr_2 pre_detr_2
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        end
end        
