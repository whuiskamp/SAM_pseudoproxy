%% Extra figures and analysis for GISS data
% 
%% 1. Simple correlation plot

load('/media/huiskamp/My Book/CM2.1/SAM/GISS_data/model_output.mat')
levels = -1:0.1:1;
levels_n = -1:0.1:-0.3;
levels_p = 0.3:0.1:1;

for i = 1:size(satG_detr,2)
    for j = 1:size(satG_detr,3)
        GISS_SAT_corr(i,j) = double(corr(SAM_G,squeeze(satG_detr(:,i,j))));
        GISS_pre_corr(i,j) = double(corr(SAM_G,squeeze(precipG_detr(:,i,j))));
    end
end

save('DataFiles/GISS_corrs.mat','GISS_SAT_corr','GISS_pre_corr')

for i = 1:2
    figure(i)
    axesm('MapProjection','pcarree','MapLatLimit',[-90 -0],'maplonlim',[0 360])
    framem;
    gridm;
    mlabel;
    plabel;
    tightmap
end

% model SAT
figure(1)
contourfm(latG(1:45),lonG,model_SAT_corr(1:45,:),levels,'linestyle','none')
hold on;
contourm(latG(1:45),lonG,model_SAT_corr(1:45,:),levels_n,'color','k')
contourm(latG(1:45),lonG,model_SAT_corr(1:45,:),levels_p,'color','k')
colormap(b2r(-1,1))

% model Precip
figure(2)
contourfm(latG(1:45),lonG,model_pre_corr(1:45,:),levels,'linestyle','none')
hold on;
contourm(latG(1:45),lonG,model_pre_corr(1:45,:),levels_n,'color','k')
contourm(latG(1:45),lonG,model_pre_corr(1:45,:),levels_p,'color','k')
colormap(b2r(-1,1))

for i = 1:4
    figure(i)
    load coast
    plotm(lat,long,'k')
end
figure(1)
print(gcf,'figures/GISS_SAM_SAT_corr','-dpdf','-bestfit')

figure(2)
print(gcf,'figures/GISS_SAM_pre_corr','-dpdf','-bestfit')






