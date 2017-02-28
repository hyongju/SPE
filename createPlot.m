clear all;close all;clc
load('ordern_result_fault.mat')

pos = savData{2,1};
particleWgt = savData{2,2};
normExpWgt = gridDisplay([sampPos particleWgt'],nGrid,max(particleWgt));
fig4 = figure('position',[100 100 600 600],'Color',[1 1 1]);

imshow(normExpWgt','InitialMagnification','fit','Colormap',jet(255));hold on;
plot(pos(:,1)*nGrid,pos(:,2)*nGrid,'p','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',16);hold on;
        plot(pos(3,1)*nGrid,pos(3,2)*nGrid,'-o',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerSize',16); hold on;
        plot(pos(8,1)*nGrid,pos(8,2)*nGrid,'-o',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerSize',16);hold on;


set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlabel([]);
ylabel([]);
axis([0 1*nGrid 0 1*nGrid]);
set(gca,'FontSize',16);
set(gca,'Ydir','reverse');



pos = savData{11,1};
particleWgt = savData{11,2};
normExpWgt = gridDisplay([sampPos particleWgt'],nGrid,max(particleWgt));
fig4 = figure('position',[100 100 600 600],'Color',[1 1 1]);

imshow(normExpWgt','InitialMagnification','fit','Colormap',jet(255));hold on;
plot(pos(:,1)*nGrid,pos(:,2)*nGrid,'p','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',16);hold on;
        plot(pos(3,1)*nGrid,pos(3,2)*nGrid,'-o',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerSize',16);hold on
        plot(pos(8,1)*nGrid,pos(8,2)*nGrid,'-o',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerSize',16);hold on;


set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlabel([]);
ylabel([]);
axis([0 1*nGrid 0 1*nGrid]);
set(gca,'FontSize',16);
set(gca,'Ydir','reverse');


% normExpWgt = gridDisplay([sampPos particleWgt'],nGrid,particleWgt(601)/0.3);



