clear all;close all;clc
load('order2_result_fault.mat')


% normExpWgt = gridDisplay([sampPos particleWgt'],nGrid,particleWgt(601)/0.3);

fig4 = figure('position',[100 100 600 600],'Color',[1 1 1]);
% for j = 1:2
    k_p = 5;
    res = 60;
    thIn = 2*pi*rand(size(pos,1),1);
    k = 0;
for j = 1:size(savData,1)-1
    j
    particleWgt = savData{j,2};

    pos = savData{j,1};
%     plot(pos(:,1)*nGrid,pos(:,2)*nGrid,'p','MarkerFaceColor','w','MarkerEdgeColor','c','MarkerSize',10);hold on;
    posAug{1} = pos;
    posAug{2} = savData{j+1,1};
    if j == 1
        normExpWgt = gridDisplay([sampPos 0.3*particleWgt'],nGrid,max(particleWgt));            
    else
        normExpWgt = gridDisplay([sampPos particleWgt'],nGrid,max(particleWgt));    
    end
    imshow(normExpWgt','InitialMagnification','fit','Colormap',jet(255));hold on;
    [path,thOut] = aniFuncPath(thIn,posAug,res,k_p);
    for i = 1:length(path)
        imshow(normExpWgt','InitialMagnification','fit','Colormap',jet(255));hold on;
        plot(path{i}(:,1)*nGrid,path{i}(:,2)*nGrid,'p','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',16);hold on;
        plot(path{i}(3,1)*nGrid,path{i}(3,2)*nGrid,'-o',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerSize',16)
        plot(path{i}(8,1)*nGrid,path{i}(8,2)*nGrid,'-o',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerSize',16);hold on;
        hold off;
        axis([0 1*nGrid 0 1*nGrid]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'FontSize',16);
        set(gca,'Ydir','reverse');
        if numel( num2str(i)) == 1
            convSec = ['0' num2str(i)];
        else
            convSec = num2str(i);
        end  
        if j == size(savData,1)-1 && i ==  length(path)
            strng = ['Time Step: ' num2str(j) ' [' convSec '/60]'];
            else
            strng = ['Time Step: ' num2str(j-1) ' [' convSec '/60]'];
        end
        xlabel(strng);
        ylabel([]);
        k = k+1
        M(k) = getframe(fig4);
    end
    hold on;
    thIn = thOut;


end

fig5 = figure('position',[100 100 600 600],'Color',[1 1 1]);
% axis('square')
% axis([0 1 0 1 0 1])
% axis('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
movie(fig5,M,1,30);
v = VideoWriter('test_flt2_color.mp4','MPEG-4');
open(v)
writeVideo(v,M)
close(v)

