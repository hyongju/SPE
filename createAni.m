clear all;close all;clc
load('most_recent_data2.mat')


% normExpWgt = gridDisplay([sampPos particleWgt'],nGrid,particleWgt(601)/0.3);

fig4 = figure('position',[100 100 600 600],'Color',[1 1 1]);

% for j = 1:2

for j = 1:size(savData,1)
    j
    particleWgt = savData{j,2};
    normExpWgt = gridDisplay([sampPos particleWgt'],nGrid,particleWgt(601)/0.3);
    pos = savData{j,1};
    imshow(normExpWgt','InitialMagnification','fit');hold on;
    plot(pos(:,1)*nGrid,pos(:,2)*nGrid,'p','MarkerFaceColor','w','MarkerEdgeColor','c','MarkerSize',10);hold on;
    for i = 1:size(b_Poly,1)
        plot(b_Poly{i}(:,1)*nGrid,b_Poly{i}(:,2)*nGrid,'b-');hold on;
    end
    convMinTmp = floor((j-1)/60);
    if convMinTmp < 14
        convMin = 46+convMinTmp;
        convHur = 21;
    else
        convMin = convMinTmp - 14;
        convHur = 22;
    end
    convSecTmp = mod((j-1),60);
    if numel( num2str(convSecTmp)) == 1
        convSec = ['0' num2str(convSecTmp)];
    else
        convSec = num2str(convSecTmp);
    end 
    if numel( num2str(convMin)) == 1
        convMin = ['0' num2str(convMin)];
    else
        convMin = num2str(convMin);
    end  
    strng = ['11-Aug-2014 ' num2str(convHur) ':' num2str(convMin) ':' num2str(convSec)];
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'FontSize',16);
    set(gca,'Ydir','normal');
%     xlabel([]);
    xlabel(strng)    
    ylabel([]);
    M(j) = getframe(fig4);
end

fig5 = figure('position',[100 100 600 600],'Color',[1 1 1]);
% axis('square')
% axis([0 1 0 1 0 1])
% axis('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
movie(fig5,M,1,30);
v = VideoWriter('test_prep.mp4','MPEG-4');
open(v)
writeVideo(v,M)
close(v)

