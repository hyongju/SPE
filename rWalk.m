% random walk...
% clear all;close all;clc
% 
% boundary 

function out = rWalk(pos,u)

xLim = [0 1];
yLim = [0 1];

% pos = rand(10,2);


% u = 0.2;

p{1} = pos;
for i = 2:2
    for j = 1:size(pos,1)
        rng('shuffle')
        phi = rand(1) * 2*pi;
        while ~(p{i-1}(j,1) + u*cos(phi) >=0 && p{i-1}(j,1) + u*cos(phi) <= 1 ...,
                && p{i-1}(j,2) + u*sin(phi) >= 0 && p{i-1}(j,2) + u*sin(phi) <= 1)
    %     rnd('shuffle')
            phi = rand(1) * 2*pi;
        end
        p{i}(j,:) = p{i-1}(j,:) + u*[cos(phi) sin(phi)];
    end
end


% 
% figure,
% for i = 1:size(pos,1)
%     aug_p{i} = []
% end
% for j = 1:size(pos,1)
%     for i = 1:2
%         aug_p{j} = [aug_p{j};p{i}(j,:)];
%     end
% end
% for j = 1:size(pos,1)
%     plot(aug_p{j}(:,1),aug_p{j}(:,2),'-.');hold on;
% end
% axis('equal');
% axis([0 1 0 1])


out = p{2};
