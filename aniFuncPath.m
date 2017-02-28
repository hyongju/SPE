%%% animation code 
%%% input p: cell array of position set
%%% stage x 1
%%% example: p{1} = [0.1 0.2;0.3 0.4;0.5 0.5]; read as at stage 1, agent 1
%%% is at (0.1, 0.2), agent 2 is at (0.3 0.4), and agent 3 is at (0.5 0.5)
%%% resolution - how many steps between each move e.g. 30
%%% adv_idx - the heterogeneous agents' index. e.g. given p{1} above, if
%%% agent 1, and 3 is the hetero-agents, adv_idx is a row vector [1 3]
% clear all;close all;clc
function [path,thOut] = aniFuncPath(thIn,p,resolution,k_p)
% load('psav_test1.mat');
% resolution = 120;
% n = size(p_sav{1},1);
% p = p_sav';
% k_p = 5;
th{1}{1} = thIn;
% th{1}{1}=2*pi*rand(size(p{1},1),1);
% f=0;
% pos_init = [23.5425   25.9122];
% pos_init2 = repmat(pos_init,size(p{1},1),1);
% color_code = jet(100);
% adv_idx = 1;
N = 2;
% N = 5;

adv_idx = [];
% n_adv = f;
vmax = -0.03;
n_adv = size(adv_idx,2);
n_cop = size(p{1},1) - n_adv;
cop_idx= setdiff(1:size(p{1},1),adv_idx);
for u1 = 1:N
    pos_int{u1} = mean(p{u1}(cop_idx,:),1);
end
delta_p = cell(N,1);
p_appended = cell(N,1);

for i = 1:N-1
    i
%    delta_p{i} = p{i+1} - p{i};
   if i >= 2
       th{i}{1} = th{i-1}{resolution};
   end
   for j = 1:resolution
     if j >= 2
         delta_p1{i}{j} = p_appended{i}{j-1} - p{i+1};
         v{i}{j} = -k_p* dota([cos(th{i}{j-1}) sin(th{i}{j-1})],delta_p1{i}{j});
         w{i}{j} = 2*k_p * arcta(th{i}{j-1},delta_p1{i}{j});
         p_appended{i}{j} = p_appended{i}{j-1} + [v{i}{j} .* cos(th{i}{j-1}) v{i}{j} .*sin(th{i}{j-1})] / resolution;
         th{i}{j} = th{i}{j-1} + 1/resolution * w{i}{j};
     else
         p_appended{i}{j} = p{i};
     end
   end
end
thOut = th{i}{j};
k = 0; k2 = 0;
for i = 1:N-1
    for j = 1:resolution
        k = k + 1;
        path{k} = p_appended{i}{j};
    end
end
