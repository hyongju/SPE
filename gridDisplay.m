% grid_display
function at_ij_norm = gridDisplay(sample,num_grid,mpint)
% input: sample n x (d+1) matrix
%        num_grid N
% output: NxN matrix with intensity

% clear all;close all;clc
% clear all;close all;clc
% sample = rand(100,3);
% 
% 
% num_grid = 1000;
grid_size = 1/num_grid;
for i = 1:num_grid
    for l = 1:num_grid
        k=0;
        for j = 1:size(sample,1)
        % at grid i,l
            if sample(j,1) >= grid_size*(i-1) && sample(j,1) <= grid_size*(i) && sample(j,2) >= grid_size*(l-1) && sample(j,2) <= grid_size*(l)
                k = k +1;
                at_il{i,l}(k)=sample(j,3);
            end
        end
    end
end
for i = 1:size(at_il,1)
    for j = 1:size(at_il,2)
        if isempty(at_il{i,j})
            at_ij2(i,j) = 0;
        else
            at_ij2(i,j) = mean(at_il{i,j});
        end
    end
end
% [row,col]=find(at_ij2);
% for i = 1:size(row,1)
%     at_tot(i) =at_ij2(row(i),col(i));
% end
% at_ij_norm = (at_ij2-min(at_tot))/(max(max(at_ij2))-min(at_tot));
at_ij_norm = at_ij2/mpint;