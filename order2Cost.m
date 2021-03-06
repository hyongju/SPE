function [sum2,indx] = order2Cost(neib4,voronoi_rg,pos,p2,coef,n1,adv,type,prob_int)
% input: active, neib3, voronoi_rg, pos,p2
% output: sum2
sum2 = 0;           % going to be the total cost...
indx = 0;    
cnt = 0;
for i = 1:size(neib4,2)             % for each i 
    for j = 1:length(neib4{i})      % for each j associated with i
%         clear clcv;
        clcv = voronoi_rg{i,neib4{i}(j)};     % for each order-2 voronoi region guarted by (i,j)  
        if ~isempty(clcv)                     % check if the region is empty. If not, proceed...
            in1 = inhull(p2,clcv,[],1e-15);   % obtain all sample points inside the region
            cl = find(in1);                   % find the index
            q1 = p2(cl,:);                    % return the positions of the sample points  
            p_int1 = prob_int(cl,:);
            if type == 1 || type == 3         % 
                if ~ismember(i,adv) && ~ismember(neib4{i}(j),adv)
                    for l = 1:size(q1,1);
                        sum2 = sum2 + (1-f_exp(q1(l,:),pos(neib4{i}(j),:),coef)/f_exp([0 0],[0 0],coef)) * (1-f_exp(q1(l,:),pos(i,:),coef)/f_exp([0 0],[0 0],coef)) *p_int1(l,:);
                        indx(cl(l)) = (1-f_exp(q1(l,:),pos(neib4{i}(j),:),coef)) * (1-f_exp(q1(l,:),pos(i,:),coef)) *p_int1(l,:);
                    end
                elseif ismember(i,adv) && ~ismember(neib4{i}(j),adv)
                    for l = 1:size(q1,1);                
                        sum2 = sum2 + (1-f_exp(q1(l,:),pos(neib4{i}(j),:),coef))*p_int1(l,:);
                        indx(cl(l)) = (1-f_exp(q1(l,:),pos(neib4{i}(j),:),coef))*p_int1(l,:);
                    end
                elseif ~ismember(i,adv) && ismember(neib4{i}(j),adv)
                    for l = 1:size(q1,1);                
                        sum2 = sum2 + (1-f_exp(q1(l,:),pos(i,:),coef))*p_int1(l,:);
                        indx(cl(l)) = (1-f_exp(q1(l,:),pos(i,:),coef))*p_int1(l,:);
                    end
                else
                    for l = 1:size(q1,1);                
                        sum2 = sum2 + p_int1(l,:);
                        indx(cl(l)) = p_int1(l,:);
                    end
                end
            else
                for l = 1:size(q1,1);           % for all sample points...
                    sum2 = sum2 + (1-f_exp(q1(l,:),pos(neib4{i}(j),:),coef)) * (1-f_exp(q1(l,:),pos(i,:),coef))*p_int1(l,:);
                    indx(cl(l)) = (1-f_exp(q1(l,:),pos(neib4{i}(j),:),coef)) * (1-f_exp(q1(l,:),pos(i,:),coef))*p_int1(l,:);
                    cnt = cnt + 1;
                end
            end
        end
    end
end
% cnt
