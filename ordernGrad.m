   function sum3 = ordernGrad(bnd_pnts,p2,pos,n1,adv,prob_int,coef,type)
    for i = 1:size(pos,1)
        sum1 = [0 0];
        sum2 = 0;
%         in1 = inhull(p2,vorvx{i},[],1e-15);
%         cl = find(in1);
%         q1 = p2(cl,:);
%         p_int1 = prob_int(cl,:);
        for l = 1:size(p2,1)
%             sum1 = sum1 + q1(l,:)*p_int1(l,:);  
%             sum2 = sum2 + p_int1(l,:); 
            
             prod = 1;
             for j=1:size(pos,1)
                 if i ~=j
                     prod = prod * (1-f_exp(p2(l,:),pos(j,:),coef)/f_exp([0 0],[0 0],coef));
                 else
                     prod = prod;
                 end
             end
            sum1 = sum1- (p2(l,:)-pos(i,:))* prod * ...
                 f_exp(p2(l,:),pos(i,:),coef) * prob_int(l,:) / ...
                 (coef*f_exp([0 0],[0 0],coef));            
             
                    
        end
%         if ~ismember(i,adv)
        sum4(i,:) = sum1/norm(sum1);
%         sum3(i,:) = pos(i,:) + sum1;
%         else
%             sum3(i,:) = pos(i,:);
%         end
    end
    
    
    
    
alph = 1;
for i =1:size(pos,1)
    pos_tmp(i,:) = pos(i,:) - alph*sum4(i,:);
end
% tl = randsample(1:size(pos,1),1);
% pos_tmp = pos;
% pos_tmp(tl,:) = pos(tl,:) - alph*sum4(tl,:);
% cst1 = calcCostExp(neib2,v1,pos,p2,coef,[],adv,type,prob_int);


cst1 = ordernCost(bnd_pnts,p2,pos,n1,coef,adv,type,prob_int);

% [~,v_1_r_1,~,~] = polybnd_voronoi(pos_tmp,bnd_pnts);
% [v_1_r,~,neib3] = polybnd_order2voronoi(pos_tmp,bnd_pnts);
% cst2 = calcCostExp(neib3,v_1_r,pos_tmp,p2,coef,[],adv,type,prob_int);


cst2 = ordernCost(bnd_pnts,p2,pos_tmp,n1,coef,adv,type,prob_int);


% cst1 = 0;
% cst2 = 1;
% max_Iter = 30;
k = 0;
while(cst2 > cst1)
    k = k+1;
    alph = alph * 0.9;
%     pos_tmp(tl,:) = pos(tl,:) - alph*sum4(tl,:);
    for i =1:size(pos,1)
        pos_tmp(i,:) = pos(i,:) - alph*sum4(i,:);
    end
%     [~,v_1_r_1,~,~] = polybnd_voronoi(pos_tmp,bnd_pnts);
%     [~,v_1_r_1,~,~] = polybnd_voronoi(pos_tmp,bnd_pnts);
%     cst2 = lloyd_cost_fin_exp(v_1_r_1,bnd_pnts,p2,pos_tmp,n1,coef,adv,type,prob_int);  
    cst2 = ordernCost(bnd_pnts,p2,pos_tmp,n1,coef,adv,type,prob_int);

end
for i = 1:size(pos,1)
    sum3(i,:) = pos_tmp(i,:);
end