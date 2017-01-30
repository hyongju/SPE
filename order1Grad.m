   function sum3 = order1Grad(vorvx,bnd_pnts,p2,pos,n1,adv,prob_int,coef,type)
    for i = 1:size(pos,1)
        sum1 = [0 0];
        sum2 = 0;
        in1 = inhull(p2,vorvx{i},[],1e-15);
        cl = find(in1);
        q1 = p2(cl,:);
        p_int1 = prob_int(cl,:);
        for l = 1:size(q1,1)
            sum1 = sum1- (q1(l,:)-pos(i,:))* ...
                 f_exp(q1(l,:),pos(i,:),coef) * p_int1(l,:) / ...
                 (coef*f_exp([0 0],[0 0],coef));
                    
                    
        end

        sum4(i,:) = sum1;

    end
    
    
    
    
alph = 1;
for i =1:size(pos,1)
    pos_tmp(i,:) = pos(i,:) - alph*sum4(i,:);
end

cst1 = order1Cost(vorvx,bnd_pnts,p2,pos,n1,coef,adv,type,prob_int);

[~,v_1_r_1,~,~] = pVoronoi(pos_tmp,bnd_pnts);

cst2 = order1Cost(v_1_r_1,bnd_pnts,p2,pos_tmp,n1,coef,adv,type,prob_int);

k = 0;
while(cst2 > cst1)
    k = k+1;
    alph = alph * 0.9;
    for i =1:size(pos,1)
        pos_tmp(i,:) = pos(i,:) - alph*sum4(i,:);
    end
    [~,v_1_r_1,~,~] = pVoronoi(pos_tmp,bnd_pnts);
    cst2 = order1Cost(v_1_r_1,bnd_pnts,p2,pos_tmp,n1,coef,adv,type,prob_int);   
end
for i = 1:size(pos,1)
    sum3(i,:) = pos_tmp(i,:);
end