function [cst,out] = ctrlOrder2(pos, p2,p_int,coef)

vmax = 1/eps;       % velocity constraint is RELAXED....
stage = 30;

p2_0 = [0 0;1 0;1 1;0 1];
bnd_idx = convhull(p2_0);
bnd_pnts = p2_0(bnd_idx,:);


n1 = length(p_int);
%%

p_sav{1} = pos;
adv = [];           % index set of faulty nodes
type = 3;
% type 1 - sensor failure
% type 2 - servo failure
% type 3 - sensor/servo failure

p_sav{1} = pos;

%% call function
for t = 1:stage
    t
    [voronoi_rg{t},neib1{t},neib2{t}] = p2Voronoi(pos,bnd_pnts);
    
    %   order2 Voronoi regions - {cell} data structure
    %   neib1: neighbors indice can be repeated
    %   neib2: neighbors indice cannot be repeated

    active{t}  = chooseSset(t,size(pos,1),neib1{t});
    l_min{t} = order2Grad(active{t},neib1{t},voronoi_rg{t},pos,p2,size(pos,1),coef,p_int,adv,type,bnd_pnts,neib2{t});   % l_min: local minimizer
    [cst(t),~] = order2Cost(neib2{t},voronoi_rg{t},pos,p2,coef,n1,adv,type,p_int);
    idx{t} = find(active{t});
    if type == 2 || type == 3
        for y = 1:length(idx{t})
            if ~ismember(idx{t}(y),adv)
                if norm(l_min{t}{idx{t}(y)}- pos(idx{t}(y),:)) <= vmax
                    pos(idx{t}(y),:) = l_min{t}{idx{t}(y)};
                else
                    pos(idx{t}(y),:) = pos(idx{t}(y),:) + vmax * (l_min{t}{idx{t}(y)}- pos(idx{t}(y),:))/norm(l_min{t}{idx{t}(y)}- pos(idx{t}(y),:));
                end
            end
        end
    else
        for y = 1:length(idx{t})
            if norm(l_min{t}{idx{t}(y)}- pos(idx{t}(y),:)) <= vmax
                pos(idx{t}(y),:) = l_min{t}{idx{t}(y)};
            else
                pos(idx{t}(y),:) =  pos(idx{t}(y),:) + vmax* (l_min{t}{idx{t}(y)}- pos(idx{t}(y),:))/norm(l_min{t}{idx{t}(y)}- pos(idx{t}(y),:));
            end
       end
    end
    p_sav{t+1} = pos;
end

out = p_sav{t+1};

