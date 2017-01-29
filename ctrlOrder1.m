function [cst3,out] = ctrlOrder1(pos, p2,p_int,coef)
stage = 30;
bnd_pnts = [0 1;1 1;1 0;0 0];
%%
p_sav{1} = pos;
adv = [];
type = 3;
%% call function
for t = 1:stage
    t
    [~,vorvx{t},~,~] = pVoronoi(pos,bnd_pnts);
    sum3{t} = order1Grad(vorvx{t},bnd_pnts,p2,pos,length(p_int),adv,p_int,coef,type);
    [cst3(t),indx{t}]= order1Cost(vorvx{t},bnd_pnts,p2,pos,length(p_int),coef,adv,type,p_int);
    if type == 2 || type == 3
        for i = 1:size(pos,1)
            if ~ismember(i,adv)
                pos(i,:) = sum3{t}(i,:);
            end
        end
    else
        pos = sum3{t};
    end
    p_sav{t+1} = pos;
end
out = p_sav{t+1};
