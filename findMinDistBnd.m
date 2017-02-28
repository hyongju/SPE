function out = findMinDistBnd(q1,pos,k,rEff)
% find_dist_min(q1,pos)
tmp_pl = (repmat(q1,size(pos,1),1) -pos).^2;
[val, idx] = sort(sqrt(tmp_pl(:,1) + tmp_pl(:,2)),1);


if isempty(idx)
    out = idx;
else
    j = 0;
    for i = 1:length(val)
        if val(i) <= rEff
            j = j+1;
            outTmp(j,1) = idx(i);
        end
    end
    if exist('outTmp')
        if size(outTmp,1) >=k
            out = outTmp(1:k,1);
        else
            out = outTmp;
        end
    else
        out = [];
    end
end
% out = idx(1:k,1);


