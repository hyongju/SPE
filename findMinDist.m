function out = findMinDist(q1,pos,k)
% find_dist_min(q1,pos)
tmp_pl = (repmat(q1,size(pos,1),1) -pos).^2;
[~, idx] = sort(sqrt(tmp_pl(:,1) + tmp_pl(:,2)),1);
out = idx(1:k,1);
