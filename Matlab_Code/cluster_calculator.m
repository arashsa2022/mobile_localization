function [nodes_cluster, dist, BTS_dist, BTS_AOA] = cluster_calculator(Rc)
global xy;
global BTS_xy;
dist = pdist2(xy,xy);
nodes_cluster = (dist < Rc);
BTS_dist = pdist2(xy,BTS_xy);
c = repmat(xy,[1 1 size(BTS_xy,1)]);
f = repmat(BTS_xy,[1 1 size(xy,1)]);
f = permute(f,[3 2 1]);
g = c-f;
g = bsxfun(@plus,g(:,1,:),1i*g(:,2,:));
g = permute(g,[1 3 2]);
BTS_AOA = angle(g);
end