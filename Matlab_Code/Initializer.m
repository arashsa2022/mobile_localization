function Initializer(nodes_list, n, Rc)
global nodes_cluster
global dist
global BTS_AOA
global BTS_dist
global xy_est_BTS
global sigma_BTS
global x_est_BTS
global y_est_BTS
global ind
for i=1:length(nodes_list)
    node = nodes_list(i);
    node.reset();
end
[nodes_cluster, dist, BTS_dist, BTS_AOA] = cluster_calculator(Rc);
[xy_est_BTS, sigma_BTS, x_est_BTS, y_est_BTS, ind] = xy_estimator(0,0,0,0,n);
avr_simulator(nodes_list);
for i=1:length(nodes_list)
    node = nodes_list(i);
    node.a_hat = node.a;
    node.P_hat = zeros(size(node.a,1));
    node.v_past = node.v;
end
end