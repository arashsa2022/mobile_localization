function mobility_simulator(k,Ts,LOS_TOA_sigma,NLOS_TOA_sigma,LOS_AOA_sigma,NLOS_AOA_sigma,n,Rc)
global xy;
global vxy;
global nodes_cluster
global dist
global BTS_dist
global BTS_AOA
global xy_est_BTS
global sigma_BTS
global x_est_BTS
global y_est_BTS
global ind
% xy = xy + k*Ts*vxy + gauss_rnd(zeros(120,2), 0.04);
xy = xy + k*Ts*vxy ;
[nodes_cluster, dist, BTS_dist, BTS_AOA] = cluster_calculator(Rc);
[xy_est_BTS, sigma_BTS, x_est_BTS, y_est_BTS, ind] = xy_estimator(LOS_TOA_sigma,NLOS_TOA_sigma,LOS_AOA_sigma,NLOS_AOA_sigma,n);
end