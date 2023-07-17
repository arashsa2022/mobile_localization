clear;
clc;
%%
lane_width = 3.66;
Rc = 50;    %communication range
nodes_list = [];
Ts = 1;
T = 90;
N = 120;    %number of nodes
N_BTS = 9;  %number of BTS
k = N/12;   %number of nodes per road
LOS_TOA_sigma = 20;
NLOS_TOA_sigma = 50;
% LOS_AOA_sigma = pi/100;
LOS_AOA_sigma = pi/80;
NLOS_AOA_sigma = pi/30;
% NLOS_AOA_sigma = pi/60;
n = 6;  %number of NLOS BTS
global xy
global vxy
global nodes_cluster
global dist
global BTS_dist
global BTS_AOA
global BTS_xy;
global xy_est_BTS
global sigma_BTS
global x_est_BTS
global y_est_BTS
global ind
xy = zeros(N,2);    %position of nodes
vxy = zeros(N,2);   %velocity of nodes
BTS_xy = zeros(N_BTS,2);    %position of BTS
nodes_cluster = [];
dist = [];  %distnace between nodes
BTS_dist = [];  %distance between nodes and BTS
BTS_AOA = [];   %AOA from each node to each BTS
xy_est_BTS = [];    %position estimated by BTS
sigma_BTS = []; %standard deviation of BTS estimation
x_est_BTS = []; %x estimation of all BTS
y_est_BTS = []; %y estimation of all BTS
ind = [];   %indices of selected BTS
dist_sigma = 3;
roads_list = [];
BTS_list = [];
R = 500*ones(3,2);
C = [400*ones(3,1) 600*ones(3,1)];
points = [0 0;R(1,1) 0;R(1,1)+R(1,2) 0;0 C(1,1);R(2,1) C(2,1);R(2,1)+R(2,2) C(3,1);0 C(1,1)+C(1,2);R(3,1) C(2,1)+C(2,2);R(3,1)+R(3,2) C(3,1)+C(3,2)];
%%  R11 Road

roads_list = [roads_list;road(points(2,:),points(1,:),R(1,1),lane_width)];

%%  R12 Road

roads_list = [roads_list;road(points(2,:),points(3,:),R(1,2),lane_width)];

%%  R21 Road

roads_list = [roads_list;road(points(5,:),points(4,:),R(2,1),lane_width)];

%%  R22 Road

roads_list = [roads_list;road(points(5,:),points(6,:),R(2,2),lane_width)];

%%  R31 Road

roads_list = [roads_list;road(points(8,:),points(7,:),R(3,1),lane_width)];

%%  R32 Road

roads_list = [roads_list;road(points(8,:),points(9,:),R(3,2),lane_width)];

%%  c11 Road

roads_list = [roads_list;road(points(4,:),points(1,:),C(1,1),lane_width)];

%%  c12 Road

roads_list = [roads_list;road(points(7,:),points(4,:),C(1,2),lane_width)];

%%  c21 Road

roads_list = [roads_list;road(points(5,:),points(2,:),C(2,1),lane_width)];

%%  c22 Road

roads_list = [roads_list;road(points(5,:),points(8,:),C(2,2),lane_width)];

%%  c31 Road

roads_list = [roads_list;road(points(3,:),points(6,:),C(3,1),lane_width)];

%%  c32 Road

roads_list = [roads_list;road(points(6,:),points(9,:),C(3,2),lane_width)];

%%
spacingR11 = exprnd(20,[k 1]);
locR11 = cumsum(spacingR11);
for i=1:k  %R11
    y0 = roads_list(1).starting_point(2);
    x0 = roads_list(1).starting_point(1) - locR11(i);
    vx = -50*1000/3600;
    vy = 0;
    vx_sigma = 1;
    vy_sigma = 0.2;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(1),vx_sigma,vy_sigma,dist_sigma)];
end

spacingR12 = exprnd(20,[k 1]);
locR12 = cumsum((spacingR12));
for i=(k+1):(2*k)   %R12
    y0 = roads_list(2).starting_point(2);
    x0 = roads_list(2).starting_point(1) + locR12(i-k);
    vx = 50*1000/3600;
    vy = 0;
    vx_sigma = 1;
    vy_sigma = 0.2;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(2),vx_sigma,vy_sigma,dist_sigma)];
end

spacingR21 = exprnd(35,[k 1]);
locR21 = cumsum((spacingR21));
for i=(2*k+1):(3*k)    %R21
    y0 = roads_list(3).starting_point(2);
    x0 = roads_list(3).starting_point(1) - locR21(i-2*k);
    vx = -80*1000/3600;
    vy = 0;
    vx_sigma = 1;
    vy_sigma = 0.2;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(3),vx_sigma,vy_sigma,dist_sigma)];
    
end


spacingR22 = exprnd(35,[k 1]);
locR22 = cumsum((spacingR22));
for i=(3*k+1):(4*k) %R22
    y0 = roads_list(4).starting_point(2);
    x0 = roads_list(4).starting_point(1) + locR22(i-3*k);
    vx = 80*1000/3600;
    vy = 0;
    vx_sigma = 1;
    vy_sigma = 0.2;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(4),vx_sigma,vy_sigma,dist_sigma)];
    
end


spacingR31 = exprnd(20,[k 1]);
locR31 = cumsum((spacingR31));
for i=(4*k+1):(5*k)  %R31
    y0 = roads_list(5).starting_point(2);
    x0 = roads_list(5).starting_point(1) - locR31(i-4*k);
    vx = -60*1000/3600;
    vy = 0;
    vx_sigma = 1;
    vy_sigma = 0.2;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(5),vx_sigma,vy_sigma,dist_sigma)];
    
end

spacingR32 = exprnd(20,[k 1]);
locR32 = cumsum((spacingR32));
for i=(5*k+1):(6*k)   %R32
    y0 = roads_list(6).starting_point(2);
    x0 = roads_list(6).starting_point(1) + locR32(i-5*k);
    vx = 60*1000/3600;
    vy = 0;
    vx_sigma = 1;
    vy_sigma = 0.2;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(6),vx_sigma,vy_sigma,dist_sigma)];
    
end

spacingC11 = exprnd(15,[k 1]);
locC11 = cumsum((spacingC11));
for i=(6*k+1):(7*k)  %C11
    y0 = roads_list(7).starting_point(2) - locC11(i-6*k);
    x0 = roads_list(7).starting_point(1);
    vx = 0;
    vy = -40*1000/3600;
    vx_sigma = 0.2;
    vy_sigma = 1;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(7),vx_sigma,vy_sigma,dist_sigma)];
    
end

spacingC12 = exprnd(15,[k 1]);
locC12 = cumsum((spacingC12));
for i=(7*k+1):(8*k)    %C12
    y0 = roads_list(8).starting_point(2) - locC12(i-7*k);
    x0 = roads_list(8).starting_point(1);
    vx = 0;
    vy = -40*1000/3600;
    vx_sigma = 0.2;
    vy_sigma = 1;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(8),vx_sigma,vy_sigma,dist_sigma)];
    
end

spacingC21 = exprnd(40,[k 1]);
locC21 = cumsum((spacingC21));
for i=(8*k+1):(9*k)    %C21
    y0 = roads_list(9).starting_point(2) - locC21(i-8*k);
    x0 = roads_list(9).starting_point(1);
    vx = 0;
    vy = -80*1000/3600;
    vx_sigma = 0.2;
    vy_sigma = 1;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(9),vx_sigma,vy_sigma,dist_sigma)];
    
end


spacingC22 = exprnd(40,[k 1]);
locC22 = cumsum((spacingC22));
for i=(9*k+1):(10*k) %C22
    y0 = roads_list(10).starting_point(2) + locC22(i-9*k);
    x0 = roads_list(10).starting_point(1);
    vx = 0;
    vy = 80*1000/3600;
    vx_sigma = 0.2;
    vy_sigma = 1;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(10),vx_sigma,vy_sigma,dist_sigma)];
    
end


spacingC31 = exprnd(30,[k 1]);
locC31 = cumsum((spacingC31));
for i=(10*k+1):(11*k)    %C31
    y0 = roads_list(11).starting_point(2) + locC31(i-10*k);
    x0 = roads_list(11).starting_point(1);
    vx = 0;
    vy = 60*1000/3600;
    vx_sigma = 0.2;
    vy_sigma = 1;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(11),vx_sigma,vy_sigma,dist_sigma)];
    
end

spacingC32 = exprnd(30,[k 1]);
locC32 = cumsum((spacingC32));
for i=(11*k+1):(12*k)   %C32
    y0 = roads_list(12).starting_point(2) + locC32(i-11*k);
    x0 = roads_list(12).starting_point(1);
    vx = 0;
    vy = 60*1000/3600;
    vx_sigma = 0.2;
    vy_sigma = 1;
    nodes_list = [nodes_list;node(i,x0,y0,vx,vy,Rc,roads_list(12),vx_sigma,vy_sigma,dist_sigma)];
end
%%  BTS
for i=1:N_BTS
    BTS_list = [BTS_list;BTS(i,points(i,:))];
end

%%  Simulation

%%  Initialize
Initializer(nodes_list,n,Rc);
%%  next step
mobility_simulator(1,Ts,LOS_TOA_sigma,NLOS_TOA_sigma,LOS_AOA_sigma,NLOS_AOA_sigma,n,Rc);
avr_simulator(nodes_list);
%%  simulation
error = zeros(T,3);

for i=1:T
    location_simulation(nodes_list,Ts)
    [error(i,1),error(i,2)] = RMSE(nodes_list,roads_list,'OFF');
    [error(i,3),~] = RMSE(nodes_list,roads_list,'ON');
    mobility_simulator(1,Ts,LOS_TOA_sigma,NLOS_TOA_sigma,LOS_AOA_sigma,NLOS_AOA_sigma,n,Rc);
    avr_simulator(nodes_list);
    arash(i)=sigma_BTS;
end
%%  Error Diagram vs Time
fig = figure;
plot(1:T,error(:,1));
hold on
plot(1:T,error(:,3));
plot(1:T,error(:,2));
xlabel('Time (s)');
ylabel('Error (m)');
legend('Kalman Error without map matching','Kalman Error with map matching','BTS Error');
title('Error Diagram vs Time, Number of Users = ' + string(N));
%%  Real Position Simulation
node_id = 10;
node = nodes_list(node_id);
real_xy = zeros(T,2);
est_xy_without_map_matching = zeros(T,2);
est_xy_with_map_matching = zeros(T,2);
est_BTS = zeros(T,2);
%%  Initialize
Initializer(nodes_list,n,Rc)
%%  next step
mobility_simulator(1,Ts,LOS_TOA_sigma,NLOS_TOA_sigma,LOS_AOA_sigma,NLOS_AOA_sigma,n,Rc);
avr_simulator(nodes_list);
%%  simulation
for i=1:T
    location_simulation(nodes_list,Ts)
    real_xy(i,:) = xy(node_id,:);
    est_xy_without_map_matching(i,:) = [node.a_hat(1,2) node.a_hat(1+size(node.a_hat,1)/2,2)];
    [x_matched,y_matched] = node.map_matching(roads_list);
    est_xy_with_map_matching(i,:) = [x_matched,y_matched];
    est_BTS(i,:) = xy_est_BTS(node_id,:);
    mobility_simulator(1,Ts,LOS_TOA_sigma,NLOS_TOA_sigma,LOS_AOA_sigma,NLOS_AOA_sigma,n,Rc);
    avr_simulator(nodes_list);
end
%%  Diagram
figure();
plot(est_xy_without_map_matching(:,1),est_xy_without_map_matching(:,2),'-x');
hold on
plot(est_xy_with_map_matching(:,1),est_xy_with_map_matching(:,2),'-+');
plot(real_xy(:,1),real_xy(:,2),'-*');
plot(est_BTS(:,1),est_BTS(:,2),'--');
xlabel('x');
ylabel('y');
legend('without map matching','with map matching','Real Position','BTS');

%%  Simulation for BTS Diagram

Desired_Time_Steps = randi(T,[4 1]);
Nodes = 1:k:6*k;
colormap = ['r';'g';'b';'y';'c';'m'];
%%  Initialize
Initializer(nodes_list,n,Rc)
%%  next step
mobility_simulator(1,Ts,LOS_TOA_sigma,NLOS_TOA_sigma,LOS_AOA_sigma,NLOS_AOA_sigma,n,Rc);
avr_simulator(nodes_list);
%%  simulation
for i=1:T
    location_simulation(nodes_list,Ts)
    [~,~] = RMSE(nodes_list,roads_list,'OFF');
    
    if (sum(i == Desired_Time_Steps))
        fig = figure();
        plt = zeros(5,1);
        plt(1) = scatter(BTS_xy(:,1),BTS_xy(:,2),100,'k','d','filled');
        hold on
        x_est_sel = x_est_BTS(ind);
        y_est_sel = y_est_BTS(ind);
        for node = Nodes
        plt(2) = scatter(x_est_BTS(node,:),y_est_BTS(node,:),50,colormap(floor(node/k)+1));
        plt(3) = scatter(x_est_sel(node,:),y_est_sel(node,:),30,colormap(floor(node/k)+1),'s','filled');
        plt(4) = scatter(xy(node,1),xy(node,2),50,colormap(floor(node/k)+1,:),'*');
        end
        for r = 1:length(roads_list)
            x = [roads_list(r).starting_point(1);roads_list(r).ending_point(1)];
            y = [roads_list(r).starting_point(2);roads_list(r).ending_point(2)];
            plt(5) = line(x,y);
        end
        title('Time = ' + string(Desired_Time_Steps(Desired_Time_Steps == i)) + ', Number of Users = ' + string(N));
        legend(plt,'BTS','BTS Estimation','Max Density Selection','Real Location','Roads');
        xlimit = xlim;
        ylimit = ylim;
        xlim([xlimit(1)-100,xlimit(2)+100]);
        ylim([ylimit(1)-100,ylimit(2)+100]);
    end
    
    mobility_simulator(1,Ts,LOS_TOA_sigma,NLOS_TOA_sigma,LOS_AOA_sigma,NLOS_AOA_sigma,n,Rc);
    avr_simulator(nodes_list);  
end


