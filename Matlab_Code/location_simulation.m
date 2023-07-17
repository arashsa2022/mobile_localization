function location_simulation(nodes_list,Ts)
global sigma_BTS
for i=1:length(nodes_list)
    node = nodes_list(i);
    [~,~]= node.estimate_location(nodes_list,Ts);
end
end