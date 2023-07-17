function avr_simulator(nodes_list)
for i=1:length(nodes_list)
    node = nodes_list(i);
    [~,~,~] = node.avr();
end

end