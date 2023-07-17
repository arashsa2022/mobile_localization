function [error_est, error_est_BTS] = RMSE(nodes_list,roads_list,map_matching_indicator)
global xy
global xy_est_BTS
xy_est = zeros(size(xy));
switch map_matching_indicator
    case 'ON'
        for i=1:length(nodes_list)
            node = nodes_list(i);
            [x_matched,y_matched] = node.map_matching(roads_list);
            xy_est(i,:) = [x_matched,y_matched];
        end
    case 'OFF'
        for i=1:length(nodes_list)
            node = nodes_list(i);
            x_matched = node.a_hat(1,2);
            y_matched = node.a_hat(1+size(node.a_hat,1)/2,2);
            xy_est(i,:) = [x_matched,y_matched];
        end   
    otherwise
        for i=1:length(nodes_list)
            node = nodes_list(i);
            x_matched = node.a_hat(1,2);
            y_matched = node.a_hat(1+size(node.a_hat,1)/2,2);
            xy_est(i,:) = [x_matched,y_matched];
        end
end
error_est = sum(diag(pdist2(xy,xy_est)))/sqrt(length(nodes_list));
error_est_BTS = sum(diag(pdist2(xy,xy_est_BTS)))/sqrt(length(nodes_list));
end