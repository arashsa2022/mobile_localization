function [xy_est_BTS, sigma_BTS, x_est_BTS, y_est_BTS, ind] = xy_estimator(LOS_TOA_sigma,NLOS_TOA_sigma,LOS_AOA_sigma,NLOS_AOA_sigma,n)
global BTS_xy
global BTS_dist
global BTS_AOA
global xy
[~,idx] = maxk(BTS_dist,n,2);
ind = sub2ind(size(BTS_dist), repmat((1:size(BTS_dist,1))',[1 n]), idx);
noisy_dist = BTS_dist + LOS_TOA_sigma*randn(size(BTS_dist));
noisy_dist(ind) = noisy_dist(ind) + NLOS_TOA_sigma*randn(size(ind));

noisy_angle = BTS_AOA + LOS_AOA_sigma*randn(size(BTS_AOA));
noisy_angle(ind) = noisy_angle(ind) + NLOS_AOA_sigma*randn(size(ind));

x_est_BTS = repmat(BTS_xy(:,1)',[size(BTS_dist,1) 1]) + noisy_dist.*cos(noisy_angle);
y_est_BTS = repmat(BTS_xy(:,2)',[size(BTS_dist,1) 1]) + noisy_dist.*sin(noisy_angle);
xy_est = zeros(size(BTS_dist,1),size(BTS_dist,2),2);
xy_est(:,:,1) = x_est_BTS;
xy_est(:,:,2) = y_est_BTS;

temp = reshape(permute(xy_est,[2 1 3]),[numel(BTS_dist) 2]);
temp_dist = pdist2(temp,temp);
temp_dist = temp_dist .* kron(eye(size(BTS_dist,1)),ones(size(BTS_dist,2)));

merit = sum(temp_dist,2);
merit = reshape(merit,[size(BTS_dist,2) size(BTS_dist,1)]);
[~,idx] = mink(merit,floor(size(BTS_dist,2)/2)+1,1);
ind = sub2ind(size(x_est_BTS), repmat((1:size(xy_est,1))',[1 floor(size(BTS_dist,2)/2)+1]), idx');
x_est_sel = x_est_BTS(ind);
x_est_sel = mean(x_est_sel,2);
y_est_sel = y_est_BTS(ind);
y_est_sel = mean(y_est_sel,2);
xy_est_BTS = [x_est_sel y_est_sel];
error = xy_est_BTS-xy;
sigma_BTS = std(error(:));
end