classdef node < handle
    properties
        id
        x0
        y0
        x
        y
        vx
        vy
        Rc
        x_min
        x_max
        y_min
        y_max
        vx_sigma
        vy_sigma
        dist_sigma
        cluster_id
        a
        a_hat
        v
        v_past
        r
        P_hat
    end
    methods
        function obj = node(id,x0,y0,vx,vy,Rc,road,vx_sigma,vy_sigma,dist_sigma)
            global xy
            global vxy;
            obj.id = id;
            obj.x_min = min(road.starting_point(1),road.ending_point(1));
            obj.y_min = min(road.starting_point(2),road.ending_point(2));
            obj.x_max = max(road.starting_point(1),road.ending_point(1));
            obj.y_max = max(road.starting_point(2),road.ending_point(2));
            obj.x0 = x0;
            obj.y0 = y0;
            obj.x = x0;
            obj.y = y0;
            obj.vx = vx;
            obj.vy = vy;
            obj.Rc = Rc;
            obj.vx_sigma = vx_sigma;
            obj.vy_sigma = vy_sigma;
            obj.dist_sigma = dist_sigma;
            obj.cluster_id = [];
            obj.a = [];
            obj.v = [];
            obj.r = [];
            obj.P_hat = [];
            obj.a_hat = [];
            obj.v_past = [];
            xy(obj.id,:) = [obj.x obj.y];
            vxy(obj.id,:) = [obj.vx obj.vy];
        end
             
        function [a,v,r] = avr(obj)
            global vxy
            global nodes_cluster
            global dist
            global xy_est_BTS
            
            cluster_id_temp = (find(nodes_cluster(obj.id,:)))';
            cluster_id_temp(cluster_id_temp == obj.id) = [];
            cluster_id_temp = [obj.id;cluster_id_temp];
            
            a = xy_est_BTS(cluster_id_temp,:);
            a = a(:);
            a = [repmat(cluster_id_temp,[2 1]) a];
            obj.a = a;
            
            v = vxy(cluster_id_temp,:) + [obj.vx_sigma*randn(size(cluster_id_temp,1),1) obj.vy_sigma*randn(size(cluster_id_temp,1),1)];
            v = v(:);
            v = [repmat(cluster_id_temp,[2 1]) v];
            obj.v = v;
            
            cluster_dist = dist(cluster_id_temp,cluster_id_temp);
            mask = tril(ones(size(cluster_dist)),-1) .* (cluster_dist<(obj.Rc));
            r = cluster_dist(logical(mask));
            r = r + obj.dist_sigma*randn(size(r));
            [m,n] = find(mask);
            r = [cluster_id_temp(n) cluster_id_temp(m) r];
            obj.r = r;
            obj.cluster_id = cluster_id_temp;
        end
        
        function [a_k_hat_plus,P_k_plus]= estimate_location(obj,nodes_list,Ts)   %Remember to run avr before estimate_location
            global sigma_BTS
            global vxy
            
            a_k_past_hat_plus = obj.a_hat;
            P_k_past_plus = obj.P_hat;
            r_k = obj.r;
            a_k = obj.a;
            v_k_past = obj.v_past;
            
            [~,idx] = ismember(a_k_past_hat_plus(:,1),v_k_past(:,1));
            idx = [idx(1:size(idx,1)/2);idx(size(idx,1)/2+1:end)+size(idx,1)/2];
            v_k_past = v_k_past(idx,:);
            
            n = size(a_k,1)/2;
            n_past = size(a_k_past_hat_plus,1)/2;
            
            R_k = zeros(2*n+size(r_k,1));
            R_k(1:2*n,1:2*n) = (sigma_BTS^2)*eye(2*n);
            R_k(2*n+1:end,2*n+1:end) = (obj.dist_sigma^2)*eye(size(r_k,1));
            
            Gamma_k_past = [(obj.vx_sigma^2)*eye(n) zeros(n);zeros(n) (obj.vy_sigma^2)*eye(n)];
            
            a_k_hat_minus = zeros(2*n,2);
            a_k_hat_minus(:,1) = repmat(obj.cluster_id,[2 1]);
            cluster_id_past = a_k_past_hat_plus(1:n_past,1);
            
            ind = ismember(obj.cluster_id,cluster_id_past);
            add = find(ind == 0);
            
            if (~isempty(add))
                for i=1:size(add,1)
                    missing_node = nodes_list(obj.cluster_id(add(i)));
                    x_missing = missing_node.a_hat(1,2);
                    y_missing = missing_node.a_hat(1+(size(missing_node.a_hat,1)/2),2);
                    vx_missing = missing_node.v(1,2);
                    vy_missing = missing_node.v(1+(size(missing_node.v,1)/2),2);
                    P_missing = [missing_node.P_hat(1,1),missing_node.P_hat(1,1+(size(missing_node.P_hat,2)/2));missing_node.P_hat(1+(size(missing_node.P_hat,1)/2),1),missing_node.P_hat(1+(size(missing_node.P_hat,1)/2),1+(size(missing_node.P_hat,2)/2))];
                    
                    a_k_past_hat_plus = [a_k_past_hat_plus(1:size(a_k_past_hat_plus,1)/2,:);[missing_node.id x_missing];a_k_past_hat_plus(size(a_k_past_hat_plus,1)/2+1:end,:);[missing_node.id y_missing]];
                    v_k_past = [v_k_past(1:size(v_k_past,1)/2,:);[missing_node.id vx_missing];v_k_past(size(v_k_past,1)/2+1:end,:);[missing_node.id vy_missing]];
                    P_k_past_plus = [P_k_past_plus(1:size(P_k_past_plus,1)/2,:);zeros(1,size(P_k_past_plus,2));P_k_past_plus(size(P_k_past_plus,1)/2+1:end,:);zeros(1,size(P_k_past_plus,2))];
                    P_k_past_plus = [P_k_past_plus(:,1:size(P_k_past_plus,2)/2),zeros(size(P_k_past_plus,1),1),P_k_past_plus(:,size(P_k_past_plus,2)/2+1:end),zeros(size(P_k_past_plus,1),1)];
                    P_k_past_plus(size(P_k_past_plus,1)/2,size(P_k_past_plus,2)/2) = P_missing(1,1);
                    P_k_past_plus(size(P_k_past_plus,1)/2,size(P_k_past_plus,2)) = P_missing(1,2);
                    P_k_past_plus(size(P_k_past_plus,1),size(P_k_past_plus,2)/2) = P_missing(2,1);
                    P_k_past_plus(size(P_k_past_plus,1),size(P_k_past_plus,2)) = P_missing(2,2);
                end
            end
            
            ind = ismember(cluster_id_past,obj.cluster_id);
            omit = find(ind == 0);
            if (~isempty(omit))
                a_k_past_hat_plus(omit,:) = [];
                a_k_past_hat_plus(n_past+size(add,1)+omit-size(omit,1),:) = [];
                v_k_past(omit,:) = [];
                v_k_past(n_past+size(add,1)+omit-size(omit,1),:) = [];
                P_k_past_plus(omit,:) = [];
                P_k_past_plus(:,omit) = [];
                P_k_past_plus(n_past+size(add,1)+omit-size(omit,1),:) = [];
                P_k_past_plus(:,n_past+size(add,1)+omit-size(omit,1)) = [];
            end
            [~,idx] = ismember(a_k_hat_minus(:,1),a_k_past_hat_plus(:,1));
            idx = [idx(1:size(idx,1)/2);idx(size(idx,1)/2+1:end)+size(idx,1)/2];
            a_k_past_hat_plus = a_k_past_hat_plus(idx,:);
            P_k_past_plus = P_k_past_plus(idx,:);
            P_k_past_plus = P_k_past_plus(:,idx);
            v_k_past = v_k_past(idx,:);
            
            
            
            if (prod(a_k_hat_minus(:,1) == a_k_past_hat_plus(:,1)) && prod(a_k_hat_minus(:,1) == v_k_past(:,1)))
                a_k_hat_minus(:,2) = a_k_past_hat_plus(:,2) + Ts*v_k_past(:,2);    %Eq.10
%                 P_k_minus = (Ts^2)*Gamma_k_past + P_k_past_plus;    %Eq.11
                   P_k_minus = (Ts^2)*Gamma_k_past + P_k_past_plus;
                   
%                    A1=eye(size(a_k_past_hat_plus,1));
%                    A=cat(2,A1,A1);
%                    
%             mmm=cat(1,a_k_past_hat_plus(:,2),v_k_past(:,2));
%             mm=cat(1,a_k_past_hat_plus(:,2),v_k_past(:,2));
%             pp=cat(1,P_k_past_plus,P_k_past_plus);
%             kk=cat(1,Gamma_k_past,Gamma_k_past);
%             ppp=repmat(P_k_past_plus,2,2);
%             sisi=A*ppp*A';
                   
                   
                   
%                    [mm,PP] = ekf_predict1(cat(1,a_k_past_hat_plus(:,2),v_k_past(:,2)),repmat(P_k_past_plus,2,2),A,Gamma_k_past);
           
%                      mmm=cat(1,a_k_past_hat_plus(:,2),v_k_past(:,2));
%             mm=cat(1,a_k_past_hat_plus(:,2),v_k_past(:,2));
%             pp=cat(1,P_k_past_plus,P_k_past_plus);
%             kk=cat(1,Gamma_k_past,Gamma_k_past);
%             ppp=repmat(P_k_past_plus,2,2);
%             sisi=A*ppp*A';
            else
                error('Mismatch!');
            end
            
            
            H_k1 = eye(2*n);
            H_k2 = zeros(size(r_k,1),2*n);
            for i=1:size(r_k,1)
                loc = [find(obj.cluster_id == r_k(i,1)) find(obj.cluster_id == r_k(i,2))];
                H_k2(i,loc(1)) = (a_k_hat_minus(loc(1),2)-a_k_hat_minus(loc(2),2))/sqrt((a_k_hat_minus(loc(1),2)-a_k_hat_minus(loc(2),2))^2+(a_k_hat_minus(n+loc(1),2)-a_k_hat_minus(n+loc(2),2))^2);
                H_k2(i,loc(2)) = (-1)*H_k2(i,loc(1));
                H_k2(i,n+loc(1)) = (a_k_hat_minus(n+loc(1),2)-a_k_hat_minus(n+loc(2),2))/sqrt((a_k_hat_minus(loc(1),2)-a_k_hat_minus(loc(2),2))^2+(a_k_hat_minus(n+loc(1),2)-a_k_hat_minus(n+loc(2),2))^2);
                H_k2(i,n+loc(2)) = (-1)*H_k2(i,n+loc(1));
            end
            Hk = [H_k1;H_k2];
            
            K_k = P_k_minus*Hk'/(Hk*P_k_minus*Hk'+R_k); %Eq.14
            
            
            if (~isempty(r_k))
                z_k = [a_k(:,2);r_k(:,end)];  %Eq.9
            else
                z_k = a_k(:,2); %Eq.9
            end
            
%             [m,P] = ekf_update1(mm,PP,z_k,Hk,Gamma_k_past,h);
            
            a_k_hat_plus = a_k_hat_minus;
            
            a_k_hat_plus(:,2) = a_k_hat_minus(:,2) + K_k*(z_k-Hk*a_k_hat_minus(:,2));  %Eq.12
            P_k_plus = P_k_minus - K_k*Hk*P_k_minus;    %Eq.13
            
            obj.a_hat = a_k_hat_plus;
            obj.P_hat = P_k_plus;
            
            v_temp = vxy(obj.cluster_id,:) + [obj.vx_sigma*randn(size(obj.cluster_id,1),1) obj.vy_sigma*randn(size(obj.cluster_id,1),1)];
            v_temp = v_temp(:);
            v_temp = [repmat(obj.cluster_id,[2 1]) v_temp];
            obj.v_past = v_temp;
            
        end
        
        function [x_matched,y_matched] = map_matching(obj,roads_list)
            x_est = obj.a_hat(1,2);
            y_est = obj.a_hat(1+size(obj.a_hat,1)/2,2);
            vx_meas = obj.v(1,2);
            vy_meas = obj.v(1+size(obj.v,1)/2,2);
            
            metric = zeros(length(roads_list),1);
            for i=1:length(roads_list)
                x1 = roads_list(i).starting_point(1);
                y1 = roads_list(i).starting_point(2);
                x2 = roads_list(i).ending_point(1);
                y2 = roads_list(i).ending_point(2);
                al = y2-y1;
                bl = x1-x2;
                cl = (x2-x1)*y1 - (y2-y1)*x1;
                d = abs(al*x_est + bl*y_est + cl)/sqrt(al^2+bl^2);
                cos = dot([vx_meas,vy_meas],[x2-x1,y2-y1])/(sqrt(vx_meas^2+vy_meas^2) * sqrt((y2-y1)^2 + (x2-x1)^2));
                metric(i) = cos/d;
            end
            [~,idx] = max(metric);
            x1 = roads_list(idx).starting_point(1);
            y1 = roads_list(idx).starting_point(2);
            x2 = roads_list(idx).ending_point(1);
            y2 = roads_list(idx).ending_point(2);
            al = y2-y1;
            bl = x1-x2;
            cl = (x2-x1)*y1 - (y2-y1)*x1;
            x_matched = (bl*(bl*x_est - al*y_est) - al*cl)/(al^2+bl^2);
            y_matched = (al*(-bl*x_est + al*y_est) - bl*cl)/(al^2+bl^2);
            
        end
        
        function reset(obj)
            global xy
            global vxy;
            obj.cluster_id = [];
            obj.a = [];
            obj.v = [];
            obj.r = [];
            obj.x = obj.x0;
            obj.y = obj.y0;
            obj.P_hat = [];
            obj.a_hat = [];
            obj.v_past = [];
            xy(obj.id,:) = [obj.x obj.y];
            vxy(obj.id,:) = [obj.vx obj.vy];
        end
    end
end
