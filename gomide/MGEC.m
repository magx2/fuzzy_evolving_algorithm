function [ v, SIGMA, cluster_created, cluster_merged, c, idx, o ] = MGEC( k, x_k, v, SIGMA, lambda, omega, SIGMA_init, c, o )
% setting to false
cluster_created = 0;
cluster_merged = 0;
Ta = 1 - lambda;

if k == 1
    % innitialize first cluster
    v(:, 1) = x_k;
    SIGMA(1) = SIGMA_init;
    c = 1;
end
% compute p_i and a_i for all clusters
p = zeros(c, 1);
a = zeros(c, 1);
for i=1:c,
    
    g = gomide_9(x_k, v(:, i), SIGMA(i))
    p(:, i) = gomide_9(x_k, v(:, i), SIGMA(i));
    p_i = p(:, i);
    
    % calculate threshold
    Tp = gomide_threshold();
    if p_i < Tp
        o(i, k) = 1;
    else
        o(i, k) = 0;
    end
    
    if k > omega
        nv_i = sum(c(i,:));
        a(i) = 1; % TODO - will always add new cluster if p_i < Tp
    else
       a(i) = 0; 
    end
end
idx = max(p);
if p_i < Tp && a(idx) > Ta
    % create new cluster 
    c = c + 1;
    v(c) = x_k;
    SIGMA(c) = SIGMA_init;
    idx = c;
    cluster_created = c;
else
    % update existing cluster
    v_idx =  v(:, idx);
    x_v = (x_k - v_idx);
    alpha = 0.1; % alpha = learning rate; alpha = [0, 1]
    G_idx = gomide_G(alpha, p(idx), a(idx));
    v(:, idx) = v_idx + G_idx * x_v;
    SIGMA(idx) = (1 - G_idx) * (SIGMA(idx) - G_idx * x_v' * x_v);
end
% check for reduntand clusters
v_idx = v(idx);
for i=1:c,
    if gomide_15(v_idx, v(i), SIGMA(idx)) > Tp ||  gomide_15(v(i), v_idx, SIGMA(i))
        v(idx) = (v(i) + v_idx) ./ 2;
        SIGMA(idx) = SIMGA_init;
        c = c - 1;
        cluster_merged = [idx, i]; % is cluster merging only onces in this loop?
    end
end
end

function [ p_i ] = gomide_9( x, v_i,  SIGMA_i )
tmp = - gomide_7(x, v_i, SIGMA_i) / 2;
p_i = exp(tmp);
end

function [ M ] = gomide_7( x, v_i,  SIGMA_i )
M = (x - v_i)' * (1 / SIGMA_i) * (x - v_i);
end

function [ M ] = gomide_15( v_i, v_j,  SIGMA_i )
M = gomide_7(v_i, v_j, SIGMA_i);
end

function [ Tp ] = gomide_threshold()
Tp = 0.5;%TODO
end

function [ G ] = gomide_G( alpha, p_idx, a_idx )
G = alpha * p_idx ^ (1 - a_idx);
end