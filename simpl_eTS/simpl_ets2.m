function [ y_daszek, R_w_czasie, opis, S, S_min, S_max, S_podmiana, S_nowy ] = simpl_ets2( x, y, r, OMEGA, opis )
                                                                            debug=1;
                                                                            opis{2} = 'simpl eTS2';
    % A                                                                       
    n = length(x(:,1));
    K = length(x);
    k = 1;
    R = 1;
    
    x_g = cell(R, 1);
    x_g{R} = x(:, 1);
    
    z = cell(K, 10);
    z{k} = [x(:, k)' y(:, k)']';
    
    THETA = cell(R, 1);
    THETA{R} = zeros(n+1, 1);
    
    C = cell(R, 1);
    C{R} = OMEGA * eye(n+1);
    
    P = cell(R, 1);
    P{R} = 1;
    
	psi = cell(R,1);
    psi{R} = [ 1 x(:, 1)' ]';
    
    y_daszek = cell(K, 1);
    y_daszek{k} = 0;
    y_daszek{k+1} = 0;
    
                                                                            R_w_czasie = zeros(1,K);
                                                                            R_w_czasie(k) = R;
                                                                            
                                                                            S = cell(1, K);
                                                                            S{1}=P{R};

                                                                            S_min = cell(K,1);
                                                                            S_min{1}=P{R};

                                                                            S_max = cell(K,1);
                                                                            S_max{1}=P{R};

                                                                            S_podmiana= cell(1,K);
                                                                            S_podmiana{1}=0;

                                                                            S_nowy = cell(K,1);
                                                                            S_nowy{1}=1;
    for k=2:K
                                                                            R_w_czasie(k) = R;
       % B Reading the next data point
       x_k = x(:, k);
       y_k = y(:, k);
       x_k_e = [1 x_k']';
       z{k} = [x_k' y_k']';
       
       
        % PSI
        mi = cell(R,1);
        for i=1:R,
            mi{i} = exp( -(4) * sqrt( sumsqr( ((x_k - x_g{i})/r ) ) ) ); % (2)
        end
        suma_mi = sum( cell2mat(mi) );
        lambda = cell(R,1);
        for i=1:R,
           lambda{i} = mi{i} / suma_mi; % 5 
        end
        
        psi_k_minus_one = psi;
        psi = cell(R,1);
        for i=1:R,
            psi{i} = lambda{i} * x_k_e;
        end
        
        % C calculation of the potential of new data point
        % (8)
        suma=0;
        for i=1:(k-1),
            suma=suma + sqrt( sumsqr(z{i}-z{k}));
        end
        P_k = 1 / ( 1 + suma/(k-1) );
        
        % D updating potencials of centers
        % (11)
        for i=1:R,
            P{i} = (k-1) * P{i} / ( k - 2 + P{i} + P{i} * sqrt( sumsqr(z{k} - z{k-1}) ) );
        end
        
        % E Rule base evolution
        is_new_potencial_higher_than_others = P_k >= max(cell2mat(P));
        
        min_idx = 1;
        min_val = sqrt(sumsqr(z{k}-z{min_idx}));
        for i=1:R,
            val = sqrt(sumsqr(z{k}-z{i}));
            if val < min_val
                min_val = val;
                min_idx = i;
            end
        end
        is_new_z_close_to_old_one = min_val < (r/2);
        
        if is_new_potencial_higher_than_others && is_new_z_close_to_old_one
           % MODIFY - replace (min_idx) cluster
           if debug
                disp([ 'podmien klaster min_idx=', num2str(min_idx), ' k=', num2str(k) ]);
           end
           % (12)
           x_g{min_idx} = x_k;
           P{min_idx} = P_k;
        elseif is_new_potencial_higher_than_others
            if debug
              disp(['nowy    klaster R=', num2str(R+1), '   k=', num2str(k) ]);
            end
            % UPGRADE - add cluster
            R = R+1;
            x_g{R} = x_k;
            P{R} = P_k;            
            
            psi{R} = zeros(n+1,1);
            C{R} = OMEGA * eye(n+1);
            THETA{R} = zeros(n+1,1);
            for i=1:R-1,
                THETA{R} = lambda{i} * THETA{i};
            end
            
            for i=1:R-1,
                % C (21)
                C{i} = C{i} - ( lambda{i} * C{i} * psi_k_minus_one{i} * psi_k_minus_one{i}' * C{i} )/( 1 + psi_k_minus_one{i}' * C{i} * psi_k_minus_one{i} ); 
                % THETA (20)
                THETA{i} = THETA{i} + C{i} * psi_k_minus_one{i} * (y_k - y_daszek{k-1});
            end
        else
            
            if debug
              disp(['nie rob nic 24=',int2str(is_new_potencial_higher_than_others), ' 25=',int2str(is_new_z_close_to_old_one), ' k=', num2str(k) ]);
            end
        end
        
        % G output predict
        y_daszek{k+1} = 0;
        for i=1:R,
            y_daszek{k+1} = y_daszek{k+1} + psi{i}' * THETA{i};
        end
        
                                                                            s= cell2mat(P);
                                                                            S_min{k} = min(s);
                                                                            S_max{k} = max(s);
                                                                            S{k} = P_k;
    end
                                                                            opis{3} = R;
                                                                            opis{4} = r;
                                                                            opis{5} = OMEGA;
                                                                            %opis{6} = z_gwiazdka;
                                                                            opis{7} = x;
end