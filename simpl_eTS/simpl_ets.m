function [ y_daszek, R_w_czasie, opis ] = simpl_ets( x, y, r, OMEGA, opis )
    % x macierz [n x k] gdzie k liczba danych
    % y prawdziwe wyjscie jakie zaobserwowano [1 x k]
    % r is a positive constant, which defines the zone of influence of the i-th fuzzy rule antecedent
    opis{2} = 'simpl eTS';
    
    % psi(k) - [(n+1), R]
    % THETA(k) - [(n+1), mR]
    % y_daszek(k) - [m, 1]
    % C(k) - [m, m]
    n = length(x(:,1));
    m = length(y(:,1));
    K = length(x);
    R = 1;
    
    x_gwiazdka = cell(R, 1);
    x_gwiazdka{R} =  x(:, 1);
    
    % stage one - inicjalizacja parametrów
    z = cell(K, 1); % read z(1) 
    z{1} = [x(:,1)', y(:,1)']';
    
    S_gwiazdka = cell(K, R);
    S_gwiazdka{1, R} = 0;
    
    z_gwiazdka = cell(R, 1);
    z_gwiazdka{R} = z{1};
    
    % (15)
    PI = cell(K, 1);
    PI{1} = zeros(m, (n+1));
    
    THETA = cell(K, R);
    THETA{1, R} = PI{1}';
    psi = cell(K, R);
    psi{1, R} = zeros((n+1), 1);
        
    y_daszek = cell(K, 1);
    y_daszek{1} = zeros(1, m);
    
    C = cell(K, R);
    C{1, R} = OMEGA * eye(n+1);
    
    % w celach statystycznuch ilosc klastrow w czasie k
    R_w_czasie = zeros(1,K);
    R_w_czasie(1) = R;
    
    k = 1;
    while k <= K - 1
                
        % stage 2
        x_k = x(:, k + 1); % k-ta kolumna; k-ta dana
        
        x_k_e = [1 x_k']'; % [1, x_e]
        
        % stage 3
        % Estimate the output: y_daszek(k + 1) by (16)
        % (16) y_daszek(k + 1) = psi^T(k + 1)THETA(k)
        for i=1:R,
            mi1 = mi(x_k, x_gwiazdka, r, R);
            lambda1 = lambda( mi1, i );
            
           psi{k + 1, i} = lambda1 * x_k_e;
        end
        y_daszek{k + 1} = zeros(1, m);
        for i=1:R,
            y_daszek{k + 1} = y_daszek{k + 1} + psi{k + 1, i}' * THETA{k, i};
        end
        
        % stage 4
        k = k + 1;
        % Read y(k)
        y_k = y(:, k);
        
        % stage 5
        % Calculate S_k recursively by (19) - (20)
        % z(k) = [ x(k), y(k) ]'
        z{k} = [x_k', y_k']';

        S_k = fS_k(k, n, m, z);
     
        % stage 6
        % update S(z_i_gwiazdka) each centre by 23a
        for i=1:R,
            S_gwiazdka{k, i} = (k-2)/(k-1) * S_gwiazdka{k-1, i} + sumsqr(z{k}-z{k-1});
        end;
        
        % Compare S(k) with S_k(z_i_gw)
        % stage 7
        retu_24 = s_ets_24( S_k, [S_gwiazdka{k, :}]);
        [ retu_25, l ] = s_ets_25( x_k, x_gwiazdka, R, r);
        
        if retu_24 && retu_25
%            disp(['podmien klaster l=', num2str(l), '   k=', num2str(k)]);
           
           x_gwiazdka{l} = x_k;
           S_gwiazdka{k, l} = S_k;
        elseif retu_24
%             disp(['nowy    klaster R=', num2str(R+1), '   k=', num2str(k)]);
            % (23)
            
            R = R + 1;
            x_gwiazdka{R} =  x_k;
            S_gwiazdka{k, R} = S_k;
            
            z_gwiazdka{R} = z{k};
            
            % dla matlaba powikszam tablice    
            for i=1:k,
                C{i, R} = OMEGA * eye(n+1);
                THETA{i, R} = PI{1}';
            end
            psi{k, R} = zeros((n+1), 1);
        else
%             disp(['nie rob nic 24=',int2str(retu_24), ' 25=',int2str(retu_25), ' k=', num2str(k)]);
        end
        
        % Estimate the parameters of local sub-models by RLS (31)-(32)        
        
        for i=1:R,
            C{k, i} = C{k - 1, i} - ( C{k - 1, i} * psi{k} * psi{k}' * C{k - 1, i} ) / ( 1 + psi{k}' * C{k - 1, i} * psi{k});
            THETA{k, i} = THETA{k - 1, i} + C{k} * psi{k, i} * ( y_k - psi{k, i}' * THETA{k - 1, i} );
        end
        
       % k = K+1;
        R_w_czasie(k) = R;
    end
    opis{3} = R;
    opis{4} = r;
    opis{5} = OMEGA;
    opis{6} = z_gwiazdka;
    opis{7} = x;
end

function [ retu_25, l ] = s_ets_25( x_k, x_gwiazdka, R, r)
    min = sumsqr(x_k - x_gwiazdka{1});
    l = 1;
    for i=2:R,
        new_min = sumsqr(x_k - x_gwiazdka{i});
        if new_min < min
            min = new_min;
            l = i;
        end
    end 
    
    retu_25 = min < ( r / 2 );
end

function [ retu_24 ] = s_ets_24( S_k, S_k_gwiazdka)
    maxx = max(S_k_gwiazdka);
    minn = min(S_k_gwiazdka);
    
    retu_24 = S_k <= minn || S_k >= maxx;
end

function [ S_k ] = fS_k(k, n, m, z )
    % (19)-(20)
    a = 1 / ( (k - 1)*(n + m) );
    
    new_gamma = gamma_k(z, k);
    
    new_z_k = 0;
    for j=1:n+m,
        beta_j =  beta_j_k(z, k, j);
        new_z_k = new_z_k + z{k}(j) * beta_j;
    end
    
    b = (k-1) * sumsqr(z{k}) - 2 * new_z_k + new_gamma;

    S_k = a * b;
end

function [ beta_j_k ] = beta_j_k( z, k, j )
    sum = 0;
    for l=1:(k-1),
        sum = sum + z{l}(j);
    end
    beta_j_k = sum;
end

function [ gamma_k ] = gamma_k( z, k )
    sum = 0;
    for l=1:(k-1),
        sum = sum + sumsqr(z{l});
    end
    gamma_k = sum;
end

function [ lambda ] = lambda( mi, i )
    % lamda - skalar
    lambda = mi{i} / sum([mi{:}]);
end

function [ mi ] = mi(x, x_gwiazdka, r, R)
    % mi - skalar
    mi = cell(R, 1);
    for i=1:R,
        out = 1;
        for j=1:length(x),
            out = out * ( 1 + ( (2/r) * (x(j) - x_gwiazdka{i}(j)) )^2 );
        end
        mi{i} = 1 / out;
    end
end