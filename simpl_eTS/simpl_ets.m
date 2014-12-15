function [ y_s ] = simpl_ets( x, y, r, OMEGA )
    % x macierz [n x k] gdzie k liczba danych
    % y prawdziwe wyjscie jakie zaobserwowano [1 x k]
    % r is a positive constant, which defines the zone of influence of the i-th fuzzy rule antecedent

    % psi(k) - [(n+1), 1]
    % THETA(k) - [(n+1), m]
    % y_daszek(k) - [m, 1]
    % C(k) - [m, m]
    n = length(x(:,1));
    m = length(y(:,1));
    K = length(x);
    R = 1;
    n_m_K = [n, m, K]
    
    x_gwiazdka = zeros(n, 1);
    x_gwiazdka(:,1) =  x(:, 1);
    
    % stage one - inicjalizacja parametrów
    z = zeros(n+m, K); % read z(1) 
    z(:, 1) = [x(:,1)', y(:,1)']'
    
    old_S_k_z_gwiazdka = zeros(R, 1);
    z_gwiazdka(1) = z(1);
    
    old_gamma = 1;
    old_beta = 1;
    
    % (15)
    PI = cell(K, 1);
    PI{1} = zeros(m, (n+1));
    
    THETA = cell(K, 1);
    THETA{1} = PI{1}';
    psi = cell(K);
    psi{1} = zeros(m, 1);
        
    y_daszek = cell(K, 1);
    y_daszek{1} = x(:, 1);
    
    C = cell(K);
    C{1} = OMEGA * eye(n+1);
    
    k = 1;
    while k <= K - 1
        
        % stage 2
        k
        x_k = x(:, k + 1); % k-ta kolumna; k-ta dana
        
        x_k_e = [1 x_k']'; % [1, x_e]
        
        % stage 3
        % Estimate the output: y_daszek(k + 1) by (16)
        % (16) y_daszek(k + 1) = psi^T(k + 1)THETA(k)
        psi_tmp = zeros(R, length(x_k_e));
        for i=1:R,
            mi = mi(x_k, x_gwiazdka(:,i), r);
            lambda = lambda( mi, i );
            
            psi_tmp(i, :) = lambda * x_k_e';
        end
        psi{k + 1} = psi_tmp';
        
        p = psi{k + 1}
        t = THETA{k}
        
        y_daszek{k + 1} = psi{k + 1}' * THETA{k};
        huj = y_daszek{k + 1}
        
        % stage 4
        k = k + 1;
        % Read y(k)
        y_k = y(:, k);
        
        % stage 5
        % Calculate S_k recursively by (19) - (20)
        % z(k) = [ x(k), y(k) ]'
        z(:, k) = [x_k', y_k']';

        S_k = fS_k(k, n, m, z);
     
        % stage 6
        % update S(z_i_gwiazdka) each centre by 23a
        for i=1:R,
            S_k_z_gwiazdka = (k-2)/(k-1) * old_S_k_z_gwiazdka(:, i) + ( z(:,k) - z(:, k-1) ).^2;
        end;
        old_S_k_z_gwiazdka = S_k_z_gwiazdka;
        
        % Compare S(k) with S_k(z_i_gw)
        % stage 7
        retu_24 = s_ets_24( S_k, S_k_z_gwiazdka);
        [ retu_25, l ] = s_ets_25( x_k, x_gwiazdka, R);
        
        if retu_24 && retu_25
           disp(['podmien klaster l=', num2str(l), ' k=', num2str(k)]);
           
        elseif retu_24
            disp(['nowy klaster R=', num2str(R+1), ' k=', num2str(k)]);
        end
        
        % Estimate the parameters of local sub=models by RLS (31)-(32)
        pis_k = psi{k}
        C{k} = C{k - 1} - ( C{k - 1} * psi{k} * psi{k}' * C{k - 1} ) * ( 1 + psi{k}' * C{k - 1} * psi{k});
        THETA{k} = THETA{k - 1} + C{k} * psi{k} * ( y_k - psi{k}' * THETA{k - 1 } );
    end 
end

function [ retu_25, l ] = s_ets_25( x_k, x_gwiazdka, R)
    min = sumsqr(x_k - x_gwiazdka(:,1));
    l = 1;
    for i=2:R,
        new_min = sumsqr(x_k - x_gwiazdka(:,i));
        if new_min < min
            min = new_min;
            l = i;
        end
    end 
    
    retu_25 = min < ( R / 2 );
end

function [ retu_24 ] = s_ets_24( S_k, S_k_z_gwiazdka)
sk = S_k
s = S_k_z_gwiazdka
    maxx = max(S_k_z_gwiazdka)
    minn = min(S_k_z_gwiazdka)
    
    retu_24 = S_k < minn || S_k > maxx;
end

function [ S_k ] = fS_k(k, n, m, z )
    % (19)-(20)
    a = 1 / ( (k - 1)*(n + m) );
    
    new_gamma = gamma_k(z, k);
    
    new_z_k = 0;
    for j=1:n+m,
        beta_j =  beta_j_k(z, k, j);
        new_z_k = new_z_k + z(j, k) * beta_j;
    end
    new_z_k = 2 * new_z_k;
    
    b = (k-1) * sumsqr(z(:,k)) - 2 * sum(new_z_k) * new_gamma

    S_k = a * b;
end

function [ beta_j_k ] = beta_j_k( z, k, j )
    sum = 0;
    for l=1:(k-1),
        sum = sum + z(j, l);
    end
    beta_j_k = sum;
end

function [ gamma_k ] = gamma_k( z, k )
    sum = 0;
    for l=1:(k-1),
        sum = sum + sumsqr(z(:, l));
    end
    gamma_k = sum;
end

function [ lambda ] = lambda( mi, i )
    % lamda - skalar
    lambda = mi(i) / sum(mi);
end

function [ mi ] = mi(x, x_gwiazdka, r)
    % mi - skalar
    out = 1;
    for j=1:length(x),
        out = out * ( 1 + ( (2/r) * (x(j) - x_gwiazdka(j)) )^2 );
    end

    mi = 1 / out;
end