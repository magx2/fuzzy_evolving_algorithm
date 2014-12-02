function [ y_s ] = eMG( x, y_d, lambda, omega, SIGMA_init )
% 
% y_d - pierwszy w
%

liczba_kolumn_x = length(x(1, :)); % liczba kolumn
liczba_wierszy_x = length(x);
c = 0; % c to  liczba klastrow

v = zeros(c + 1, liczba_kolumn_x);
v(1, :) = x(1, :);

SIGMA = cell(c + 1);
SIGMA{1} = SIGMA_init;

gamma = zeros(length(y_d), liczba_kolumn_x); % tutaj chyba trzeba dodac y_1 na sam poczatek a nie tylko x
gamma(:,1) = y_d;
% wymiery gamma: [liczba klastrow x liczba kolumn x]

Q_init = omega * eye(liczba_kolumn_x);
Q = cell(1,1);
Q{1,1} = Q_init;

o = zeros(liczba_kolumn_x, liczba_kolumn_x);

p = zeros(1,1);
y = zeros(1, liczba_kolumn_x); % wkladamy tutaj wiersze o dlugosci wiersza x
for k=1:liczba_wierszy_x, % lenght(x) - liczba wierszy
    x_k = x(k, :) % k-ty wiersz
    
    % compute the output
    for i=1:c,
        p(i) = gomide_9(x_k, v(c, :), SIGMA{i});
        y(:, i) = x_k * gamma(:, i);
    end
    y_s = weights_multiple_y(p, y, c) / weights(p); % wynik
    
%     [ v, SIGMA, cluster_created, cluster_merged, c, idx, o ] = MGEC(k, x(:, k), v, SIGMA, lambda, omega, SIGMA_init, c, o);
cluster_created=1;c=c+1
    if (cluster_created) && (c ~= 1) % dla pierwszego klastra mamy juz policzona gamme
       % jezeli zostal stworzony klaster to juz mamy c zwiekszone o 1!
       % create new rule
       gamma(:, c) = weights_multiple_y(p, gamma, c-1) / weights(p)
       Q{1, c} = Q_init;
    end
    
    % to jest robione w MGEC
    SIGMA{c} = SIGMA_init; 
    v(c, :) = x_k
    
    
    
%     for i=1:c, 
%        update consequent parameters
%        gom =  gomide_21_PSI(x(:, k), v(:, k), SIGMA, i) * ( y(i, k) - (x(:, k)' * gamma(i)) )
%        q0 = x(:, k)
%        q1 =  Q{i} * q0
%        q2 =  q1' * gom'
%        
%         gg = gamma(i) + q2
%         
%        gamma(i) = gamma(i) + (Q{i} * x(:, k))' * gomide_21_PSI(x(:, k), v(:, k), SIGMA, i) * ( y(i, k) - (x(:, k)' * gamma(i))' );
%        Q(i) = 
%     end
%     if cluster_merged
%        merge two clusetrs
%        
%     end
end
end

function [ PSI ] = gomide_21_PSI(x_k, v_k, SIGMA, i_param)
PSI = exp( (x_k- v_k(i_param))' * (1 / SIGMA(i_param)) * (x_k- v_k(i_param)) );

sum = 0;
for i=1:length(v_k(i_param)),
    sum = sum + gomide_21_helper(x_k, v_k(i), SIGMA(i));
end

PSI = PSI / sum;
end

function [ out ] = gomide_21_helper(x_k, v_i_k, SIGMA_i_k)
out = exp( (x_k - v_i_k)' * inv(SIGMA_i_k) * (x_k - v_i_k) );
end

function [ weights_multiple_y ] = weights_multiple_y( p, gamma, num )
    weights_multiple_y = zeros(length(gamma(:, 1)), 1);
    for i=1:num
        weights_multiple_y = weights_multiple_y + (p(i) * gamma(:, i));
    end
end

function [ weights ] = weights( p )
weights = sum(p);
end

function [ p_i ] = gomide_9( x, v_i,  SIGMA_i )
    tmp = - gomide_7(x, v_i, SIGMA_i) / 2;
    p_i = exp(tmp);
end

function [ M ] = gomide_7( x, v_i,  SIGMA_i )
    a = (x - v_i);
    b = inv(SIGMA_i);
    c = (x - v_i)';
        
    M = a * b * c;
end
