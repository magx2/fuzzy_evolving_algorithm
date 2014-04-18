function [ y_s ] = eMG( x, y_d, lambda, omega, SIGMA_init )
x_length = length(x(1, :));
m = x_length;

gamma = zeros(1, x_length);

Q_init = omega * eye(m + 1);
Q = cell(1,1);
Q{1,1} = Q_init;

c = 0;

for k=1:length(x),
    p = zeros(1,c);
    y = zeros(1, c);
    x_k = x(k);
    
    % compute the output
    for i=1:c,
        p(i) = gomide_9(x_k, v_i, SIGMA(i));
        y(i) = x_k * gamma(i);
    end
end
end