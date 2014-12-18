function run_simpl_ets
format long
    % y = 2x_1 + 3x_2 + 4
%     x = [
%         1 2 3 4  5  6  7   8   9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30;
%         1 2 4 8 16 32 64 128 256 512 15 15 65 44 52 65 26 45 56 7 4 14  48 35 20 20 69 554 24 54
%     ];

    K=10;
    for k=1:K,
        x(:, k) = [ k k*2 ]';
    end
    
    y = zeros(1, length(x));
    for k=1:length(x),
        x_k = x(:, k);
      % y(k) =  x_k(1) * 2 + x_k(2) * 3 + 4;
      y(k) = sin(x_k(1) +x_k(2));
    end
    %plot3(x(1,:)', x(2,:)', y')

    r = 11;
    
    OMEGA = 1;
    
    simpl_ets( x, y, r, OMEGA )
end