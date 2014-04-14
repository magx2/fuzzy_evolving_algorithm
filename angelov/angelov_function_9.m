function [ y ] = angelov_function_9( x, tau, R, pi, T )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y = 0;
for i=1:R,
    y = y + lambda(tau, i) * xe(x)' * pi;
end
end



% MI
function [ mi ] = mi( r, x, x_gwiazdka )
    mi = exp(1)^(-4/r^2 * norm(x - x_gwiazdka)^2);
end
% END MI


% EPSILON
%function [ epsilon ] = epsilon( x, tau, R )
%    xe2 = xe(x)'
%    epsilon = [];
%    for i=1:R,
%        epsilon = [epsilon, lambda(tau, i) * xe2]
%    end
%end
%
function [ lambda ] = lambda( tau, i )
    suma = sum(tau);
    lambda = tau(i) / suma;
end

function [ xe ] = xe( x )
    xe = [1 ; x];
end
%% END EPSILON
%
%% THETA
%function [ theta ] = theta( pi, T, R)
%theta = zeros(R, 1);
%for i=1:R,
%   theta(i) = pi * T / i; 
%end
%end
% END THETA