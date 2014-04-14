function [ p ] = angelov_function_20( k, z )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    mianownik = ( (k - 1) + ( ipsilon(k, z) + 1 ) + sigma(k, z) - 2 * v(k) );
    p = ( k - 1 ) / mianownik;

end

function [ out ] = ipsilon( k, z )
    z_k = z(k,:);
    out = 0;
    for j=1:length(z_k),
        out = out + z_k(k, j)^2;
    end
end

function [ out ] = sigma( k, z )
    out = 0;
    for l=1:k-1,
        for j=1:length(z),
            out = out + z(l, k);
        end
    end
end

function [ out ] = v( k, z )
    z_k = z(k,:);
    out = 0;
    for j=1:length(z),
        out = out + z_k(j) * beta(k, j, z);
    end
end

function [ out ] = beta( k, j, z )
    out = 0;
    for l=1:k - 1,
        out = out + z(l, j);
    end
end