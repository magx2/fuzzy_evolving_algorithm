function [ p ] = angelov_function_21( k, z )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    if k > 1
        licznik = ( k - 1 ) * angelov_function_21(k - 1, z );
        mianownik = ( k - 2 ) + angelov_function_21(k - 1, z ) + angelov_function_21(k - 1, z) * sum(k, z);
    
        p = licznik / mianownik;
    else 
        p = 0;
    end
end

function [ out ] = sum( k, z )
    out = 0;
    for j=1:length(z)
       out = out + ( z(k, j) - z(k - 1, j) )^2;
    end
end