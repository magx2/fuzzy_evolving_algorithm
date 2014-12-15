function [ out ] = run_eMG() 
dlugosc_wektora = 2;
ilosc_danych = 8;
omega = 50;
SIGMA_init = eye(dlugosc_wektora)
x = rand(ilosc_danych, dlugosc_wektora)
y_d = rand(dlugosc_wektora, 1)

out = eMG(x, y_d, 0, omega, SIGMA_init);
end

