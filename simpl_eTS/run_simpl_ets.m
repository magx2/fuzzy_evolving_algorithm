function run_simpl_ets
format long

    [x, y, r, OMEGA, opis] = prosta_funkcja();
    wykonaj(x, y, r, OMEGA, opis);
    
    [x, y, r, OMEGA, opis] = nieruchomosci();
    wykonaj(x, y, r, OMEGA, opis);
end

function wykonaj( x, y, r, OMEGA, opis ) 
    k=1;
    wiersze = zeros(1, 4);
    for r_val=r,
        for OMEGA_val=OMEGA,
            [ y_przewidywane, R_w_czasie, opis ] = simpl_ets( x, y, r_val, OMEGA_val, opis );
            RMSE = wyniki(cell2mat(y_przewidywane), y', R_w_czasie, opis);
        
            R = opis{3};
            wiersze(k, :) = [ r_val OMEGA_val RMSE R ];
            
            k=k+1;
        end
    end
    
    dane = opis{1};
    algorytm = opis{2};
    
    csvFileName = ['G:\mgr\csv\' dane '\' dane '-' algorytm '-' TimeStamp '.csv'];
    
    fid = fopen(csvFileName, 'w');
    fprintf(fid, 'r,OMEGA,RMSE,R\n');
    fclose(fid);
    
    %csvwrite(csvFileName,wiersze);
    dlmwrite(csvFileName, wiersze, '-append', 'precision', '%.6f', 'delimiter', ',');
end

function [ x, y, r, OMEGA, opis ] = prosta_funkcja() 
    % y = 2x_1 + 3x_2 + 4
    
    K=50;
    for k=1:K,
        x(:, k) = [ k k*2 ]';
    end
    
    y = zeros(1, length(x));
    for k=1:length(x),
        x_k = x(:, k);
        y(k) =  x_k(1) * 2 + x_k(2) * 3 + 4;
    end
    
    r = [ 1 20 ];
    
    OMEGA = [ 1 2 5 7 11 13 ];
    
    opis = cell(1,1);
    opis{1} = 'y=2x+3x+4';
end

function [ x, y, r, OMEGA, opis ] = nieruchomosci() 
    filename = 'G:\mgr\dane\wart_nier_niezab_wg_czasu_07.csv';
 
     csv = csvread(filename,1,0);
     
     x = csv(:, 1:13)';
     y = csv(:, 14)';
    
    r = [ 0.1 0.5 1 20 ];
    
    OMEGA = [ 1 2 5 7 11 13 ];
    
    opis = cell(1,1);
    opis{1} = 'Nieruchomosci';
end