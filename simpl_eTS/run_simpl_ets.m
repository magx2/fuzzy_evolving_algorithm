function run_simpl_ets
format long

    [x, y, r, OMEGA, opis] = prosta_funkcja();
    %wykonaj(x, y, r, OMEGA, opis);
    
    [x, y, r, OMEGA, opis] = nieruchomosci();
    wykonaj(x, y, r, OMEGA, opis);
    
    [x, y, r, OMEGA, opis] = sml();
    wykonaj(x, y, r, OMEGA, opis);
    
    [x, y, r, OMEGA, opis] = bike_data_day();
    wykonaj(x, y, r, OMEGA, opis);
    
    [x, y, r, OMEGA, opis] = bike_data_hour();
    wykonaj(x, y, r, OMEGA, opis);
end

function wykonaj( x, y, r, OMEGA, opis ) 
    disp([ 'Zaczynam [' opis{1} '] x '  num2str(length(r) * length(OMEGA)) ]);
    
    dane = opis{1};
    algorytm = 'simpl eTS';
    
    tismp = TimeStamp;
    
    csvFileName = ['G:\mgr\csv\' dane '\' dane '-' algorytm '-' tismp '.csv'];
    
    fid = fopen(csvFileName, 'w');
    fprintf(fid, 'nr,r,OMEGA,RMSE,R,czas\n');
    fclose(fid);
    
    wynikiFolder = [ tismp  '\'];
    mkdir(['G:\mgr\wyniki\' dane '\' wynikiFolder]);
    
    wynikiFileName = ['G:\mgr\wyniki\' dane '\' dane '-' algorytm '-' tismp '.csv'];
    fid = fopen(wynikiFileName, 'w');
    fprintf(fid, 'nr,y,y^,R\n');
    fclose(fid);
    
    k=1;
    for r_val=r,
        for OMEGA_val=OMEGA,
            
            wynikiFileName = ['G:\mgr\wyniki\' dane '\' wynikiFolder dane '-' algorytm '-' tismp '-' num2str(k) '.csv'];
            fid = fopen(wynikiFileName, 'w');
            fprintf(fid, 'y,y^,R\n');
            fclose(fid);
            
            tic;
            [ y_przewidywane, R_w_czasie, opis ] = simpl_ets( x, y, r_val, OMEGA_val, opis );
            czas = toc;
            RMSE = wyniki(cell2mat(y_przewidywane), y', R_w_czasie, opis, k);
            
            R = opis{3};
            wiersz = [ k r_val OMEGA_val RMSE R czas ];
            
            disp([ 'Koniec k=' num2str(k) ' R=' num2str(R) ' OMEGA=' num2str(OMEGA_val) ' r=' num2str(r_val)  ' RMSE='  num2str(RMSE) ]);
            
            k=k+1;
            
            %fid = fopen(csvFileName, 'a');
            dlmwrite(csvFileName, wiersz, '-append', 'precision', '%.6f', 'delimiter', ',');
            %fclose(fid);
            
            %fid = fopen(wynikiFileName, 'a');
            dlmwrite(wynikiFileName, [ y', cell2mat(y_przewidywane), R_w_czasie' ], '-append', 'precision', '%.6f', 'delimiter', ',');
            %fclose(fid);
        end
    end    

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
    
    r = [ 0.001 0.01 0.1 1 10 100 ];
    
    OMEGA = [ 0.001 0.01 0.1 1 10 100 ];
    
    opis = cell(1,1);
    opis{1} = 'y=2x+3x+4';
end

function [ x, y, r, OMEGA, opis ] = nieruchomosci() 
    filename = 'G:\mgr\dane\wart_nier_niezab_wg_czasu_07.csv';
 
     csv = csvread(filename,1,0);
     
     x = csv(:, 1:13)';
     y = csv(:, 14)';
    
    r = [ 0.001 0.01 0.1 1 10 100 ];
    
    OMEGA = [ 0.001 0.01 0.1 1 10 100 ];
    
    opis = cell(1,1);
    opis{1} = 'Nieruchomosci';
end

function [ x, y, r, OMEGA, opis ] = bike_data_day() 
    filename = 'G:\mgr\dane\BikeSharingDatasetDataSet\day.csv';
 
    csv = csvread(filename,1,0);
    
    x = csv(:, 3:13)';
    y = csv(:, 16)';
    
    r = [ 0.001 0.01 0.1 1 10 100 ];
    
    OMEGA = [ 0.001 0.01 0.1 1 10 100 ];
    
    opis = cell(1,1);
    opis{1} = 'Bike Sharing Data - day';
end

function [ x, y, r, OMEGA, opis ] = bike_data_hour() 
    filename = 'G:\mgr\dane\BikeSharingDatasetDataSet\hour.csv';
 
    csv = csvread(filename,1,0);
     
    x = csv(:, 3:13)';
    y = csv(:, 16)';
    
    r = [ 0.001 0.01 0.1 1 10 100 ];
    
    OMEGA = [ 0.001 0.01 0.1 1 10 100 ];
    
    opis = cell(1,1);
    opis{1} = 'Bike Sharing Data - hour';
end

function [ x, y, r, OMEGA, opis ] = sml() 
    filename = 'G:\mgr\dane\SML2010DataSet\NEW-DATA-1.T15.txt.csv';
 
    %csv = csvread(filename,1,0);
     
    %x1 = csv(:, 3:21);
    %y1 = csv(:, 22);
    
    filename = 'G:\mgr\dane\SML2010DataSet\NEW-DATA-2.T15.txt.csv';
 
    csv = csvread(filename,1,0);
    
    %x2 = csv(:, 3:21);
    %y2 = csv(:, 22);
    x = csv(:, 13:21)';
    y = csv(:, 22)';
    
    %x = [x1; x2]';
    %y = [y1;y2]';
    
    r =  fliplr([ 0.001 0.01 0.1 1 10 100 ]);
    
    OMEGA =  fliplr([ 0.001 0.01 0.1 1 10 100 ]);
    
    opis = cell(1,1);
    opis{1} = 'sml';
end