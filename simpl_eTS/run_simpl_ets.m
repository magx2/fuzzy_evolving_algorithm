function run_simpl_ets
format long

%     [x, y, r, OMEGA, opis] = prosta_funkcja();
%      %r=1000;
%     wykonaj(x, y, r, OMEGA, opis);
%     
%     [x, y, r, OMEGA, opis] = gas();
%     wykonaj(x, y, r, OMEGA, opis);
%     
%     [x, y, r, OMEGA, opis] = nieruchomosci2();
%     wykonaj(x, y, r, OMEGA, opis);
%     
    [x, y, r, OMEGA, opis] = nieruchomosci();
    wykonaj(x, y, r, OMEGA, opis);
% %     
%     [x, y, r, OMEGA, opis] = sml();
%     wykonaj(x, y, r, OMEGA, opis);
%     
%     [x, y, r, OMEGA, opis] = bike_data_day();
%     wykonaj(x, y, r, OMEGA, opis);
     
%     [x, y, r, OMEGA, opis] = bike_data_hour();
%     wykonaj(x, y, r, OMEGA, opis);
end

function [ r, OMEGA ] = wejscie()
%     r = [ 20000 50000 10000 5000 100000  ];
%     OMEGA = [ 2000000 20 2000 10000 100000 ] * 1000;
    r = [20000, 40000];
    OMEGA = 2000000;
end

function wykonaj( x, y, r, OMEGA, opis ) 
    disp([ 'Zaczynam [' opis{1} '] x '  num2str(length(r) * length(OMEGA)) ]);
   
    times = 3;
    
    dane = opis{1};
    algorytm = 'simpl eTS';
    
    tismp = TimeStamp;
    
    csvFileName = ['G:\mgr\csv\' dane '\' dane '-' algorytm '-' tismp '.csv'];
    
    fid = fopen(csvFileName, 'w');
    fprintf(fid, 'nr,r,OMEGA,RMSE,R,czas\n');
    fclose(fid);
    
    wynikiFolder = [ tismp  '\'];
    mkdir(['G:\mgr\wyniki\' dane '\' wynikiFolder]);
    
    RMSE_all = zeros(length(r)*length(OMEGA),1);
    
    k=1;
    for r_val=r,
        for OMEGA_val=OMEGA,
            
            wynikiFileName = ['G:\mgr\wyniki\' dane '\' wynikiFolder dane '-' algorytm '-' tismp '-' num2str(k) '.csv'];
            fid = fopen(wynikiFileName, 'w');
            fprintf(fid, 'y,y^,y-y^,R\n');
            fclose(fid);
            
            tic;
            [ y_przewidywane, R_w_czasie, opis, S, S_min, S_max, S_podmiana, S_nowy ] = simpl_ets( x, y, r_val, OMEGA_val, opis );
            czas = toc;
             
%             fi = length(y)/times;
%             
%             y = y((times-1)*fi+1:end);
%             y_przewidywane = y_przewidywane((times-1)*fi+1:end);
%             R_w_czasie = R_w_czasie(1:fi);
%             
%             S = S(1:fi); 
%             S_min = S_min(1:fi); 
%             S_max = S_max(1:fi); 
%             S_podmiana = S_podmiana(1:fi); 
%             S_nowy  = S_nowy(1:fi); 
            
%             y = y(2:end);
%             y_przewidywane = y_przewidywane(2:end);
%             R_w_czasie = R_w_czasie(2:end);
            
            RMSE = wyniki(cell2mat(y_przewidywane), y', R_w_czasie, opis, k , S, S_min, S_max, S_podmiana, S_nowy);
            RMSE_all(k) = RMSE;
            
            
            R = opis{3};
            wiersz = [ k r_val OMEGA_val RMSE R czas ];
            
            disp([ 'Koniec k=' num2str(k) ' R=' num2str(R) ' czas=' num2str(czas/60) ' (OMEGA|r)=(' num2str(OMEGA_val) '|' num2str(r_val)  ') RMSE='  num2str(RMSE) ]);
            
            k=k+1;
            
            %fid = fopen(csvFileName, 'a');
            dlmwrite(csvFileName, wiersz, '-append', 'precision', '%.6f', 'delimiter', ',');
            %fclose(fid);
            
            %fid = fopen(wynikiFileName, 'a');
            dlmwrite(wynikiFileName, [ y', cell2mat(y_przewidywane),y'-cell2mat(y_przewidywane), R_w_czasie' ], '-append', 'precision', '%.6f', 'delimiter', ',');
            %fclose(fid);
        end
    end
    
    Wyniki = RMSE_all
    najmniejsza = min(RMSE_all)
end

function [ x, y, r, OMEGA, opis ] = nieruchomosci() 
    filename = 'G:\mgr\dane\Datastream 1998-2011 uw 4 zm Area Storeys Centre Age PriceTrans.csv';
 
     csv = csvread(filename,1,0);
     
     x = csv(:, 1:4)';
     y = csv(:, 5)';
     
     x=x(:,1:1000);
     y=y(:,1:1000);
     
    [r, OMEGA] = wejscie();
    
    opis = cell(1,1);
    opis{1} = 'Nieruchomosci';
end













function [ x, y, r, OMEGA, opis ] = prosta_funkcja() 
    % y = 2x_1 + 3x_2 + 4
    
    K=250;
    for k=1:K,
        x(:, k) = [ k k*2 ]';
    end
    
    y = zeros(1, length(x));
    for k=1:length(x),
        x_k = x(:, k);
        y(k) =  x_k(1) * 2 + x_k(2) * 3 + 4;
    end
    
    
    x = [  x x x];
    y = [  y y y];
    
    [r, OMEGA] = wejscie();
    
    opis = cell(1,1);
    opis{1} = 'y=2x+3x+4';
end

function [ x, y, r, OMEGA, opis ] = sinus() 
    % y = 2x_1 + 3x_2 + 4
    
    K=50;
    for k=1:K,
        x(:, k) = [ k ]';
    end
    
    y = zeros(1, length(x));
    for k=1:length(x),
        x_k = x(:, k);
        y(k) =  sin(x_k(1));
    end
    
    [r, OMEGA] = wejscie();
    
    opis = cell(1,1);
    opis{1} = 'y=2x+3x+4';
end

function [ x, y, r, OMEGA, opis ] = nieruchomosci_stare() 
    filename = 'G:\mgr\dane\wart_nier_niezab_wg_czasu_07.csv';
 
     csv = csvread(filename,1,0);
     
     x = csv(:, 1:13)';
     y = csv(:, 14)';
    
    x = [  x x x];
    y = [  y y y];
     
    [r, OMEGA] = wejscie();
    
    opis = cell(1,1);
    opis{1} = 'Nieruchomosci';
end

function [ x, y, r, OMEGA, opis ] = nieruchomosci2() 
    filename = 'G:\mgr\dane\wart_nier_niezab_wg_czasu_07-2.csv';
 
     csv = csvread(filename,1,0);
     
     x = csv(:, 1:13)'*1000;
     y = csv(:, 14)'*1000;
    
     
    x = [  x x x];
    y = [  y y y];
     
    [r, OMEGA] = wejscie();
    
    opis = cell(1,1);
    opis{1} = 'Nieruchomosci2';
end

function [ x, y, r, OMEGA, opis ] = bike_data_day() 
    filename = 'G:\mgr\dane\BikeSharingDatasetDataSet\day.csv';
 
    csv = csvread(filename,1,0);
    
    x = mat2gray(csv(:, 3:13)');
    y = mat2gray(csv(:, 16)');
    
    [r, OMEGA] = wejscie();
    
    opis = cell(1,1);
    opis{1} = 'Bike Sharing Data - day';
end

function [ x, y, r, OMEGA, opis ] = bike_data_hour() 
    filename = 'G:\mgr\dane\BikeSharingDatasetDataSet\hour.csv';
 
    csv = csvread(filename,1,0);
     
    x =  mat2gray(csv(:, 3:13)');
    y =  mat2gray(csv(:, 16)');
    
    [r, OMEGA] = wejscie();
    
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
    
    [r, OMEGA] = wejscie();
    
    opis = cell(1,1);
    opis{1} = 'sml';
end

function [ x, y, r, OMEGA, opis ] = gas()     
    filename = 'G:\mgr\dane\gas-furnace.csv';
 
    csv = csvread(filename,0,0);

%     x = mat2gray(csv(:, 1)');
%     y = mat2gray(csv(:, 2)');
    x = (csv(:, 1)');
    y = (csv(:, 2)');
    
    x = [  x x x];
    y = [  y y y];
    
%     x=x(1,1:50);
%     y=y(1,1:50);
    
    [r, OMEGA] = wejscie();
    
    opis = cell(1,1);
    opis{1} = 'gas';
end