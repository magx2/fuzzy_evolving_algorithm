% dane wej?ciowe
liczba_danych = 5
dlugosc_danych = 10

% przygotowanie danych matlaba


% inicjalizcacja algorytmu
k = 1;
R = 1;
theta = 0;
pi = 0;
x_gwiazdka = wczytaj_porcje_danych(k);

% algorytm
for i=1:liczba_danych,
        k = k + 1;
        
        % wczytaj porcje danych
        x_k = wczytaj_porcje_danych(k);
        y = 
        
        % policz rekursywnie z (20)
        % w tym miejscu jest liczony wczytanych danych
        p = angelov_function_20(k, z);
        
        % odswiez P_k(z*l) wed?ug (21)
        for l=1:R, 
            z_gwiazdka(l) = angelov_function_21(k, z_teraz);
        end
         
        % sprawdzamy warunek czy nowa regula jest lepsza od tych co mamy
        jestLepsza = angelov_function_nowa_regula_lepsza(p, z_gwiazdka);
        
        if p > z_gwiazdka(R) && angelov_function_22(0,0,0,0,0,0,0)
            
        elseif p > z_gwiazdka(R)
            % ADD bew rule
            R = R + 1;
            z_gwiazdka = [z_gwiazdka; p];
            % !!! ??? x_gwiazdka = [x_gwiazdka; x_k?!];
            % TODO reset by (27)-(28) or (32)
        end
        
        % rekursywnie odswierz parametry poprzez RLS (24)-(26) 
        % lub poprze wRLS (29)-(31)
end
