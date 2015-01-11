function [ y_daszek, R_w_czasie, opis, S, S_min, S_max, S_podmiana, S_nowy ] = simpl_ets( x, y, r, OMEGA, opis )

    debug=1;

    % x macierz [n x k] gdzie k liczba danych
    % y prawdziwe wyjscie jakie zaobserwowano [1 x k]
    % r is a positive constant, which defines the zone of influence of the i-th fuzzy rule antecedent
    opis{2} = 'simpl eTS';
    
    % psi(k) - [(n+1), R]
    % THETA(k) - [(n+1), mR]
    % y_daszek(k) - [m, 1]
    % C(k) - [m, m]
    n = length(x(:,1));
    m = length(y(:,1));
    K = length(x);
    R = 1;
    k = 1;
    
    x_gwiazdka = cell(R, 1);
    x_gwiazdka{R} =  x(:, 1);
    
    ilosc_w_klastrze = cell(R,1);
    ilosc_w_klastrze{R} = 1;
    
    % stage one - inicjalizacja parametrów
    z = cell(K, 1); % read z(1) 
    z{1} = [x(:,1)', y(:,1)']';
    z_gwiazdka = cell(R,1);
    z_gwiazdka{R} = z{1};
    
    z_sr = cell(K,1);
    z_sr{1} = zeros((m+n),1);
    sigma = cell(K,1);
    sigma{1} = zeros((m+n),1);
    
    created=0.3;
    
    S_gwiazdka = cell(K, R);
    S_gwiazdka{1, R} = 1;
                                                                            S = cell(1, K);
                                                                            S{1}=S_gwiazdka{1, R};
                                                                            S_min = cell(1, K);
                                                                            S_min{1}=S_gwiazdka{1, R};
                                                                            S_max = cell(1,K);
                                                                            S_max{1}=S_gwiazdka{1, R};
                                                                            S_podmiana= cell(1,K);
                                                                            S_podmiana{1}=0;
                                                                            S_nowy= cell(1,K);
                                                                            S_nowy{1}=1;
    
    
    % (15)
    PI = cell(K, 1);
    PI{1} = zeros((n+1), m);
    
    THETA = cell(K, R);
    THETA{1, R} = PI{1};
    psi = cell(K, R);
    psi{1, R} = zeros((n+1), 1);
        
    y_daszek = cell(K, 1);
    y_daszek{1} = zeros(1, m);
    
    C = cell(K, R);
    C{1, R} = OMEGA * eye(n+1);
    
    % w celach statystycznuch ilosc klastrow w czasie k
    R_w_czasie = zeros(1,K);
    R_w_czasie(1) = R;
 
    
                    x_k =x(:, 2);
                    x_k_e = [1 x_k']';
                     mi1 = mi(x_k, x_gwiazdka, r, R);
                        lambda1=cell(R,1);
                        lam_0=0;
                        lam_z_0=zeros(1,R);
                        for i=1:R,
                            lambda1{i} = lambda( mi1, i );

                            if ilosc_w_klastrze{i} / 2 < 0.01
                                lam_0=lam_0+lambda1{i};
                                lambda1{i}=0;
                                lam_z_0(i)=i;
                            end

                            psi{2, i} = lambda1{i} * x_k_e;
                        end
                        lam_0=lam_0/(R-nnz(lam_z_0));
                        if lam_0 > 0
                            % dodawanie z zer
                            for i=1:R,
                                if ~ismember([i],lam_z_0)
                                    lambda1{i}=+lambda1{i}+1;
                                end
                            end
                        end
    
    while k <= K - 1
        
        for i=1:R,
           psi{k+1, i} =  psi{k, i}; 
        end
        % stage 2
        x_k = x(:, k + 1); % k-ta kolumna; k-ta dana
        
        x_k_e = [1 x_k']'; % [1, x_k]
        
        % stage 3
        % Estimate the output: y_daszek(k + 1) by (16)
        % (16) y_daszek(k + 1) = psi^T(k + 1)THETA(k)
       
        y_daszek{k + 1} = zeros(1, m);
        for i=1:R,
            y_daszek{k + 1} = y_daszek{k + 1} + psi{k + 1, i}' * THETA{k, i};
        end
        
        % stage 4
        k = k + 1;
        % Read y(k)
        y_k = y(:, k);
        
        % stage 5
        % Calculate S_k recursively by (19) - (20)
        z{k} = [x_k', y_k']';
        
        z_sr{k} =  ( (k-1)/k ) * z_sr{k-1} + z{k}/k; % (11)
        sigma{k} = (k-1)/k * sigma{k-1} + 1/(k-1) * (z{k}-z_sr{k}).^2; %( 12)
        
        for j=1:(m+n),
            z{k}(j) = ( z{k}(j)-z_sr{k}(j) ) / sigma{k}(j); % (13)
        end

        S_k = fS_k(k, n, m, z);
        
        S_k = S_ang(k, n, m, z);
     
        % stage 6
        % update S(z_i_gwiazdka) each centre by 23a
%         z_tmp = sumsqr(z{k}-z{k-1});
%         for i=1:R,
%             S_gwiazdka{k, i} = (k-2)/(k-1) * S_gwiazdka{k-1, i} + z_tmp;
%         end;
        z_tmp=sumsqr(z{k}-z{k-1});
        for i=1:R,
%             S_gwiazdka{k, i} = (k-2)/(k-1) * S_gwiazdka{k-1, i} + sumsqr(z{k}-z_gwiazdka{i});
%             S_gwiazdka{k, i} = (k-2)/(k-1) * S_gwiazdka{k-1, i} + sumsqr(z{k}-z{k-1});
            a = (k-1)*S_gwiazdka{k-1, i};
            b = k-2+S_gwiazdka{k-1, i}+ ( S_gwiazdka{k-1, i} * z_tmp );
            S_gwiazdka{k,i} = a/b;
        end;
    
        % ang(21)
        for i=1:R,
            a = (k-1)*S_gwiazdka{k-1, i};
            b = k - 2 + S_gwiazdka{k-1, i} + (S_gwiazdka{k-1, i} * z_tmp);
            S_gwiazdka{k,i}=a/b;
        end

                                                                            s= [S_gwiazdka{k, :}];
                                                                            S_min{k} = min(s);
                                                                            S_max{k} = max(s);
                                                                            S{k} = S_k;
        
        % Compare S(k) with S_k(z_i_gw)
        % stage 7
        retu_24 = s_ets_24( S_k, [S_gwiazdka{k, :}]);
        [ retu_25, l ] = s_ets_25( x_k, x_gwiazdka, R, r);
        byl_dodany_klaster=0;
        
        if retu_24 && retu_25 
% podmien klaster
            if debug
                disp([ 'podmien klaster l=', num2str(l), '   k=', num2str(k) 'R=', num2str(R),  ' y-y^=' num2str(y_k-y_daszek{k}) ]);
            end
           
           x_gwiazdka{l} = x_k;
           S_gwiazdka{k, l} = S_k;
           z_gwiazdka{l} = z{k};
           
           ilosc_w_klastrze{l} = ilosc_w_klastrze{l} +1;    

                                                                            S_podmiana{k} = S_k;
                                                                            S_nowy{k} = -created;
% END podmien klaster
        elseif retu_24 
% nowy klaster
            if debug
              disp(['nowy    klaster R=', num2str(R+1), '   k=', num2str(k) ' y-y^=' num2str(y_k-y_daszek{k}) ]);
            end
            % (23)
            
            R = R + 1;
            x_gwiazdka{R} =  x_k;
            S_gwiazdka{k, R} = S_k;
            ilosc_w_klastrze{R} = 1;
                                                                            S_podmiana{k} = -created;
                                                                            S_nowy{k} = S_k;
            
            z_gwiazdka{R} = z{k};
            
            
            
%                         mi1 = mi(x_k, x_gwiazdka, r, R);
%                         lambda1=cell(R,1);
%                         lam_0=0;
%                         lam_z_0=zeros(1,R);
%                         for i=1:R,
%                             lambda1{i} = lambda( mi1, i );
% 
%                             if ilosc_w_klastrze{i} / k < 0.01
%                                 lam_0=lam_0+lambda1{i};
%                                 lambda1{i}=0;
%                                 lam_z_0(i)=i;
%                             end
% 
%                             psi{k, i} = lambda1{i} * x_k_e;
%                         end
%                         lam_0=lam_0/(R-nnz(lam_z_0));
%                         if lam_0 > 0
%                             % dodawanie z zer
%                             for i=1:R,
%                                 if ~ismember([i],lam_z_0)
%                                     lambda1{i}=+lambda1{i}+1;
%                                 end
%                             end
%                         end
            
                        mi1 = mi(x_k, x_gwiazdka, r, R);
                        lambda1=cell(R,1);
                        for i=1:R,
                            lambda1{i} = lambda( mi1, i );


                            psi{k, i} = lambda1{i} * x_k_e;
                        end
                       
            
            
            % dla matlaba powikszam tablice  
            C{k, R} = OMEGA * eye(n+1);  
            C{k-1, R} = C{k, R};
            
            suma=0;
            for i=1:(R-1),
                suma=suma+THETA{k-1, i}*lambda1{i};
            end
            THETA{k, R} = suma;% / sum(cell2mat(lambda1)) ;
            THETA{k-1, R} = THETA{k, R};
            
            psi{k, R} = zeros((n+1), 1);
            
            
%             for i=1:(R-1),
%                 C{k, i} = C{k - 1, i} - ( C{k - 1, i} * psi{k,i} * psi{k,i}' * C{k - 1, i} ) / ( 1 + psi{k,i}' * C{k - 1, i} * psi{k,i});
%                 THETA{k, i} = THETA{k - 1, i} + C{k} * psi{k, i} * ( y_k - psi{k, i}' * THETA{k - 1, i} );
%             end
%             for i=1:(R-1),
%                 C{k, i} = C{k - 1, i};
%                 THETA{k, i} = THETA{k - 1, i};
%             end
            byl_dodany_klaster=1;
% END nowy klaster
        else
% nic nie robienie z klastrem
            if debug
              disp(['nie rob nic R=', num2str(R), '24=',int2str(retu_24), ' 25=',int2str(retu_25), ' k=', num2str(k) ' y-y^=' num2str(y_k-y_daszek{k}) ]);
            end
             
            % dodaje ilosc do klastra
            min_idx = 1;
            min_val = ( sumsqr((x_k - x_gwiazdka{min_idx})) );
            for i=1:R,
                val = ( sumsqr((x_k - x_gwiazdka{i})) );
                if val < min_val
                    min_val=val;
                    min_idx = i;
                end
            end
            ilosc_w_klastrze{min_idx} = ilosc_w_klastrze{min_idx} + 1;
            
             for i=1:(R),
                C{k, i} = C{k - 1, i} ;
                THETA{k, i} = THETA{k - 1, i};
             end
            
             
            S_podmiana{k} = -created;
            S_nowy{k} = -created;
% END nic nie robienie z klastrem
        end

       for i=1:(R),
            C{k, i} = C{k - 1, i} - ( C{k - 1, i} * psi{k,i} * psi{k,i}' * C{k - 1, i} ) / ( 1 + psi{k,i}' * C{k - 1, i} * psi{k,i});
            THETA{k, i} = THETA{k - 1, i} + C{k} * psi{k, i} * ( y_k - y_daszek{k} );
       end
        
                                                                            R_w_czasie(k) = R;


                                                                            %s= [S_gwiazdka{k, :}];
                                                                            %S_min{k} = min(s);
                                                                            %S_max{k} = max(s);
                                                                            %S{k} = S_k;
    end
    
    opis{3} = R;
    opis{4} = r;
    opis{5} = OMEGA;
    %opis{6} = z_gwiazdka;
    opis{7} = x;
end

function [ retu_25, l ] = s_ets_25( x_k, x_gwiazdka, R, r)
    min = sumsqr(x_k - x_gwiazdka{1});
    l = 1;
    for i=2:R,
        new_min = sumsqr(x_k - x_gwiazdka{i});
        if new_min < min
            min = new_min;
            l = i;
        end
    end 
    
    retu_25 = min < ( r / 2 );
end

function [ retu_24 ] = s_ets_24( S_k, S_k_gwiazdka)
    maxx = max(S_k_gwiazdka);
    minn = min(S_k_gwiazdka);
    
    retu_24 = S_k <= minn || S_k >= maxx;
%     retu_24 = S_k >= maxx;
end

function [ S_k ] = fS_k(k, n, m, z )
    % (19)-(20)
    a = ( (k - 1)*(n + m) );
    
    new_gamma = gamma_k(z, k);
    
    new_z_k = 0;
    for j=1:n+m,
        beta_j =  beta_j_k(z, k, j);
        new_z_k = new_z_k + z{k}(j) * beta_j;
    end
    
    b = (k-1) * sumsqr(z{k}) - 2 * new_z_k + new_gamma;

    S_k = b / a;
end

function [ S_k ] = globalS(k, n, m, z )
    % (17)
    suma=0;
    for l=1:k-1,
        suma=suma+sumsqr(z{l} - z{k});
    end
    
    S_k = suma / ( (k-1)*(n+m) );
end

function [ S_k ] = S_ang(k, n, m, z )
    % (19)
    suma=0;
    for l=1:k-1,
        suma=suma+sumsqr(z{l} - z{k});
    end
    
    S_k = 1 / ( 1 + (suma/(k-1)) );
end

function [ beta_j_k ] = beta_j_k( z, k, j )
    sum = 0;
    for l=1:(k-1),
        sum = sum + z{l}(j);
    end
    beta_j_k = sum;
end

function [ gamma_k ] = gamma_k( z, k )
    sum = 0;
    for l=1:(k-1),
        sum = sum + sumsqr(z{l});
    end
    gamma_k = sum;
end

function [ lambda ] = lambda( mi, i )
    % lamda - skalar
    lambda = mi{i} / sum([mi{:}]);
end

function [ mi1 ] = mi(x, x_gwiazdka, r, R)
    % mi - skalar
    mi1 = cell(R, 1);
    for i=1:R,
        out = 1;
        for j=1:length(x),
            out = out * ( 1 + ( 2 * (x(j) - x_gwiazdka{i}(j))  / r )^2);
        end
        mi1{i} = 1 / out;
    end
end