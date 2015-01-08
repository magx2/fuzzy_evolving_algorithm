function [ y_daszek, R_w_czasie, opis, S, S_min, S_max, S_podmiana, S_nowy ] = ets( x, y, r, OMEGA, opis )

    debug=1;

    % x macierz [n x k] gdzie k liczba danych
    % y prawdziwe wyjscie jakie zaobserwowano [1 x k]
    % r is a positive constant, which defines the zone of influence of the i-th fuzzy rule antecedent
    opis{2} = 'eTS';
    
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
    
    % stage one - inicjalizacja parametrów
    z = cell(K, 1); % read z(1) 
    z{1} = [x(:,1)', y(:,1)']';
    z_gwiazdka = cell(R,1);
    z_gwiazdka{R} = z{1};
        
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
    
    
    
    THETA = cell(K, R);
    psi = cell(K, R);
    psi{1, R} = zeros((n+1), 1);
        
    y_daszek = cell(K, 1);
    y_daszek{1} = zeros(1, m);
    y_daszek{2} = zeros(1, m);
    
    C = cell(K, R);
    C{1, R} = OMEGA * eye(n+1);
    
    THETA{1, R} = C{1, R}*[1 x(:, 1)']'*(y(:, 1));
    
    % w celach statystycznuch ilosc klastrow w czasie k
    R_w_czasie = zeros(1,K);
    R_w_czasie(1) = R;
 
    sigma_k = 0;
    beta_k = 0;
    
    while k <= K - 1
        k=k+1;
       
        % stage 2
        x_k = x(:, k); % k-ta kolumna; k-ta dana
        x_k_e = [1 x_k']'; % [1, x_k]
        y_k = y(:, k);
        
        % read z_k
        z{k} = [x_k', y_k]';
        z_k = z{k};
        
        % stage 3
        % calck P_k(z_k) by (20)
%         theta_k = sumsqr(z_k);
%         sigma_k = sigma_k + sumsqr(z{k-1});
%         beta_k = beta_k - z{k-1};
%         v_k = 0;
%         for j=1:(n+1),
%             v_k = v_k + z_k(j) * beta_k(j);
%         end
%         S_k = (k-1) / ( (k-1) * (theta_k + 1) + sigma_k - 2 * v_k );
        % (19)
        suma=0;
        for l=1:k,
           suma = suma + sumsqr(z{l}-z_k); 
        end
        S_k = 1 / ( 1 + 1/(k-1) * suma);
        
        % stage 4
        % Up-date P_k(z*) (21)
        z_tmp = sumsqr(z{k}-z{k-1})*4.5;
        for i=1:R,
            a = (k-1)*S_gwiazdka{k-1, i};
            b = k - 2 + S_gwiazdka{k-1, i} + (S_gwiazdka{k-1, i} * z_tmp);
            S_gwiazdka{k,i}=a/b;
%             tmp = [ z_tmp a b S_gwiazdka{k,i} S_k]
        end

                                                                            % statystyka
                                                                            s= [S_gwiazdka{k, :}];
                                                                            S_min{k} = min(s);
                                                                            S_max{k} = max(s);
                                                                            S{k} = S_k;
                                                                            % END statystyka
        % wylicz psi
        lambda=zeros(1,R);
        for i=1:R,
           lambda(i) = exp(-sumsqr(x_k - x_gwiazdka{i}));
           psi{k,i} =  x_k_e * lambda(i);
        end
        % recursily update RLS by (24)-(26)
        for i=1:R,
            % (25)
            C{k, i} = C{k-1, i} - ( C{k-1, i} * psi{k-1,i} * psi{k-1,i}' * C{k-1,i} )/( 1 +  psi{k-1,i}' * C{k-1,i} * psi{k-1,i});
            % (24)
            THETA{k, i} = THETA{k-1, i} + C{k, i}*psi{k-1,i}*(y_k - y_daszek{k});
        end
                                                                            
        % stage 5
        % Compare S(k) with S_k(z_i_gw)
        maxSg = max([S_gwiazdka{k, :}]);
        ifP = S_k > maxSg;
        
        minn = sumsqr(z_k - z_gwiazdka{1});
        l = 1;
        for i=2:R,
            new_min = sumsqr(z_k - z_gwiazdka{i});
            if new_min < minn
                minn = new_min;
                l = i;
            end
        end 
        if22 = (S_k / maxSg) - (minn / r) >= 1;
        
        if ifP && if22
% podmien klaster
            if debug
                disp([ 'podmien klaster l=', num2str(l), '   k=', num2str(k) ' y-y^=' num2str(y_k-y_daszek{k}) ]);
            end
           
            z_gwiazdka{l} = z_k;
            x_gwiazdka{l} = x_k;
            S_gwiazdka{k, l} = S_k;
            % PI?!
            % C?!
            % wylicz psi
            lambda=zeros(1,R);
            for i=1:R,
               lambda(i) = exp(-sumsqr(x_k - x_gwiazdka{i}));
               psi{k,i} =  x_k_e * lambda(i);
            end
            % recursily update RLS by (24)-(26)
            for i=1:R,
                % (25)
                C{k, i} = C{k-1, i} - ( C{k-1, i} * psi{k-1,i} * psi{k-1,i}' * C{k-1,i} )/( 1 +  psi{k-1,i}' * C{k-1,i} * psi{k-1,i});
                % (24)
                THETA{k, i} = THETA{k-1, i} + C{k, i}*psi{k-1,i}*(y_k - y_daszek{k});
            end

                                                                            S_podmiana{k} = S_k;
                                                                            S_nowy{k} = -created;
% END podmien klaster
        elseif ifP
% nowy klaster
            if debug
              disp(['nowy    klaster R=', num2str(R+1), '   k=', num2str(k) ' y-y^=' num2str(y_k-y_daszek{k}) ]);
            end
            
            R = R + 1;
            
            z_gwiazdka{R} = z{k};
            x_gwiazdka{R} =  x_k;
            S_gwiazdka{k, R} = S_k;
            
            % reste by (32)
            C{k, R} = OMEGA * eye(n+1);
            THETA{k, R}=zeros(n+1,1);
            for i=1:(R-1),
               THETA{k, R} =  THETA{k, R} + THETA{k-1, i}*lambda(i);
            end
            % wylicz psi
%             lambda=zeros(1,R);
%             for i=1:R,
%                lambda(i) = exp(-sumsqr(x_k - x_gwiazdka{i}));
%                psi{k,i} =  x_k_e * lambda(i);
%             end
            psi{k,R}=zeros(1,(n+1))';
                % DO SPRAWDZENIA JAKI WYMIAR TABLICY!!!
            % consequence parameters and co-variance matrices are
            % reset by (27)–(28) or (32), respectively, for the global or
            % local estimation.
            
                                                                            S_podmiana{k} = -created;
                                                                            S_nowy{k} = S_k;
% END nowy klaster
        else
% nic nie robienie z klastrem
            if debug
              disp(['nie rob nic ifP=',int2str(ifP), ' if22=',int2str(if22), ' k=', num2str(k) ' y-y^=' num2str(y_k-y_daszek{k}) ]);
            end
                                                                            S_podmiana{k} = -created;   
                                                                            S_nowy{k} = -created;
% END nic nie robienie z klastrem
        end
        
        % predict next output y(k+1) from (23)
        y_daszek{k+1} = 0;
        for i=1:R,
            y_daszek{k+1} = y_daszek{k+1} + psi{k, i}' * THETA{k, i};
        end
        
                                                                            R_w_czasie(k) = R;
        
                                                                            s= [S_gwiazdka{k, :}];
                                                                            S_min{k} = min(s);
                                                                            S_max{k} = max(s);
                                                                            S{k} = S_k;
    end
    y_daszek=y_daszek(1:(end-1),1);
    
    opis{3} = R;
    opis{4} = r;
    opis{5} = OMEGA;
    %opis{6} = z_gwiazdka;
    opis{7} = x;
end