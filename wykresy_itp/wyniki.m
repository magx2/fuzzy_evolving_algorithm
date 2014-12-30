function [ RMSE ] = wyniki( y_przewidywane, y_prawdziwe, R_w_czasie, opis )
% opis = [ nazwa_danych algorytm R r OMEGA ]
    dane = opis{1};
    algorytm = opis{2};
    R = opis{3};
    r = opis{4};
    OMEGA = opis{5};

    K = length(y_prawdziwe);

    out = [ y_przewidywane, y_prawdziwe, sqrt((y_przewidywane - y_prawdziwe).^2)];
    
    RMSE = sum(out(:,3)) / K;
    
    filename = ['G:\mgr\wykresy\' dane '\' dane '-' algorytm '-' TimeStamp];
    
    %plot1 = figure;
    plot([1:K], out(:,1), 'g',[1:K], out(:,2), '--b',[1:K], out(:,3), 'r');
    set(gcf, 'visible','off');
    title(['Wykres wartosci przewidywanej, prawdziwej i RMSE (' dane ' - ' algorytm ')']);
    xlabel('nr danej przychodzacej');
    ylabel('Wyniki');
    legend('y przewidywane','y', 'RMSE');
    
    MyBox = text(0.5,0.5,['RMSE=', num2str(RMSE') char(10) 'R=' num2str(R) char(10) 'r=' num2str(r) char(10) '\Omega=' num2str(OMEGA)]);
    set(MyBox,'Units','Normalized');
    set(MyBox,'Position',[0.02,0.9]);
    
    set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
    saveas(gcf, filename, 'png') %Save figure
    
    %plot([1:K], out(:,1), 'g',[1:K], out(:,2), '--b'); % wykres bez bledow
end