function [ RMSE ] = wyniki( y_przewidywane, y_prawdziwe, R_w_czasie, opis, nr, S, S_min, S_max, S_podmiana, S_nowy )
% opis = [ nazwa_danych algorytm R r OMEGA z_gwiazdka x ]
    dane = opis{1};
    algorytm = opis{2};
    R = opis{3};
    r = opis{4};
    OMEGA = opis{5};

    K = length(y_prawdziwe);

    out = [ y_przewidywane, y_prawdziwe, (y_przewidywane - y_prawdziwe)];
    
    RMSE = sqrt(mean(out(:,3).^2));
    
    filename = ['G:\mgr\wykresy\' dane '\' dane '-' algorytm '-' TimeStamp '-' num2str(nr)];
    
    %plot1 = figure;
    %plot([0:K-1], out(:,1), 'g',[0:K-1], out(:,2), '--b',[0:K-1], out(:,3), 'r');
    plot([0:K-1], out(:,1), 'g',[0:K-1], out(:,2), '--b');
    set(gcf, 'visible','off');
    title(['Wykres wartosci przewidywanej, prawdziwej i RMSE (' dane ' - ' algorytm ')']);
    xlabel('nr danej przychodzacej');
    ylabel('Wyniki');
    %legend('y przewidywane','y', 'RMSE');
    legend('y przewidywane','y');
    
    MyBox = text(0.5,0.5,['algorytm=', algorytm, char(10), 'RMSE=', num2str(RMSE') char(10) 'R=' num2str(R) char(10) 'r=' num2str(r) char(10) '\Omega=' num2str(OMEGA)]);
    set(MyBox,'Units','Normalized');
    set(MyBox,'Position',[0.02,0.9]);
    
    set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
    saveas(gcf, filename, 'png') %Save figure
    
    
    % R w czasie
    filename = ['G:\mgr\wykresy\' dane '\' dane '-' algorytm '-' TimeStamp 'R-w-czasie' '-' num2str(nr)];
    
    %plot1 = figure;
    %plot([0:K-1], out(:,1), 'g',[0:K-1], out(:,2), '--b',[0:K-1], out(:,3), 'r');
    plot([0:K-1], R_w_czasie, 'g');
    set(gcf, 'visible','off');
    title(['Liczba klastórw w czasie (' dane ' - ' algorytm ')']);
    xlabel('nr danej przychodzacej');
    ylabel('Liczba klastrów');
    %legend('y przewidywane','y', 'RMSE');
    legend('R');
    
    MyBox = text(0.5,0.5,['algorytm=', algorytm]);
    set(MyBox,'Units','Normalized');
    set(MyBox,'Position',[0.02,0.9]);
    
    set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
    saveas(gcf, filename, 'png') %Save figure
    
    
    % roznica
    filename = ['G:\mgr\wykresy\' dane '\' dane '-' algorytm '-' TimeStamp '-roznica'  '-' num2str(nr)];
    
    %plot1 = figure;
    plot([0:K-1], out(:,3), 'r');
    set(gcf, 'visible','off');
    title(['Róznica miedzy przewidywanym wyj?ciem a prawdziwym wyj?ciem (' dane ' - ' algorytm ')']);
    xlabel('nr danej przychodzacej');
    ylabel('Róznica');
    %legend('y przewidywane','y', 'RMSE');
    legend('y - y_p');
    
    MyBox = text(0.5,0.5,['algorytm=', algorytm]);
    set(MyBox,'Units','Normalized');
    set(MyBox,'Position',[0.02,0.9]);
    
    set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
    saveas(gcf, filename, 'png') %Save figure
  
    % S min max
    filename = ['G:\mgr\wykresy\' dane '\' dane '-' algorytm '-' TimeStamp '-s-min-max'  '-' num2str(nr)];
    
    %plot1 = figure;
    smi=cell2mat(S_min);
    sma= cell2mat(S_max);
    s=cell2mat(S);
    p=cell2mat(S_podmiana);
    n=cell2mat(S_nowy);
    plot([0:K-1], smi, '-*r',[0:K-1], sma, ':.g',[0:K-1], s, '--ob',[0:K-1], p,'+c',[0:K-1], n,'xm');
    set(gcf, 'visible','off');
    title([' (' dane ' - ' algorytm ')']);
    xlabel('nr danej przychodzacej');
    ylabel('S');
    %legend('y przewidywane','y', 'RMSE');
    legend('S min','S max','S');
    
    MyBox = text(0.5,0.5,['algorytm=', algorytm]);
    set(MyBox,'Units','Normalized');
    set(MyBox,'Position',[0.02,0.9]);
    
    set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
    saveas(gcf, filename, 'png') %Save figure
    
    
end