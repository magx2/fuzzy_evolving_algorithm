function [ RMSE ] = wyniki( y_przewidywane, y_prawdziwe, R_w_czasie, opis, nr, S, S_min, S_max, S_podmiana, S_nowy, tismp )
% opis = [ nazwa_danych algorytm R r OMEGA z_gwiazdka x ]

    warning('off','MATLAB:Axes:NegativeDataInLogAxis')

    
    
    dane = opis{1};
    katalog = opis{100};
    algorytm = opis{2};
    R = opis{3};
    r = opis{4};
    OMEGA = opis{5};

    
    
%     tismp = TimeStamp;
    wynikiFolder = [  tismp  '\' num2str(nr) '\' ];
%     mkdir(['G:\mgr\wykresy\' dane '\']);
    mkdir(['G:\mgr\wykresy\' katalog '\' wynikiFolder]);
    
    
    
    K = length(y_prawdziwe);

    out = [ y_przewidywane, y_prawdziwe, (y_przewidywane - y_prawdziwe)];
    
    RMSE = sqrt(mean(out(:,3).^2));
    
    % normlny
        filename = ['G:\mgr\wykresy\' katalog '\' wynikiFolder katalog '-' algorytm '-' TimeStamp '-' num2str(nr)];

        plot([0:K-1], out(:,1), 'g',[0:K-1], out(:,2), '--b');
        set(gcf, 'visible','off');
        title('{Wykres warto$\acute{s}$ci przewidywanej, prawdziwej}','interpreter','latex');
        xlabel('Nr danej');
        ylabel('{Warto$\acute{s}\acute{c}$ wyj$\acute{s}$cia}', 'Interpreter', 'latex');
        h = legend('${\hat y}$','$y$');
        set(h, 'Interpreter', 'latex')

        MyBox = text(0.5,0.5,[ 'dane: ' dane char(10) 'algorytm: ', algorytm, char(10), 'RMSE=', num2str(RMSE) char(10) 'R=' num2str(R) char(10) 'r=' num2str(r) char(10) '\Omega=' num2str(OMEGA) ]);
        set(MyBox,'Units','Normalized');
        set(MyBox,'Position',[0.02,0.9]);

        set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
        saveas(gcf, filename, 'png') %Save figure
	
     % normlny i R w czasie
        filename = ['G:\mgr\wykresy\' katalog '\' wynikiFolder katalog '-' algorytm '-' TimeStamp '-R-i-y' num2str(nr)];

        
        min1 = min(out(:,1));
        max1 = max(out(:,1));
        
        min2 = min(out(:,2));
        max2 = max(out(:,2));
        
        min3 = min([min1, min2]);
        max3 = max([max1, max2]);
        
        ax = plotyy([[0:K-1]',[0:K-1]'], [out(:,1), out(:,2)],[0:K-1], R_w_czasie );
        set(gcf, 'visible','off');
        % title('Liczba klastrów i warto?ci wyj?ciowe');
        title('{Liczba klastr$\acute{o}$w $\hat{y}$ i $y$}','interpreter','latex')
        xlabel('Nr danej');
        ylabel(ax(1), '{Warto$\acute{s}\acute{c}$ wyj$\acute{s}$cia}', 'Interpreter', 'latex');
        ylabel(ax(2), '{Ilo$\acute{s}\acute{c}$ klastr$\acute{o}$w}', 'Interpreter', 'latex');
        h=legend('$y$','${\hat y}$','$R$');
        set(h, 'Interpreter', 'latex')
%         set(ax(1),'XScale','linear')
%         set(ax(1),'YScale','linear')
         set(ax(1),'YLim',[min3, max3])
        
        
        MyBox = text(0.5,0.5,[ 'dane: ' dane char(10) 'algorytm: ', algorytm, char(10), 'RMSE=', num2str(RMSE) char(10) 'R=' num2str(R) char(10) 'r=' num2str(r) char(10) '\Omega=' num2str(OMEGA) ]);
        set(MyBox,'Units','Normalized');
        set(MyBox,'Position',[0.02,0.9]);

        set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
        saveas(gcf, filename, 'png') %Save figure
    
        
    % roznica
        filename = ['G:\mgr\wykresy\' katalog '\' wynikiFolder katalog '-' algorytm '-' TimeStamp '-roznica'  '-' num2str(nr)];

        plot([0:K-1], out(:,3), 'r');
        set(gcf, 'visible','off');
%         title('{Róznica miedzy $\hat{y}$ i $y$}', 'interpreter', 'latex')
        title('{R$\acute{o}\dot{z}$nica miedzy $\hat{y}$ i $y$}','interpreter','latex')
        xlabel('Nr danej');
        ylabel('Róznica');
        h=legend('$y - \hat{y}$');
        set(h, 'Interpreter', 'latex')

        MyBox = text(0.5,0.5,[ 'dane: ' dane char(10) 'algorytm: ', algorytm, char(10), 'RMSE=', num2str(RMSE) char(10) 'R=' num2str(R) char(10) 'r=' num2str(r) char(10) '\Omega=' num2str(OMEGA) ]);
        set(MyBox,'Units','Normalized');
        set(MyBox,'Position',[0.02,0.9]);

        set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
        saveas(gcf, filename, 'png') %Save figure
    
    % R w czasie
%         filename = ['G:\mgr\wykresy\' katalog '\' wynikiFolder katalog '-' algorytm '-' TimeStamp 'R-w-czasie' '-' num2str(nr)];
% 
%         plot([0:K-1], R_w_czasie, 'g');
%         set(gcf, 'visible','off');
%         title(['Liczba klastórw w czasie (' dane ' - ' algorytm ')']);
%         xlabel('nr danej przychodzacej');
%         ylabel('Liczba klastrów');
%         h=legend('$R$');
%         set(h, 'Interpreter', 'latex')
% 
%         MyBox = text(0.5,0.5,['algorytm=', algorytm, char(10), 'RMSE=', num2str(RMSE) char(10) 'R=' num2str(R) char(10) 'r=' num2str(r) char(10) '\Omega=' num2str(OMEGA)]);
%         set(MyBox,'Units','Normalized');
%         set(MyBox,'Position',[0.02,0.9]);
% 
%         set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
%         saveas(gcf, filename, 'png') %Save figure
%     
	% R w czasie - delta y
%         filename = ['G:\mgr\wykresy\' katalog '\' wynikiFolder katalog '-' algorytm '-' TimeStamp 'R-w-czasie-delta' '-' num2str(nr)];
% 
%         plotyy([0:K-1], R_w_czasie,[0:K-1], out(:,3),'semilogy','plot');
%         set(gcf, 'visible','off');
%         title(['Liczba klastórw w czasie (' dane ' - ' algorytm ')']);
%         xlabel('nr danej przychodzacej');
%         ylabel('Liczba klastrów');
%         legend('R');
% 
%         MyBox = text(0.5,0.5,['algorytm=', algorytm]);
%         set(MyBox,'Units','Normalized');
%         set(MyBox,'Position',[0.02,0.9]);
% 
%         set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
%         saveas(gcf, filename, 'png') %Save figure
    
  
    % S min max
%         filename = ['G:\mgr\wykresy\' katalog '\' wynikiFolder katalog '-' algorytm '-' TimeStamp '-s-min-max'  '-' num2str(nr)];
% 
%         smi=cell2mat(S_min);
%         sma= cell2mat(S_max);
%         s=cell2mat(S);
%         p=cell2mat(S_podmiana);
%         n=cell2mat(S_nowy);
%         plot([0:K-1], smi, '-*r',[0:K-1], sma, ':.g',[0:K-1], s, '--ob',[0:K-1], p,'+c',[0:K-1], n,'xm');
%         set(gcf, 'visible','off');
%         title([' (' dane ' - ' algorytm ')']);
%         xlabel('nr danej przychodzacej');
%         ylabel('S');
%         legend('S min','S max','S');
% 
%         MyBox = text(0.5,0.5,['algorytm=', algorytm, char(10), 'RMSE=', num2str(RMSE) char(10) 'R=' num2str(R) char(10) 'r=' num2str(r) char(10) '\Omega=' num2str(OMEGA)]);
%         set(MyBox,'Units','Normalized');
%         set(MyBox,'Position',[0.02,0.9]);
% 
%         set(gcf, 'PaperSize', [5,5]); %Set the paper to have width 5 and height 5.
%         saveas(gcf, filename, 'png') %Save figure
    
    
end