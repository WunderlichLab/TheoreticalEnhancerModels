
function figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize)
        figure(figStart)
        figStart = figStart + 1;

        tempFigFormat = 'svg';
        ax = gca;
        ax.FontSize = ticksFontSize; 
        xlabel(strcat('\textbf{',xLabel,'}'),'fontweight', ...
        'bold','fontsize',fs,'Interpreter','latex') 
        ylabel(strcat('\textbf{',yLabel,'}'),'fontweight', ...
            'bold','fontsize', fs,'Interpreter','latex')
        set(gcf, 'Position', get(0, 'Screensize'));
        xticks(0:5 + 1)
        % lgd = legend('T_2 sites = 0','T_2 sites = 1','T_2 sites = 2','T_2 sites = 3');
        % lgd.FontSize = lgdFontSize; 
        saveas(gcf,char(tempFigName), tempFigFormat)
        % to see if this fixes a bug where figures are being plotted too quickly
        pause(0.4)
        close gcf