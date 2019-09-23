%% Calculating running correlations
%
% Data Matrix 2 contains: Column ...
% 1 - Law Dome CO2
% 2 - Araucaria (normalised)
% 3 - Austrocedrus (normalised)
% 4 - Halocarpus (normalised)
% 5 - Ice dust (normalised)
% 6 - Nothofagus (normalised)
% 7 - Callitris (normalised)
% 8 - EOF1 (note: EOF1 and 2 only span 1:397)
% 9 - EOF2

%correlation = zeros(5,6,398,5);   %first is window, second is index column
%correlation2 = zeros(5,9,398,8);
%pval = zeros(2,7,398,6);
clear
window = [11 21 51 71 101]; %window size
load 'data_for_correlations.mat'

for i = 1:5
    for j = 1:6
        indexColumn = j;
        windowSize = window(i); 
        [correlationTS correlationTS2 PvalTS] = movingCorrelation(dataMatrix, windowSize, indexColumn);
        correlation(i,j,1:398,1:5) = correlationTS;
        pval(i,j,1:398,1:5) = PvalTS;
        correlation2(i,j,1:398,1:5) = correlationTS2;
    end
end


%% Plotting results

% Plot running correlations for Law Dome CO2

for k = 1:6
    for l = 1:5
        figure(l)
        subplot(3,2,k)
        plot(NAustro(1:398,1),squeeze(correlation(l,1,1:398,k)))
        hold on
        plot(NAustro(1:398,1),squeeze(correlation2(l,1,1:398,k)),'color','red','marker','+')
        line([1750 1750],[-1 1],'LineStyle','--','color','black')
        line([1850 1850],[-1 1],'LineStyle','--','color','black')
        axis([1600 2000 -1 1])
        hold off
    end
end

for i=1:5
    figure(i)
    subplot(3,2,1)
    title(['Law Dome CO_2 and Araucaria, ', num2str(window(i)), 'yr window (red + indicates significance at 95%)'])
    subplot(3,2,2)
    title(['Law Dome CO_2 and Austrocedrus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,3)
    title(['Law Dome CO_2 and Halocarpus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,4)
    title(['Law Dome CO_2 and Ice dust, ', num2str(window(i)), 'yr window'])
    subplot(3,2,5)
    title(['Law Dome CO_2 and Nothofagus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,6)
    title(['Law Dome CO_2 and Callitris, ', num2str(window(i)), 'yr window'])
end

% Plot running correlations for Araucaria

for k = 1:6
    for l = 1:5
        figure(l)
        subplot(3,2,k)
        plot(NAustro(1:398,1),squeeze(correlation(l,2,1:398,k)))
        hold on
        plot(NAustro(1:398,1),squeeze(correlation2(l,2,1:398,k)),'color','red','marker','+')
        line([1750 1750],[-1 1],'LineStyle','--','color','black')
        line([1850 1850],[-1 1],'LineStyle','--','color','black')
        axis([1600 2000 -1 1])
        hold off
    end
end

for i=1:5
    figure(i)
    subplot(3,2,1)
    title(['Araucaria and Law Dome CO_2, ', num2str(window(i)), 'yr window (red + indicates significance at 95%)'])
    subplot(3,2,2)
    title(['Araucaria and Austrocedrus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,3)
    title(['Araucaria and Halocarpus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,4)
    title(['Araucaria and Ice dust, ', num2str(window(i)), 'yr window'])
    subplot(3,2,5)
    title(['Araucaria and Nothofagus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,6)
    title(['Araucaria and Callitris, ', num2str(window(i)), 'yr window'])
end


% Plot running correlations for Austrocedrus

for k = 1:6
    for l = 1:5
        figure(l)
        subplot(3,2,k)
        plot(NAustro(1:398,1),squeeze(correlation(l,3,1:398,k)))
        hold on
        plot(NAustro(1:398,1),squeeze(correlation2(l,3,1:398,k)),'color','red','marker','+')
        line([1750 1750],[-1 1],'LineStyle','--','color','black')
        line([1850 1850],[-1 1],'LineStyle','--','color','black')
        axis([1600 2000 -1 1])
        hold off
    end
end

for i=1:5
    figure(i)
    subplot(3,2,1)
    title(['Austrocedrus and Law Dome CO_2, ', num2str(window(i)), 'yr window (red + indicates significance at 95%)'])
    subplot(3,2,2)
    title(['Austrocedrus and Araucaria, ', num2str(window(i)), 'yr window'])
    subplot(3,2,3)
    title(['Austrocedrus and Halocarpus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,4)
    title(['Austrocedrus and Ice dust, ', num2str(window(i)), 'yr window'])
    subplot(3,2,5)
    title(['Austrocedrus and Nothofagus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,6)
    title(['Austrocedrus and Callitris, ', num2str(window(i)), 'yr window'])
end

% Plot running correlations for Halocarpus

for k = 1:6
    for l = 1:5
        figure(l)
        subplot(3,2,k)
        plot(NAustro(1:398,1),squeeze(correlation(l,4,1:398,k)))
        hold on
        plot(NAustro(1:398,1),squeeze(correlation2(l,4,1:398,k)),'color','red','marker','+')
        line([1750 1750],[-1 1],'LineStyle','--','color','black')
        line([1850 1850],[-1 1],'LineStyle','--','color','black')
        axis([1600 2000 -1 1])
        hold off
    end
end

for i=1:5
    figure(i)
    subplot(3,2,1)
    title(['Halocarpus and Law Dome CO_2, ', num2str(window(i)), 'yr window (red + indicates significance at 95%)'])
    subplot(3,2,2)
    title(['Halocarpus and Araucaria, ', num2str(window(i)), 'yr window'])
    subplot(3,2,3)
    title(['Halocarpus and Austrocedrus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,4)
    title(['Halocarpus and Ice dust, ', num2str(window(i)), 'yr window'])
    subplot(3,2,5)
    title(['Halocarpus and Nothofagus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,6)
    title(['Halocarpus and Callitris, ', num2str(window(i)), 'yr window'])
end

% Plot running correlations for Ice dust

for k = 1:6
    for l = 1:5
        figure(l)
        subplot(3,2,k)
        plot(NAustro(1:398,1),squeeze(correlation(l,5,1:398,k)))
        hold on
        plot(NAustro(1:398,1),squeeze(correlation2(l,5,1:398,k)),'color','red','marker','+')
        line([1750 1750],[-1 1],'LineStyle','--','color','black')
        line([1850 1850],[-1 1],'LineStyle','--','color','black')
        axis([1600 2000 -1 1])
        hold off
    end
end

for i=1:5
    figure(i)
    subplot(3,2,1)
    title(['Ice Dust and Law Dome CO_2, ', num2str(window(i)), 'yr window (red + indicates significance at 95%)'])
    subplot(3,2,2)
    title(['Ice Dust and Araucaria, ', num2str(window(i)), 'yr window'])
    subplot(3,2,3)
    title(['Ice Dust and Austrocedrus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,4)
    title(['Ice Dust and Halocarpus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,5)
    title(['Ice Dust and Nothofagus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,6)
    title(['Ice Dust and Callitris, ', num2str(window(i)), 'yr window'])
end

% Plot running correlations for Nothofagus

for k = 1:6
    for l = 1:5
        figure(l)
        subplot(3,2,k)
        plot(NAustro(1:398,1),squeeze(correlation(l,6,1:398,k)))
        hold on
        plot(NAustro(1:398,1),squeeze(correlation2(l,6,1:398,k)),'color','red','marker','+')
        line([1750 1750],[-1 1],'LineStyle','--','color','black')
        line([1850 1850],[-1 1],'LineStyle','--','color','black')
        axis([1600 2000 -1 1])
        hold off
    end
end

for i=1:5
    figure(i)
    subplot(3,2,1)
    title(['Nothofagus and Law Dome CO_2, ', num2str(window(i)), 'yr window (red + indicates significance at 95%)'])
    subplot(3,2,2)
    title(['Nothofagus and Araucaria, ', num2str(window(i)), 'yr window'])
    subplot(3,2,3)
    title(['Nothofagus and Austrocedrus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,4)
    title(['Nothofagus and Halocarpus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,5)
    title(['Nothofagus and Nothofagus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,6)
    title(['Nothofagus and Callitris, ', num2str(window(i)), 'yr window'])
end


% Plot running correlations for Lake Tay record

for k = 1:6
    for l = 1:5
        figure(l)
        subplot(3,2,k)
        plot(NTay(1:398,1),squeeze(correlation(l,7,1:398,k)))
        hold on
        plot(NTay(1:398,1),squeeze(correlation2(l,7,1:398,k)),'color','red','marker','+')
        line([1750 1750],[-1 1],'LineStyle','--','color','black')
        line([1850 1850],[-1 1],'LineStyle','--','color','black')
        axis([1600 2000 -1 1])
        hold off
    end
end

for i=1:5
    figure(i)
    subplot(3,2,1)
    title(['Callitris and Law Dome CO_2, ', num2str(window(i)), 'yr window (red + indicates significance at 95%)'])
    subplot(3,2,2)
    title(['Callitris and Araucaria, ', num2str(window(i)), 'yr window'])
    subplot(3,2,3)
    title(['Callitris and Austrocedrus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,4)
    title(['Callitris and Halocarpus, ', num2str(window(i)), 'yr window'])
    subplot(3,2,5)
    title(['Callitris and Ice dust, ', num2str(window(i)), 'yr window'])
    subplot(3,2,6)
    title(['Callitris and Nothofagus, ', num2str(window(i)), 'yr window'])
end

% Plot running correlations for EOF1

for k = 1:8
    for l = 1:5
        figure(l)
        subplot(4,2,k)
        plot(EOF1(1:397,1),squeeze(correlation(l,8,1:397,k)))
        hold on
        plot(EOF1(1:397,1),squeeze(correlation2(l,8,1:397,k)),'color','red','marker','+')
        line([1750 1750],[-1 1],'LineStyle','--','color','black')
        line([1850 1850],[-1 1],'LineStyle','--','color','black')
        axis([1600 2000 -1 1])
        hold off
    end
end

for i=1:5
    figure(i)
    subplot(4,2,1)
    title(['EOF1 and Law Dome CO_2, ', num2str(window(i)), 'yr window (red + indicates significance at 95%)'])
    subplot(4,2,2)
    title(['EOF1 and Araucaria, ', num2str(window(i)), 'yr window'])
    subplot(4,2,3)
    title(['EOF1 and Austrocedrus, ', num2str(window(i)), 'yr window'])
    subplot(4,2,4)
    title(['EOF1 and Halocarpus, ', num2str(window(i)), 'yr window'])
    subplot(4,2,5)
    title(['EOF1 and Ice dust, ', num2str(window(i)), 'yr window'])
    subplot(4,2,6)
    title(['EOF1 and Nothofagus, ', num2str(window(i)), 'yr window'])
    subplot(4,2,7)
    title(['EOF1 and Callitris, ', num2str(window(i)), 'yr window'])
    subplot(4,2,8)
    title(['EOF1 and EOF2, ', num2str(window(i)), 'yr window'])
end

% Plot running correlations for EOF2

for k = 1:8
    for l = 1:5
        figure(l)
        subplot(4,2,k)
        plot(EOF1(1:397,1),squeeze(correlation(l,9,1:397,k)))
        hold on
        plot(EOF1(1:397,1),squeeze(correlation2(l,9,1:397,k)),'color','red','marker','+')
        line([1750 1750],[-1 1],'LineStyle','--','color','black')
        line([1850 1850],[-1 1],'LineStyle','--','color','black')
        axis([1600 2000 -1 1])
        hold off
    end
end

for i=1:5
    figure(i)
    subplot(4,2,1)
    title(['EOF2 and Law Dome CO_2, ', num2str(window(i)), 'yr window (red + indicates significance at 95%)'])
    subplot(4,2,2)
    title(['EOF2 and Araucaria, ', num2str(window(i)), 'yr window'])
    subplot(4,2,3)
    title(['EOF2 and Austrocedrus, ', num2str(window(i)), 'yr window'])
    subplot(4,2,4)
    title(['EOF2 and Halocarpus, ', num2str(window(i)), 'yr window'])
    subplot(4,2,5)
    title(['EOF2 and Ice dust, ', num2str(window(i)), 'yr window'])
    subplot(4,2,6)
    title(['EOF2 and Nothofagus, ', num2str(window(i)), 'yr window'])
    subplot(4,2,7)
    title(['EOF2 and Callitris, ', num2str(window(i)), 'yr window'])
    subplot(4,2,8)
    title(['EOF2 and EOF1, ', num2str(window(i)), 'yr window'])
end
