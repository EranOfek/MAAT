%This script outputs a scatter plot (variance vs mean) of the data acquired and processed in
%FlatCapsAnalysis


for m = 1:size(exptimes,2)
   for CH = 1:4
        subplot(4,1,CH)
        xdata = reshape(medianStore(:,:,CH,m),1,[]);
        ydata = reshape(varStore(:,:,CH,m),1,[]);
        %fit = polyfit(xdata,ydata,1);      %create linear fit
        %plot(polyval(fit,1:35000))         %plot fit
        hold on
        scatter(xdata,ydata,2)          %plot variance vs median with small point size 
        title(strcat('CH',num2str(CH),' Median vs Variance'))
        xlabel('Median Pixel Val')
        ylabel('Pixel Variance')
        %axis([0 30 5 9])
        hold on
        grid on
   end
end

hold off