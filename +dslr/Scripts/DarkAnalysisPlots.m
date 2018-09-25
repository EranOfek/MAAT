%Plots results from DCNew (median, variance vs. exposure time)
%in units of e and e^2 (calculated from G vector) and fits a line to each.
%
%Note: REMEMBER TO ENTER: G = [gain1, gain2, gain3, gain4]

medianSTDStore = squeeze(medianSTDStore);
medianMeanStore = squeeze(medianMeanStore);

medianSTDStore = medianSTDStore./G;          %divide by gain to convert to electrons
medianMeanStore = medianMeanStore./G;

varStore = medianSTDStore.*medianSTDStore;
meanStore = medianMeanStore;
 
figure

for CH = 1:4
    subplot(4,2,(2*CH-1))
    p = polyfit(exptimes,meanStore(:,CH)',1);
    f = polyval(p,exptimes);
    plot(exptimes,f,exptimes,meanStore(:,CH))
    legend(strcat('m=',num2str(p(1)),', c=',num2str(p(2))))
    title(strcat('CH',num2str(CH),' Median Values'))
    xlabel('Exp Time (s)')
    ylabel('Median Pix. Val. (e)')
    axis([0 30 1.5 7])
    hold on
    grid on

    subplot(4,2,(2*CH))
    p = polyfit(exptimes,varStore(:,CH)',1);
    f = polyval(p,exptimes);
    plot(exptimes,f,exptimes,varStore(:,CH))
    title(strcat('CH',num2str(CH),' Variance'))
    legend(strcat('m=',num2str(p(1)),', c=',num2str(p(2))))
    xlabel('Exp Time (s)')
    ylabel('Pixel Variance (e^2)')
    %axis([0 30 7 14])
    hold on
    grid on
end
    
%intercept gives read noise? store intercepts?

clear CH
clear A
clear B