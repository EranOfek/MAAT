%This cript plots output of DarkEvolution (Median and STD against elapsed
%time. This does NOT return in # of electrons.


for CH = 1:4                    %plot medianMean and medianSTD against elapsedTimeStack, channel-by-channel
    subplot(4,2,2*CH-1)
    plot(elapsedTimeStack,medianMean(:,CH))
    title(strcat('CH',num2str(CH),' Mean Values DC Test'))
    xlabel('Time (s)')
    ylabel('Median Val')
   % axis([0 2000 4 9])
    hold on
    grid on
    
    subplot(4,2,2*CH)
    plot(elapsedTimeStack,medianSTD(:,CH))
    title(strcat('CH',num2str(CH),' STD Values DC Test'))
     xlabel('Time (s)')
    ylabel('STD Val')
   % axis([0 2000 6 14])
    hold on
    grid on
end


