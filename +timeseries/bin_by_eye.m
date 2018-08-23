function B=bin_by_eye(Data)
% Bin data by define bins interactively
% Package: timeseries
% Description: Bin light curve by eye. The function plots the light
%              curve and let the user to define bins using the mouse.
%              The user mark (with the mouse) the beginning and end
%              points of each bin. The left and right limits of each
%              user defined bin are marked on the plot using cyan and
%              red dashed lines, respectively.
% Input  : - Data matrix, [Time, Value, Error].
% Output : - Matrix with the following columns:
%            [Bin_Center,
%             <Y>,
%             StD(Y)/sqrt(N),
%             <X>,
%             Weighted-<Y>,
%             Formal-Error<Y>,
%             N]
%            In case that no errors are given, columns 5 and 6
%            will be zeros.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2003
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: binning.m, runmean.m
% Reliable: 1
%------------------------------------------------------------------------------
import Util.stat.*


ColT = 1;
ColM = 2;
ColE = 3;

errorxy(Data);
hold on;
YLim = get(gca,'YLim');

R = 'y';

zoom on;
disp(sprintf('\n Zoom is on, use maouse to set zoom'));
R = input('  Strike y to continue, other keys to stop\n','s');

B       = zeros(0,7);
LoopInd = 0;
while (R=='y'),
   LoopInd = LoopInd + 1;
   disp(sprintf('\n'));
   disp(sprintf('\n mark beginning and end points for current bin'));
   [X,Y] = ginput(2);
   if (X(1)>X(2)),
      T1 = X(2);
      T2 = X(1);
   else
      T2 = X(2);
      T1 = X(1);
   end
   I       = find(Data(:,ColT)>=T1 & Data(:,ColT)<T2);
   SubData = Data(I,:); 

   [WM,WE] = wmean(SubData(:,[ColM, ColE]));
   N       = length(SubData(:,1));
   B(LoopInd,:) = [0.5.*(T1+T2),...
   	           mean(SubData(:,ColM)),...
		    std(SubData(:,ColM))./sqrt(N),...
		   mean(SubData(:,ColT)),...
		   WM,...
		   WE,...
		   N];   

   plot([T1,T1],YLim,'c--');
   plot([T2,T2],YLim,'r--');


   % next bin: y/n
   R = input('  Strike y to continue, other keys to stop\n','s');
end


