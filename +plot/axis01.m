function axis01(varargin)
% Plot an axis with grid on X=0 and Y=0 only.
% Package: plot
% Description: Create an axis that have only an the 0 X and Y axes.
% Input  : * Arbitrary pairs of arguments, ...,keyword,value,...
%            Options:
%            'MaxX'          - X range max, default is 1.   
%            'MaxY'          - Y range max, default is 1.   
%            'MinX'          - X range min, default is -1.   
%            'MinY'          - Y range min, default is -1.   
%            'NumberOfTickX' - Number of Ticks between 0 and MaxX,
%                              default is 10.
%            'NumberOfTickY' - Number of Ticks between 0 and MaxY,
%                              default is 10.
%            'TickSizeFrac'  - Ticks Size fraction as measured relative 
%                              to MaxX/MaxY, default is 0.025.
% Output : null
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: aplot.xis01; hold on; plot(rand(100,1),rand(100,1),'.');
%          plot.axis01('MaxX',2.5,'NumberOfTickY',5); hold on;
% Reliable: 2
%--------------------------------------------------------------------------

MaxX          = 1;       % max x-axis range
MaxY          = 1;       % max y-axis range
MinX          = -1;      % min x-axis range
MinY          = -1;      % min y-axis range
NumberOfTickX = 10;      % number of ticks in x-axis
NumberOfTickY = 10;      % number of ticks in y-axis
TickSizeFrac  = 0.025;   % size of ticks in range fraction

Nvar = length(varargin);
for I=1:2:Nvar-1
   switch varargin{I}
    case 'MaxX'
       MaxX          = varargin{I+1};
    case 'MaxY'
       MaxY          = varargin{I+1};    
    case 'MinX'
       MinX          = varargin{I+1};    
    case 'MinY'
       MinY          = varargin{I+1};    
    case 'NumberOfTickX'
       NumberOfTickX = varargin{I+1};
    case 'NumberOfTickY'
       NumberOfTickY = varargin{I+1};    
    case 'TickSizeFrac'
       TickSizeFrac  = varargin{I+1};
    otherwise
       error('Unknown keyword option');
   end
end

plot(0,0,'k.');
Haxis = gca;
NextPlot = get(Haxis,'NextPlot');
hold on;
axis([MinX MaxX MinY MaxY]);



Color = get(Haxis,'Color');
set(Haxis,'XColor',Color);
set(Haxis,'XTickLabel',[]);
set(Haxis,'YColor',Color);
set(Haxis,'YTickLabel',[]);

%--- plot central axis ---
plot([MinX;MaxX],[0;0],'k-');
plot([0;0],[MinY;MaxY],'k-');


%--- Plot ticks in x-axis ---
SemiLengthTickY = TickSizeFrac.*MaxY;
StepX = MaxX./NumberOfTickX;
for TickX=MinX:StepX:MaxX
   plot([TickX;TickX],[-SemiLengthTickY;SemiLengthTickY],'k-');
end
%--- Plot ticks in y-axis ---
SemiLengthTickX = TickSizeFrac.*MaxX;
StepY = MaxY./NumberOfTickY;
for TickY=MinY:StepY:MaxY
   plot([-SemiLengthTickX;SemiLengthTickX],[TickY;TickY],'k-');
end

%--- Plot Tick labels ---
Offset = 0.1;
text(MinX-Offset,0,num2str(MinX));
text(MaxX+Offset,0,num2str(MaxX));
text(0,MinY-Offset,num2str(MinY));
text(0,MaxY+Offset,num2str(MaxY));


%--- set NextPlot default ---
set(Haxis,'NextPlot',NextPlot);
