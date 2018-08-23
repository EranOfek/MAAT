function [H2,H1,H0]=multi_axis(AxisType,Function,varargin)
% Create additional x or y axis related by function existing axes.
% Package: plot
% Description: Create additional x or y axis which is related to
%              the exsiting axis by a given function.
% Input  : - Which axis to add: {'x' | 'y'}.
%            If 'x', then new axis is added on top.
%            If 'y', then new axis is added on right.
%          - Function (inline object or function name as string),
%            that relates the two axes.
%            (The function should be monotonic!).
%          * Arbitrary number of additional (optional) parameters
%            to pass to "Function".
% Output : - Handle for new axis.
%          - Handle for current (old) axis.
%          - Handle for the completed box without ticks around figure.
% Plot   : Add upper x-axis or left y-axis to the plot.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Sep 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: F=inline('1./x','x');   % monotonic function
%          plot(rand(10,2));
%          [H2,H1]=multi_axis('x',F);
%          F=inline('(x.*10).^2+5+a','x','a');
%          [H3,H4]=multi_axis('y',F,100);
%          xlabel(H1,'axis1')
%          xlabel(H2,'axis2')
%          ylabel(H3,'axis3')
%          ylabel(H4,'axis4')
% Reliable: 1
%--------------------------------------------------------------------------

H1 = gca;    % first axis handle
box(H1,'off');


%--- Create complete box without ticks around figure ---
H0 = axes('Position',get(H1,'Position'),...
          'OuterPosition',get(H1,'OuterPosition'),...
              'Box','on',...
              'XTick',[],...
              'YTick',[],...
              'Color','none');

switch AxisType
 case 'x'
    % create second axis
    H2 = axes('Position',get(H1,'Position'),...
              'OuterPosition',get(H1,'OuterPosition'),...
              'Box','off',...
              'XAxisLocation','top',...
              'YTick',[],...
              'Color','none');

    Lim1   = get(H1,'XLim');
    Scale1 = get(H1,'XScale');
    
    Lim2   = feval(Function,Lim1,varargin{:});
    if (Lim2(2)>Lim2(1))
       %Lim2 = Lim2;
       Dir  = 'normal';
    else 
       Lim2 = [Lim2(2), Lim2(1)];
       Dir  = 'reverse';
    end
    set(H2,'XLim',Lim2,'XDir',Dir,'XScale',Scale1);

 case 'y'
    % create second axis
    H2 = axes('Position',get(H1,'Position'),...
              'OuterPosition',get(H1,'OuterPosition'),...
              'Box','off',...
              'YAxisLocation','right',...
              'XTick',[],...
              'Color','none');

    Lim1   = get(H1,'YLim');
    Scale1 = get(H1,'YScale');

    Lim2   = feval(Function,Lim1,varargin{:});
    if (Lim2(2)>Lim2(1))
       %Lim2 = Lim2;
       Dir  = 'normal';
    else 
       Lim2 = [Lim2(2), Lim2(1)];
       Dir  = 'reverse';
    end

    set(H2,'YLim',Lim2,'YDir',Dir,'YScale',Scale1);


 otherwise
    error('Unknown AxisType Option');
end

%linkaxes([H2;H0]);
Pos=get(H0,'Position');
set(H1,'Position',Pos);
set(H2,'Position',Pos);


