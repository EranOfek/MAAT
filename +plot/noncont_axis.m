function [Hax]=noncont_axis(AxisType,Range,RangeC,varargin)
%-------------------------------------------------------------------
% noncont_axis function                                    plotting
% Description: Create non continuous axis.
% Input  : - Axis type {'x' | 'y'}, default is 'x'. Set the
%            non continuous axis.
%          - Data range in each one of the noncontinuus axes
%            [Min Max; Min Max;...].
%          - Data range of axes of continuus axis [Min Max].
%          * Aribtrary number of pairs of input arguments,
%            ...,keyword,value,...
%            'XScale'    - {'linear','log'}, default is 'linear'.
%            'YScale'    - {'linear','log'}, default is 'linear'.
%            'AxesPosX'  - position of X-axes bounderies,
%                          default is [0.130 0.905].
%            'AxesPosY'  - position of Y-axes bounderies,
%                          default is [0.110 0.925].
%            'Space'     - Spaceing between axes, default is 0.03.            
% Output : - Vector of handles for axes.
% Tested : Matlab 7.0
%     By : Eran O. Ofek        October 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: Hax=noncont_axis('x',[10   70;  872 890],...
%                               [1 100],'YScale','log');
%          % To change current axis use:
%          set(gcf,'CurrentAxes',Hax(1))
%-------------------------------------------------------------------
BoxLineColor   = [0 0 0];
BoxLineStyle   = '-';
BoxNoLineColor = [1 1 1];

%DefAllAxesPos  = [0.13 0.11 0.775 0.815];    % not used

XScale         = 'linear';
YScale         = 'linear';
AxesPosX       = [0.130 0.905];
AxesPosY       = [0.110 0.925];
Space          = 0.03;

Narg           = length(varargin);

for Iarg=1:2:Narg-1,
   switch varargin{Iarg}
    case 'XScale'
       XScale     = varargin{Iarg+1};
    case 'YScale'
       YScale     = varargin{Iarg+1};
    case 'AxesPosX'
       AxesPosX   = varargin{Iarg+1};
    case 'AxesPosY'
       AxesPosY   = varargin{Iarg+1};
    case 'Space'
       Space      = varargin{Iarg+1};
    otherwise
       error('Unknown keyword Option');
   end
end

switch AxisType
 case 'x'
    DefPosC    = AxesPosY;
    DefPos     = AxesPosX;
 case 'y'
    DefPosC    = AxesPosX;
    DefPos     = AxesPosY;
 otherwise
    error('Unknown AxisType Option');
end

PosC = DefPosC;


Nax = size(Range,1);   % Number of axes

%--- Construct axes position based on Range ---
DiffRange  = Range(:,2)-Range(:,1);    % range of each "subplot"
% available space for plot
AvailSpace = DefPos(2)-DefPos(1) - Space.*(Nax-1);

AxesWidth  = AvailSpace.*DiffRange./sum(DiffRange);

Pos       = zeros(Nax,2);

for Iax=1:1:Nax,
   if (Iax==1),
      Pos(Iax,1)  = DefPos(1);
   else
      Pos(Iax,1)   = Pos(Iax-1,2) + Space;
   end
   Pos(Iax,2) = Pos(Iax,1) + AxesWidth(Iax);
end



switch lower(AxisType)
 case 'x'
    Position = [Pos(:,1), PosC(1).*ones(Nax,1), Pos(:,2)-Pos(:,1), (PosC(2)-PosC(1)).*ones(Nax,1)];
 case 'y'
    Position = [PosC(1).*ones(Nax,1), Pos(:,1), (PosC(2)-PosC(1)).*ones(Nax,1), Pos(:,2)-Pos(:,1)];
 otherwise
    error('Unknown AxisType Option');
end


Hax = zeros(Nax,1);
for Iax=1:1:Nax,
   H = axes('Position',Position(Iax,:));
   Hax(Iax) = H;
   box off;
   
   switch lower(AxisType)
    case 'x'

       % add x-axis upper line of box
       H=line(Range(Iax,:),[RangeC(2), RangeC(2)]);
       set(H,'Color',BoxLineColor,'LineStyle',BoxLineStyle);

       if (Iax>1),
          % remove y-axis labels
          set(Hax(Iax),'YTickLabel',[],'YTick',[]);
          % remove y-axis line of box
	  H=line([Range(Iax,1),Range(Iax,1)],RangeC);
          set(H,'Color',BoxNoLineColor,'LineStyle','-');

       end
       if (Iax==Nax),
          %  add y-axis right line of box
	  H=line([Range(Iax,2),Range(Iax,2)],RangeC);
          set(H,'Color',BoxLineColor,'LineStyle',BoxLineStyle);
       end

       %--- set axis ---
       set(Hax(Iax),'XLim',Range(Iax,:));
       set(Hax(Iax),'YLim',RangeC);
       set(Hax(Iax),'XScale',XScale);
       set(Hax(Iax),'YScale',YScale);

       set(Hax(Iax),'XLimMode','manual','YLimMode','manual');
       hold on; 

       if (Iax~=Nax),
	 Ht = text(Range(Iax,2),RangeC(1),'~');

	  set(Ht,'FontSize',16,...
                 'Rotation',90,...
                 'Units','normalized',...
                 'HorizontalAlignment','center')
      end


    case 'y'

       % add y-axis right line of box
       H=line([RangeC(2), RangeC(2)],Range(Iax,:));
       set(H,'Color',BoxLineColor,'LineStyle',BoxLineStyle);

       if (Iax>1),
          % remove x-axis labels
          set(Hax(Iax),'XTickLabel',[],'XTick',[]);
          % remove x-axis line of box
          H=line(RangeC,[Range(Iax,1),Range(Iax,1)]);
          set(H,'Color',BoxNoLineColor,'LineStyle','-');          
       end
       if (Iax==Nax),
          %  add x-axis upper line of box
	  H=line(RangeC,[Range(Iax,2), Range(Iax,2)]);
          set(H,'Color',BoxLineColor,'LineStyle',BoxLineStyle);
       end

       %--- set axis ---
       set(Hax(Iax),'YLim',Range(Iax,:));
       set(Hax(Iax),'XLim',RangeC);

       set(Hax(Iax),'XLimMode','manual','YLimMode','manual');
       hold on; 

    otherwise
       error('Unknown AxisType Option');
   end
end
