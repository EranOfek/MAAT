function [H]=hist_ofplot(varargin)
% calculate histograms of the X and Y axis of the data in current figure.
% Package: plot
% Description: Given an existing plot, calculate histograms of the X and Y
%              axis of the data in current figure axis. The x-axis histogram
%              is displayed below the current plot, and the y-axis histogram
%              is dispalyed on the right side.
% Input  : * Arbitrary number of pairs of ...,keyword,value,...
%            Avaliable keywordsare:
%            'NbinX'     - Number of histograms bin in X axis.
%                          Default is 30.
%            'NbinY'     - Number of histograms bin in X axis.
%                          Default is 30.
%            'HistHight' - Fraction of figure occupying histogram hight.
%                          Default is 0.15.
%            'Gap'       - Gap between histograms and plot, default is 0.01.
%            'Scale'     - Scale of histogram N-axis.
%                          'linear' - linear scale (default).
%                          'log'    - log scale.
%                          'logdata'- plot log of N.
%            'Norm'      - Histogram normalization:
%                          'no'  - no normalization (default).
%                          'sum' - by sum.
%                          'max' - by max.
%            'Color'     - A flag indicating if to use the color of each
%                          children object or to use the same color
%                          for all the childrens.
%                          If empty matrix then will use the color of each
%                          children (default). Otherwise this is the color
%                          to use for all the histograms.
%            'Plot'      - Plot hisogram (true) or only plot axes (false).
%                          Default is true.
% Output : - Three element vectors containing handles to:
%            [main plot, X-axis hist, Y-axis hist].
% Output : - Handle to axes [main plot, x histogram, y histogram]
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Jul 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: plot(randn(10000,1).*0.2+round(rand(10000,1)),randn(10000,1).*0.5,'k.','MarkerSize',4);
%          hold on
%          plot(rand(1000,1).*0.3,randn(1000,1).*0.5,'r.','MarkerSize',4);
%          [H]=hist_ofplot('Norm','sum');
%          set(get(H(2),'Ylabel'),'String','Frac.','FontSize',16);
%          set(get(H(3),'Xlabel'),'String','Frac.','FontSize',16);  
% Web example: http://adsabs.harvard.edu/abs/2012PASP..124...62O
% Reliable: 1
%--------------------------------------------------------------------------

DefV.NbinX = 30;
DefV.NbinY = 30;
DefV.HistHight = 0.15;
DefV.Gap       = 0.01;
DefV.Scale     = 'linear'; % {'linear'|'log'|'logdata'}
DefV.Norm      = 'no';
DefV.Color     = [];
DefV.Plot      = true;
%InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

AxesPosOrig = get(gca,'Position');
Hgca = gca;
XLim = get(Hgca,'XLim');
YLim = get(Hgca,'YLim');
set(Hgca,'YLimMode','Manual');
set(Hgca,'XLimMode','Manual');

AxesPos    = AxesPosOrig;
AxesPos(3) = AxesPosOrig(3)-InPar.HistHight-InPar.Gap;
AxesPos(4) = AxesPosOrig(4)-InPar.HistHight-InPar.Gap;
set(gca,'Position',AxesPos);


%--- calculate histograms for current plot ---
Hchild = get(Hgca,'Children');
Nchild = length(Hchild);
if (~InPar.Plot)
    Nchild = 0;
end
for Ichild=1:1:Nchild
   [X,N] = Util.stat.realhist(get(Hchild(Ichild),'XData'),[XLim InPar.NbinX]);
   HistX(Ichild).X = X;
   HistX(Ichild).N = N;
   [X,N] = Util.stat.realhist(get(Hchild(Ichild),'YData'),[YLim InPar.NbinX]);
   HistY(Ichild).X = X;
   HistY(Ichild).N = N;
   if (isempty(InPar.Color))
      Color(Ichild).Color = get(Hchild(Ichild),'Color');
   else
      Color(Ichild).Color = InPar.Color;
   end
   switch lower(InPar.Scale)
    case 'logdata'
       HistX(Ichild).N = log10(HistX(Ichild).N);
       HistY(Ichild).N = log10(HistY(Ichild).N);
    otherwise
       % do nothing
   end
   switch lower(InPar.Norm)
    case 'no'
       % do nothing
    case 'max'
       HistX(Ichild).N = HistX(Ichild).N./max(HistX(Ichild).N);
       HistY(Ichild).N = HistY(Ichild).N./max(HistY(Ichild).N);
    case 'sum'
       HistX(Ichild).N = HistX(Ichild).N./sum(HistX(Ichild).N);
       HistY(Ichild).N = HistY(Ichild).N./sum(HistY(Ichild).N);
    otherwise
       error('Unknown Norm option');
   end
end


%--- plot X histhogram ---
HistXpos = [AxesPos(1), AxesPos(2)+AxesPos(4)+InPar.Gap, AxesPos(3), InPar.HistHight];
Hx = axes('position',HistXpos);


for Ichild=1:1:Nchild
   H=plot.hist_stairs(HistX(Ichild).X,HistX(Ichild).N,'k-','Type','v');
   set(H,'Color',Color(Ichild).Color);
   switch lower(InPar.Scale)
    case 'log'
       set(Hx,'YScale','log');
    otherwise
       % do nothing
   end
   hold on;
end
set(gca,'XLim',XLim);
set(gca,'XTickLabel',[]);

%--- plot Y histogram ---
HistYpos = [AxesPos(1)+AxesPos(3)+InPar.Gap, AxesPos(2), InPar.HistHight, AxesPos(4)];
Hy = axes('position',HistYpos);


for Ichild=1:1:Nchild
   H=plot.hist_stairs(HistY(Ichild).X,HistY(Ichild).N,'k-','Type','h');
   set(H,'Color',Color(Ichild).Color);
   switch lower(InPar.Scale)
    case 'log'
       set(Hy,'XScale','log');
    otherwise
       % do nothing
   end
   hold on;
end
set(gca,'YLim',YLim);
set(gca,'YTickLabel',[]);


H = [Hgca,Hx,Hy];
