function H_a=splabel(LabelType,Dim,FontSize);
%--------------------------------------------------------------------
% splabel function                                          plotting
% Description: Draw a spectral type labeles on an axis of a plot.
%              Assuming the current axis values are B-V color index.
% Input  : - Label type:
%            's' : spectral type labels. (default).
%            't' : effective temperature labels.
%          - Labels dimension.
%            1 : X axis. (default).
%            2 : Y axis.
%            3 : Z axis.
%          - Font size, (default is 12).
% Output : - axis handle.
% Tested : Matlab 5.3
%     By : Eran O. Ofek         November 1999
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html 
%--------------------------------------------------------------------
if (nargin==0),
   LabelType = 's';
   Dim       = 1;
   FontSize  = 12;
elseif (nargin==1),
   Dim       = 1;
   FontSize  = 12;
elseif (nargin==2),
   FontSize  = 12;
elseif (nargin==3),
   % no default
else
   error('Illigal number of arguments');
end


ColorIndex   = [-0.4;-0.3; +0.0; +0.27; +0.58; +0.89; +1.45]; % B-V color index.
SpectralType = ['O0';'B0';'A0';'F0';'G0';'K0';'M0'];
EffTemp      = ['50000'; ' 9900'; ' 7400'; ' 6000'; ' 4900'; ' 3500'];

H_a = gca;
if (Dim==1),
   Lim       = 'XLim';
   Tick      = 'XTick';
   TickLabel = 'XTickLabel';
elseif (Dim==2),
   Lim       = 'YLim';
   Tick      = 'YTick';
   TickLabel = 'YTickLabel';
elseif (Dim==3),
   Lim       = 'ZLim';
   Tick      = 'ZTick';
   TickLabel = 'ZTickLabel';
else
   error('Illigal dimension');
end


set(H_a,Tick,ColorIndex);
if (LabelType=='t'),
   % plot temperature
   set(H_a,TickLabel,EffTemp);
elseif (LabelType=='s'),
   % plot spectral type
   set(H_a,TickLabel,SpectralType);
else
   error('Illigal Label Type');
end
set (H_a,'FontSize',FontSize);
