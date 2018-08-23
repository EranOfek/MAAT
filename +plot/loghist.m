function [Edges,N,H]=loghist(X,LogMinX,LogMaxX,Nbin,Color,varargin)
% Plot an histogram in log space.
% Package: plot
% Description: Plot an histogram in log space.
% Input  : - Vector of values for which to plot the histogram.
%          - log10 of histogram starting value
%          - log10 of histogram ending value
%          - Number of bins.
%          - Color, default is 'w'.
%          * Additional patch properties (...,keyword,value,...)
% Output : - X Edges of bins.
%          - Number of events in each bin.
%          - Vector of handles to each patch. 
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=logspace(5,10,1000)'.*(randn(1000,1)+1).*100;
%          loghist(X,5,12,8,'r');
%          [Edges,N,H]=plot.loghist(10.^(rand(10000,1)),-5, 1,10,'w',...
%                              'FaceColor','r');
%-------------------------------------------------------------------
if (nargin==4)
   Color = 'w';
else
   % do nothing
end

Edges = logspace(LogMinX,LogMaxX,Nbin);
Nbin  = length(Edges) - 1;
N     = histc(X,Edges);
for I=1:1:Nbin
   H(I) = patch([Edges(I);Edges(I+1);Edges(I+1);Edges(I)],[0;0;N(I);N(I)],Color);
   if (length(varargin)>0)
      set(H(I),varargin{:});
   end
   hold on;
end
set(gca,'XScale','log');

