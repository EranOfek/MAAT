function [P,Res,IndUse,IndNotUse,H]=polyfit_sc(X,Y,N,SigClip,MaxIter,PlotOption);
%------------------------------------------------------------------------------
% polyfit_sc function                                                   FitFun
% Description: Fitting a polynomial, of the form
%              Y=P(1)*X^N + ... +  P(N)*X + P(N+1)
%              to data with no errors, but with sigma clipping: 
% Input  : - Vector of X data points.
%          - Vector of Y data points.
%          - Degree of polynomials, if empty then interactively ask
%            the user for polynomial degree.
%          - Sigma clipping [Low High], default is [-Inf Inf], in units
%            of residuals StD.
%          - Maximum number of sigma clipping iterations, default is 1.
%            (0 means no-sigma clipping).
%          - Plot options:
%            [0]  - Do not plot anything (default).
%            [1]  - Plot best fit line in current axis.
%            [2]  - Plot best fit line in current axis, and plot the
%                   residuals in additional figure.
%                   If Handle parameter is given, then remove
%                   the graphic handle from plot before plotting
%                   the new fit.
% Output : - Best fit polynomial coefficients [P(1),..., P(N+1)].
%          - Vector of best fit residuals.
%          - Indices of data points used in the final iteration.
%          - Indices of data points not used in the final iteration.
%          - Handle for plotted line
% See also: fitpoly.m, fit_lin.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                           April 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: X=[1:1:100]'; Y=0.1+0.02.*X-0.01.*X.^2+rand(100,1).*0.1;
%          [Par,Res] = polyfit_sc(X,Y,2,[-1 1]);
% Reliable: 2
%------------------------------------------------------------------------------
Marker        = 'r-';

DefSigClip    = [-Inf Inf];
DefMaxIter    = 1;
DefPlotOption = [0];
if (nargin==3),
   SigClip    = DefSigClip;
   MaxIter    = DefMaxIter;
   PlotOption = DefPlotOption;
elseif (nargin==4),
   MaxIter    = DefMaxIter;
   PlotOption = DefPlotOption;
elseif (nargin==5),
   PlotOption = DefPlotOption;
elseif (nargin==6),
   % do nothing
else
   error('Illegal number of input arguments');
end


if (isempty(N)==1),
  N = input('Enter degree of polynomial to fit : ');
end

SigClip = abs(SigClip);

P       = polyfit(X,Y,N);
FitY    = polyval(P,X);
Res     = Y - FitY;
StD     = std(Res);
IndUse  = find(Res>=(-StD.*SigClip(1)) & Res<=(StD.*SigClip(2)));
IndNotUse = find(Res<(-StD.*SigClip(1)) | Res>(StD.*SigClip(2)));
if (length(IndNotUse==0)),
   Cont      = 0;
else
   Cont      = 1;
end


Niter = 0;
while (Niter<MaxIter & Cont==1),
   Niter = Niter + 1;
   P = polyfit(X(IndUse),Y(IndUse),N);

   FitY    = polyval(P,X);
   Res     = Y - FitY;
   StD     = std(Res);
   IndUse  = find(Res>=(-StD.*SigClip(1)) & Res<=(StD.*SigClip(2)));
   IndNotUse = find(Res<(-StD.*SigClip(1)) | Res>(StD.*SigClip(2)));
   if (length(IndNotUse==0)),
      Cont      = 0;
   else
      Cont      = 1;
   end
end

switch PlotOption
 case 0
    % do nothing
    H = [];
 case 1
    [SX,IS] = sort(X);
    H = plot(SX,FitY(IS),Marker);
 case 2
    [SX,IS] = sort(X);
    H = plot(SX,FitY(IS),Marker);
    
    
    Hres = figure(200);
    clf;
    Hres = plot(SX,Res(IS),Marker);

    set(0,'CurrentFigure',H)
 otherwise
    error('Unknwon PlotOption option');
end
