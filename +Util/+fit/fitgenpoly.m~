function [Res]=fitgenpoly(X,Y,DelY,Deg,varargin)
%------------------------------------------------------------------------------
% fitgenpoly function                                                   FitFun
% Description: Fit general polynomials (e.g., normal, Legendre, etc) to data.
%              The fitted function has the form:
%              Y = P(N+1) + P(N)*f_1(X) + ... + P(1)*f_N(X).
%              Where f_i(X) are the polynomial functions.
%              For example,  in the case of normal polynomials
%              these are f_1(X) = X, f_2(X) = X^2, etc.
% Input  : - Column vector of the independent variable.
%          - Column Vector of the dependent variable.
%          - Vector of the std error in the dependent variable.
%            If only one value is given then, points
%            are taken to be with equal weight. and Std error
%            equal to the value given.
%            Default is 1.
%          - Degree of highest polynomial.
%            Alternatively, this can be a vector of the specific polynomial
%            degrees to fit. For example [0 3],
%            will fit Y=P(4) + P(1)*X^3, and set P(2)=P(3)=0.
%            Default is 1.
%          * Arbitrary number of pairs of input arguments, ...keyword,value,...
%            where the following keywords are available:
%            'PolyType' - polynomial type {'norm','legen'}. Default is 'norm'.
%            'Algo'     - Least square algorith {'chol','orth'}. Default is 'chol'.
%            'NormX'    - Normlize X value to be between -1 and 1 {'y' | 'n'}.
%                         Default is 'n'.
%            'MaxIter'  - Maximum number of sigma clipping iteration.
%                         0 Means runs once with no sigma clipping. Default is 0.
%            'Method'   - Sigma clipping method {'StD','StdP','Constant','MinMax','Perc'}
%                         Default is 'StD'. See clip_resid.m for details.
%            'Mean'     - Sigma clipping method to calculate the mean {'mean','median'}.
%                         Default is 'median'. See clip_resid.m for details.
%            'Clip'     - Two elements vector containing the lower and upper
%                         values for the sigma clipping.
%                         This is [Lower, Upper] number of sigmas (positive)
%                         below/above the mean. See clip_resid.m for details.
%            'StdZ'     - Add epsilon to std {'y'|'n'}, default is 'y'.
%                         This is useful in case std is zero.
%            'Plot'     - Plot data points and best fit model.
%                         The following options are available.
%                         'none'    - Do not plot (default).
%                         'plot'    - normal plot, X vs. Y.
%                         'res'     - residuals plot, X vs. Y residulas.
%                         'fitonly' - plot only the best fit line on
%                                     exsiting plot.
%            'Xplot'    - Vector of points for which to plot the model.
%                         If empty matrix, then use default.
%                         Default is [min(X):range(X)./(N-1):max(X)].';
%            'SymData'  - Data points symbol type, default is 'ko'.
%                         If empty matrix then will not plot data points.
%            'ModelLine'- Model symbol type, default is 'r-'.
%                         If empty matrix then will not plot model line.
%            'Sort'     - Sort data points by thier X value before fiiting.
%                         This is important for plotting purposes
%                         {'y'|'n'}. Default is 'y'.
% Output : - Results structure containing the following fields:
%            .Par       - Vector of best fit parameters [P(N+1),P(N),...,P(1)].
%            .ParErr    - Vector of error in best fit parameters.
%            .Chi2      - \chi^{2} of best fit.
%            .Dof       - Number of degrees of freedom.
%            .Npar      - Number of free parameters.
%            .Ndata     - Number of data points.
%            .Cov       - The covarience matrix.
%                         Note that this is the coveriance matrix of the
%                         fit, so if NormX='y' this can not be used as is.
%            .Resid     - Vector of residuals.
%            .RMS       - StD of non-clipped residuals.
%            .ErrCL     - 1,2,3 sigma percentile lower and upper confidence
%                         interval.
%            .FlagUse   - For each residual, a flag indicating if residual is
%                         was sigma clipped (0) or not (1).
%            .NormXpoly - Linear transformatin used to normzlize X,
%                         X = .NormXpoly(1)*X + .NormXpoly(2)
%            .invNormXpoly- Linear transformatin used to return X to its
%                         original state.
% Tested : Matlab R0211b.
%     By : Eran O. Ofek                    Nov 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=[1:1:100].'; Y = 0.1 + X - 0.01.*X.^3+randn(100,1).*0.01; 
%          [Res]=fitgenpoly(X,Y,0.01,3)
% Reliable: 2
%------------------------------------------------------------------------------

Def.DelY = 1;
Def.Deg  = 1;
if (nargin==2),
   DelY   = Def.DelY;
   Deg    = Def.Deg;
elseif (nargin==3),
   Deg    = Def.Deg;
elseif (nargin>3),
   % do nothing
else
   error('Illegal number of input arguments');
end

Nplot = 100;

DefV.PolyType   = 'norm';
DefV.Algo       = 'chol';
DefV.NormX      = 'y';
DefV.MaxIter    = 0;
DefV.Method     = 'StD';
DefV.Mean       = 'median';
DefV.Clip       = [];
DefV.StdZ       = 'y';
DefV.Plot       = 'none';
DefV.Xplot      = (min(X):range(X)./(Nplot-1):max(X)).';
DefV.SymData    = 'ko';
DefV.ModelLine  = 'r-';
DefV.Sort       = 'y';
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

if (isempty(DelY)),
    DelY = 1;
end
Nplot = length(InPar.Xplot);

if (size(X,1)==1),
   X = X.';
end
if (size(Y,1)==1),
   Y = Y.';
end

N_X  = length(X);
N_Y  = length(Y);
N_DY = length(DelY);
if (N_X~=N_Y),
   error('X and Y must have the same length');
end

if (N_X~=N_DY),
   if (N_DY==1)
      % take equal weights
      DelY = DelY.*ones(N_X,1);
   else
      error('Y and DelY must have the same length');
   end
end

switch lower(InPar.Sort)
 case 'y'
    [~,SI] = sort(X);
    X    = X(SI);
    Y    = Y(SI);
    DelY = DelY(SI);
 otherwise
    % do nothing
end

if (length(Deg)==1),
   Deg = (0:1:Deg);
end

% degree of freedom
Res.Ndata = N_X;
Res.Npar  = length(Deg);
Res.Dof   = Res.Ndata - Res.Npar;

Deg       = sort(Deg);
Deg       = fliplr(Deg);
MaxDeg    = max(Deg);

% normalize x to be between -1 and 1
OrigX     = X;
OrigXplot = InPar.Xplot;
switch lower(InPar.NormX)
 case 'y'
    MinX   = min(X);
    RangeX = range(X);

    Res.NormXpoly    = [2./RangeX, -MinX.*2./RangeX - 1];
    Res.invNormXpoly = [1./Res.NormXpoly(1), -Res.NormXpoly(2)./Res.NormXpoly(1)];
    X = (X - MinX)./RangeX .*2 - 1;

    %InPar.Xplot      = (InPar.Xplot - MinX)./RangeX .*2 - 1;
 otherwise
    Res.NormXpoly    = [1 0];
    Res.invNormXpoly = [1 0];
    % do nothing
end


% building the design H matrix
switch lower(InPar.PolyType)
 case 'norm'
    % [X^N... X^1 1]
    H     = ones(Res.Ndata,MaxDeg+1);
    Hplot = ones(Nplot,MaxDeg+1);
    for Ip=1:1:MaxDeg,
       H(:,Ip+1)     = X.^Ip;
       Hplot(:,Ip+1) = InPar.Xplot.^Ip;
    end 
    
    H     = H(:,Deg+1);
    %H     = fliplr(H);
    Hplot = Hplot(:,Deg+1);
    %Hplot = fliplr(Hplot);
 
 case 'legen'
    % [X^N... X^1 1]
    H     = legendre(MaxDeg,X)';
    Hplot = legendre(MaxDeg,InPar.Xplot)';
    % remove specific columns not in Deg
    H     = H(:,Deg+1);
    %H     = fliplr(H);
    Hplot = Hplot(:,Deg+1);
    %Hplot = fliplr(Hplot);

 otherwise
    error('Unknown PolyType option');
end


[Res.Par,Res.ParErr,~,Res.Cov] = lscov(H,Y,1./DelY.^2,InPar.Algo);
Res.Resid = Y - H*Res.Par;
Res.FlagUse = ones(size(Res.Resid));

Iter = 0;
while (Iter<InPar.MaxIter),
   Iter = Iter + 1;
   [Iclip,~,Res.FlagUse] = clip_resid(Res.Resid,varargin{:});
   
   [Res.Par,Res.ParErr,~,Res.Cov] = lscov(H(Iclip,:),Y(Iclip),1./DelY(Iclip).^2,InPar.Algo);
   Res.Resid = Y - H*Res.Par;
end

% UN-normzlize X scale
switch lower(InPar.NormX)
    case 'y'
        % FFU
    otherwise
        % do nothing
end

Res.RMS   = std(Res.Resid(Res.FlagUse==1));
Res.ErrCL = err_cl(Res.Resid(Res.FlagUse==1));
Res.Chi2  = sum((Res.Resid(Res.FlagUse==1)./DelY(Res.FlagUse==1)).^2);
% return X variable to original state

Res.Par   = polysubstitution(Res.Par,Res.NormXpoly);

% plot data and model 
switch lower(InPar.Plot)
 case {'none','no','n'}
    % do nothing
 case 'plot'
    % plot the data points and the fit
    if (~isempty(InPar.SymData))
       plot(X,Y,'ko','MarkerFaceColor','k');
       hold on;
    end
    if (~isempty(InPar.ModelLine))
       %plot(X,H*Res.Par,'r-');
       plot(InPar.Xplot,Hplot*Res.Par,'r-');
    end
   
 case 'fitonly'
     % plot fitted line on exsiting plot
     % assume hold on is active
     %plot(X,H*Res.Par,'r-');
     plot(InPar.Xplot,Hplot*Res.Par,'r-');
 case 'res'
    if (~isempty(InPar.SymData))
       plot(X,Res.Resid,'ko','MarkerFaceColor','k');
       hold on;
    end
    if (~isempty(InPar.ModelLine))
       %plot(X,zeros(size(X)),'r-');
       plot(InPar.Xplot,Hplot*Res.Par,'r-');
    end

 otherwise
    error('Unknown Plot option');
end
