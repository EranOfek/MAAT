function [Par,ParErr,Chi2,Dof,Resid,Cov,IndUse]=fit_lin(Fun,Y,ErrY,X,MapX,Method,varargin)
% Fit a general function which linear in its free parameters.
% Package: Util.fit
% Description: Fit a general function which linear in its free parameters.
% Input  : - Cell array of strings. Each cell containining a function
%            which multiply one of the coefficients.
%            Alternatively, this can be the design matrix.
%            Note that in this case 'Plot' is set to [0 0 0].
%          - Vector of Y (dependent variable).
%          - Vector of errors in Y.
%          - Vector of X (independent variable).
%            or matrix of independent variables if several of them
%            are needed (i.e. multi-dimensional fit; e.g. Y=f[X;Y;Z]).
%          - Cell vector in which the i-th element containing the column
%            in the X matrix which the i-th function is refering to.
%            if the i-th function has more than one dependent variable,
%            than the i-th cell should contain the column index (in the
%            X matrix) for each dependent variable,
%            by order of thie appearance.
%            If cell is empty, than function is a constant.
%            Default is assuming one-dim fit MapX{*}=[1];
%            if empty ([]), set to default.
%            In case the first parameter is a numeric matrix,
%            this parameter is not used.
%          - Solution method:
%            'Inv'   - simple inversion, default.
%            'SVD'   - Singular value decomposition.
%          * Arbitrary number of pairs of arguments:
%            ...,Keyword,Value,...
%            Where the possible keywords are:
%            - 'SigClip'        : [Lower, Upper] in sigma units,
%                                 degault is [Inf Inf].
%            - 'Nrej'           : [low, high] number of pixels
%                                 to reject in each iteration,
%                                 default is [0 0].
%                                 Note that if sum(Nrej)>0 it override
%                                 the 'SigClip' option.
%            - 'MaxIter'        : Maximum number of sigma clipping
%                                 iterations, default is 1.
%            - 'SingThresh'     : Singularity threshold for SVD,
%                                 default is 1e-3.
%            - 'Plot'           : Plot flags:
%                                 [DataPlot, ResidualsFlag, HistFlag],
%                                 where default is [0 0 0]
%                                 (do not plot anything).
%                                 If DataPlot==1, then plot Y vs. X,
%                                 best fit curve and data.
%                                 If ResidualsFlag==1, then plot
%                                 residuals vs. X.
%                                 If HistFlag, then plot residuals histogram.
%                                 This option is working only for one-dim
%                                 functions
%                                 (i.e. the number of columns of X is 1).
%            - 'LineType'       : Best fit curve line type
%                                 {':' | '.-' | '--' | '-'}, default is '-'.
%            - 'LineColor'      : Best fit curve line color, default is 'r'.
%            - 'Marker'         : Used data points marker type, default is 'o'.
%            - 'MarkerColor'    : Used data points (and hsitogram)
%                                 marker color, default is 'b'.
%            - 'OutMarker'      : Un-used (points removed in sigma-clipping)
%                                 points marker type, default is '^'.
%                                 If NaN, do not plot sigma-clipped points.
%            - 'OutMarkerColor' : Un-used (points removed in sigma-clipping)
%                                 points (and histogram) marker color,
%                                 default is 'g'.
%            - 'MarkerFace'     : Plot marker face {'y' | 'n'}, default is 'y'.
%            - 'MarkerSize'     : Marker size, default is 6.
%            - 'HistBinSigma'   : Histogram bin size, in sigma units,
%                                 default is 0.2.
% Output : - Best fit parameters.
%          - Errors in best fit parameters.
%          - Chi2.
%          - Dof.
%          - Fit residuals.
%            Resid(IndUse) contains only the residuals that participated
%            in the last iteartion.
%          - Covariance matrix.
%          - Indices of used elements in the last sigma clipping iteration.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                July 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: % Fit the function: Y=A + B*X^2 + C*sin(Y) + D*X*Y
%          % set the functions:
%          N = 1000;
%          X = [rand(N,1).*100, rand(N,1).*100];
%          Fun{1}='1'; Fun{2}='X.^2'; Fun{3}='log(Y)'; Fun{4}='X.*Y';
%          MapX{1}=[]; MapX{2}=1; MapX{3}=2; MapX{4}=[1 2];
%          ParIn = [15; 0.1; 10; 0.02];
%          ErrY = 0.8;
%          [Y,H] = lin_fun(Fun,ParIn,X,MapX); Y = Y + randn(N,1).*ErrY;
%          [Par,ParErr,Chi2,Dof,Resid,Cov,IndUse]=fit_lin(Fun,Y,ErrY,X,MapX);
%          % example with plot
%          N = 200;
%          X = [rand(N,1).*100];
%          Fun{1}='1'; Fun{2}='X.^2'; Fun{3}='log(X)'; Fun{4}='sqrt(X)';
%          ParIn = [15; 0.1; 10; 0.02];
%          ErrY = 30;
%          [Y,H] = lin_fun(Fun,ParIn,X,[]); Y = Y + randn(N,1).*ErrY;
%          [Par,ParErr,Chi2,Dof,Resid,Cov,IndUse]=fit_lin(Fun,Y,ErrY,X,[],'SVD','Plot',[1 1 1],'SigClip',[1.8 1.8]);
% Reliable: 2
%--------------------------------------------------------------------------
Npt  = 1000;  % Number of points and curve plotting

Ny   = length(Y);
Nfun = length(Fun);

DefMethod          = 'Inv';
DefSigClip         = [Inf Inf];
DefNrej            = [0 0];
DefMaxSigClipIter  = 1;
DefSingThresh      = 1e-3;
if (nargin==4)
   for I=1:1:Ny
      MapX{I} = [1];
   end
   Method         = DefMethod;
elseif (nargin==5)
   Method         = DefMethod;
else
   % do nothing
end

if (isempty(MapX)==1)
   for I=1:1:Ny
      MapX{I} = [1];
   end
end

Narg = length(varargin);
if (0.5.*Narg~=floor(0.5.*Narg))
   error('Illegal number of input arguments');
end
% set default values
SigClip           = DefSigClip;
Nrej              = DefNrej;
MaxSigClipIter    = DefMaxSigClipIter;
SingThresh        = DefSingThresh;
Plot              = [0 0 0];
LineType          = '-';
LineColor         = 'r';
Marker            = 'o';
MarkerColor       = 'b';
OutMarker         = '^';
OutMarkerColor    = 'g';
MarkerFace        = 'y';
MarkerSize        = 6;
HistBinSigma      = 0.2;

for Iarg=1:2:Narg-1
   switch varargin{Iarg}
    case 'SigClip'
       SigClip        = varargin{Iarg+1};
    case 'Nrej'
       Nrej           = varargin{Iarg+1};
    case 'MaxIter'
       MaxSigClipIter = varargin{Iarg+1};
    case 'SingThresh'
       SingThresh      = varargin{Iarg+1};
    case 'Plot'
       Plot            = varargin{Iarg+1};
    case 'LineType'
       LineType        = varargin{Iarg+1};
    case 'LineColor'
       LineColor       = varargin{Iarg+1};
    case 'Marker'
       Marker          = varargin{Iarg+1};
    case 'MarkerColor'
       MarkerColor     = varargin{Iarg+1};
    case 'OutMarker'
       OutMarker       = varargin{Iarg+1};
    case 'OutMarkerColor'
       OutMarkerColor  = varargin{Iarg+1};
    case 'MarkerFace'
       MarkerFace      = varargin{Iarg+1};
    case 'MarkerSize'
       MarkerSize      = varargin{Iarg+1};
    case 'HistBinSigma'
       HistBinSigma    = varargin{Iarg+1};
    otherwise
       error('Unknown keyword Option');
   end
end

if (length(SigClip)==1)
   SigClip = [SigClip, SigClip];
end
LowerSigClip = SigClip(1);
UpperSigClip = SigClip(2);


if (size(ErrY,1)==1)
   ErrY = ErrY.*ones(Ny,1);
end
%if (size(ErrY,1)==size(ErrY,2) & size(ErrY,1)>1),
%   CovOption = 1;
%else
%   CovOption = 0;
%end

% build design matrix
EvalPar    = ones(Nfun,1);
if (isnumeric(Fun)==1)
   H = Fun;   % Fun conatins the design matrix
   Plot = [0 0 0];
else
   [EvalY, H] = lin_fun(Fun,EvalPar,X,MapX);
end

% sigma cliping
UseFlag       = ones(Ny,1);        % 1 if used; 0 if removed by SigClip
IndUse        = find(UseFlag==1);  
SigClipStatus = 1;
Isc           = 1;
while (SigClipStatus==1 && Isc<=MaxSigClipIter)
   % current design matrix to use
   Hs    = H(IndUse,:);
   Ys    = Y(IndUse,:);
   ErrYs = ErrY(IndUse,:);

   switch Method
    case 'Inv'

       % Old - Memory consuming
       % Var  = diag(ErrYs.^2);
       % Cov     = inv(Hs'*inv(Var)*Hs);
       % Par     = Cov*Hs'*inv(Var)*Ys;

       InvVar  = diag(1./(ErrYs.^2));
       Tmp = Hs'*InvVar;
       Cov     = inv(Tmp*Hs);
       Par     = Cov*Tmp*Ys;
       clear Tmp;

       % speed consuming
       % Cov     = inv(Hs'*InvVar*Hs);
       % Par     = Cov*Hs'*InvVar*Ys;

       ParErr  = sqrt(diag(Cov));

       %'Number of degree of freedom :', Freedom
       Resid = Y - H*Par;
       Chi2  = sum((Resid(IndUse)./ErrYs).^2);
       Dof   = length(IndUse) - length(Par);

    case 'SVD'
       [U,S,V]   = svd(Hs,0);
       DS        = diag(S);
       Ising     = find(diag(S)<SingThresh);
       Inonsing  = find(diag(S)>=SingThresh);
       IW        = 1./DS;
       IW(Ising) = 0;
       SizeV     = size(V);
       SizeU     = size(U);
       %V*[diag(IW)]*(U.')
       Par       = V*[diag(IW)]*(U.'*Ys);
       Cov       =  ((ones(Nfun,1)*(IW.^2)').*V)*V';
%       Inonsing
       ParErr    = sqrt(diag(Cov));

       %'Number of degree of freedom :', Freedom
       Resid = Y - H*Par;
       Chi2  = sum((Resid(IndUse)./ErrYs).^2);
       Dof   = length(IndUse) - length(Par) + length(Ising);

    otherwise
       error('Unknown Method Option');
   end

   if (sum(Nrej)>0)
      % use Nrej rejection schme:
      [SortedResid,SortInd] = sort(Resid);
      CurFlag = zeros(size(UseFlag));
      CurFlag(SortInd(Nrej(1):end-Nrej(2))) = 1; 
   else
      % use SigClip rejection scheme:
      CurFlag    = ((Resid./ErrY)>=-LowerSigClip & (Resid./ErrY)<=UpperSigClip);
   end
   UseFlag    = (UseFlag & CurFlag);
   IndUse     = find(UseFlag==1);
   IndNotUse  = find(UseFlag==0);

   Isc = Isc + 1;
end



switch MarkerFace
 case 'y'
    FaceColor    = MarkerColor;
    OutFaceColor = OutMarkerColor;
 case 'n'
    % do nothing
 otherwise
    error('Unknown MarkerFace Option');
end


%------------------------------------
%--- Plot Data and best fit curve ---
%------------------------------------
switch Plot(1)
 case 0
    % do nothing
 case 1
    % plot best fit curves and data
    
    if (size(X,2)>1)
       warning('Can not plot best fit curve for multi-dimensional fit');
    else
       Ix = 1;
       figure;
       % plot data
       errorxy([X(IndUse,Ix),   Y(IndUse),   ErrY(IndUse)],...
               'Marker',Marker,'EdgeColor',MarkerColor,...
   	       'FaceColor',FaceColor,'MarkSize',MarkerSize);
       hold on;
       errorxy([X(IndNotUse,Ix),Y(IndNotUse),ErrY(IndNotUse)],...
               'Marker',OutMarker,'EdgeColor',OutMarkerColor,...
   	       'FaceColor',OutFaceColor,'MarkSize',MarkerSize);

       % plot curve
       MinX = min(X(:,Ix));
       MaxX = max(X(:,Ix));
       Xpt  = [MinX:(MaxX-MinX)./Npt:MaxX].';
       
       Ypt  = lin_fun(Fun,Par,Xpt,MapX);       
       plot(Xpt,Ypt,LineType,'Color',LineColor);

       % set fonts and labels
       set(gca,'FontSize',14);
       H=xlabel(sprintf('X_{%d}',Ix));
       set(H,'FontSize',18);
       H=ylabel('Y');
       set(H,'FontSize',18);
    end

 otherwise
    error('Unknown Plot(1) Option');
end



%----------------------
%--- Plot residuals ---
%----------------------
switch Plot(2)
 case 0
    % do nothing
 case 1
    % plot residuals
    
    if (size(X,2)>1)
       warning('Can not plot residuals for multi-dimensional fit');
    else
       Ix = 1;
       figure;
       % plot data
       errorxy([X(IndUse,Ix),   Resid(IndUse),   ErrY(IndUse)],...
               'Marker',Marker,'EdgeColor',MarkerColor,...
   	       'FaceColor',FaceColor,'MarkSize',MarkerSize);
       hold on;
       errorxy([X(IndNotUse,Ix),Resid(IndNotUse),ErrY(IndNotUse)],...
               'Marker',OutMarker,'EdgeColor',OutMarkerColor,...
   	       'FaceColor',OutFaceColor,'MarkSize',MarkerSize);


       % set fonts and labels
       set(gca,'FontSize',14);
       H=xlabel(sprintf('X_{%d}',Ix));
       set(H,'FontSize',18);
       H=ylabel('\Delta{Y}');
       set(H,'FontSize',18);
    end
    
 otherwise
    error('Unknown Plot(1) Option');
end



%--------------------------------
%--- Plot residuals histogram ---
%--------------------------------
switch Plot(3)
 case 0
    % do nothing
 case 1
    % plot residuals
    
    if (size(X,2)>1)
       warning('Can not plot residuals for multi-dimensional fit');
    else
       Ix = 1;
       figure;
       % plot data
       StdSize  = std(Resid);
       MaxResid = max(abs(Resid)); 
       Nbins = 2.*max(abs(Resid))./StdSize./HistBinSigma;

       [XU,NU]   = realhist(Resid(IndUse),[-MaxResid MaxResid Nbins]);
       if (isempty(IndNotUse)==1)
	  NnU = zeros(size(NU));
       else
          [XnU,NnU] = realhist(Resid(IndNotUse),[-MaxResid MaxResid Nbins]);
       end

       Hbar = bar(XU,[NU,NnU],'stacked');
       set(Hbar(1),'FaceColor',MarkerColor);
       set(Hbar(2),'FaceColor',OutMarkerColor);

       % set fonts and labels
       set(gca,'FontSize',14);
       H=xlabel('\Delta{Y}');
       set(H,'FontSize',18);
       H=ylabel('Number');
       set(H,'FontSize',18);
    end
    
 otherwise
    error('Unknown Plot(1) Option');
end
