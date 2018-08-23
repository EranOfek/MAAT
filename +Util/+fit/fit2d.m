function [PX,PY,PeX,PeY,Chi2X,Chi2Y,CovX,CovY,Freedom,ResX,ResY]=fit2d(RefXY,DataXY,DelDataXY,XY_0);
%------------------------------------------------------------------------------
% fit2d function                                                        FitFun
% Description: Find 2D geometrical transformation between two sets of
%              x/y coordinates. The first set is the reference set and
%              the second is the data. The fit is of the form:
%                      Y' = dy + ay0(y-y0)   + by0(x-x0)   +
%                                ay1(y-y0)^2 + by1(x-x0)^2 +
%                                cy(x-x0)(y-y0)
%                      X' = dx + ax0(x-x0)   + bx0(y-y0)   +
%                                ax1(x-x0)^2 + bx1(y-y0)^2 +
%                                cx(x-x0)(y-y0)
%                      Where x/y are the reference coordinates
%                      and X'/Y' are the data coordinates.
%                      x0/y0 are coordinates defined by the user.
%                      a0/a1/b0/b1/c/d are the model parameters.
% Input  : - Two column matrix of reference coordinates x and y.
%          - Two column matrix of data coordinates x and y.
%          - Two column matrix of errors in data coordinates x and y.
%            In case that two elements vector is given,
%            the first element is taken as the error in x and the
%            second element for the error in y. 
%          - [x0,y0], default is [0 0].
% Output : - Fitted parameters for X:  [dx ax0 bx0 ax1 bx1 cx]
%          - Fitted parameters for Y:  [dy ay0 by0 ay1 by1 cy]
%          - Fitted error parameters for X:  [dx ax0 bx0 ax1 bx1 cx]
%          - Fitted error parameters for Y:  [dy ay0 by0 ay1 by1 cy]
%          - Chi2 of the X fit.
%          - Chi2 of the Y fit.
%          - Degrees of freedom.
%          - The covariance matrix for X.
%          - The covariance matrix for Y.
%          - The X axis residuals vector. [calculated error for Chi^2=1,
%            can be calculated from mean(abs(Resid)) ].
%          - The Y axis residuals vector. 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                     March 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------------------

Deg = 6;
if (nargin==3),
   XY_0 = [0 0];
elseif (nargin==4),
   % do nothing
else
   error('Illigal number of input arguments');
end

X_0 = XY_0(1);
Y_0 = XY_0(2);


N_R  = length(RefXY(:,1));
N_D  = length(DataXY(:,1));
N_DD = length(DelDataXY(:,1));

if (N_DD==1),
   DelData(:,1) = DelDataXY(1).*ones(N_D,1);
   DelData(:,2) = DelDataXY(2).*ones(N_D,1);
else
   DelData = DelDataXY;
end

if (N_R~=N_D),
   error('number of reference point and data point should be the same');
end


% degree of freedom
Freedom = N_D - Deg;
   
% building the H matrix for: X
HX = zeros(N_D,Deg);
HX(:,1) = ones(N_D,1);
HX(:,2) = RefXY(:,1)-X_0;
HX(:,3) = RefXY(:,2)-Y_0;
HX(:,4) = (RefXY(:,1)-X_0).^2;
HX(:,5) = (RefXY(:,2)-Y_0).^2;
HX(:,6) = (RefXY(:,1)-X_0).*(RefXY(:,2)-Y_0);

   
% building the Covariance matrix for: X
VX = diag(DelData(:,1).^2);
   
CovX   = inv(HX.'*inv(VX)*HX);
PX     = CovX*HX.'*inv(VX)*DataXY(:,1);
PeX    = sqrt(diag(CovX));

%'Number of degree of freedom :', Freedom
ResX   = DataXY(:,1) - HX*PX;
Chi2X  = sum((ResX./DelData(:,1)).^2);

%Chi2/Freedom
%sqrt(2/Freedom)




% building the H matrix for: Y
HY = zeros(N_D,Deg);
HY(:,1) = ones(N_D,1);
HY(:,2) = RefXY(:,2)-Y_0;
HY(:,3) = RefXY(:,1)-X_0;
HY(:,4) = (RefXY(:,2)-Y_0).^2;
HY(:,5) = (RefXY(:,1)-X_0).^2;
HY(:,6) = (RefXY(:,1)-X_0).*(RefXY(:,2)-Y_0);

   
% building the Covariance matrix for: X
VY = diag(DelData(:,2).^2);
   
CovY   = inv(HY.'*inv(VY)*HY);
PY     = CovY*HY.'*inv(VY)*DataXY(:,2);
PeY    = sqrt(diag(CovY));

%'Number of degree of freedom :', Freedom
ResY   = DataXY(:,2) - HY*PY;
Chi2Y  = sum((ResY./DelData(:,2)).^2);

%Chi2/Freedom
%sqrt(2/Freedom)



