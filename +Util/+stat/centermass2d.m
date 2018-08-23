function M=centermass2d(X,Y,Data);
%------------------------------------------------------------------------------
% centermass2d function                                                 FitFun
% Description: Calculate the center of mass and second moments of
%              2-dimensional matrix.
% Input  : - Matrix or vector containing X coordinates of data.
%          - Matrix or vector containing Y coordinates of data.
%          - Matrix or vector containing data.
% Output : - Structure containing the moments of data:
%            .X     - first moment (barycenter) in X
%            .Y     - first moment (barycenter) in Y
%            .X2    - second moment in X
%            .Y2    - second moment in Y
%            .XY    - second moment in XY
%            .ErrX  - Error in first momet in X.
%            .ErrY  - Error in first momet in Y.
%            .Theta - ellipsoid rotation
%            .R     - radius
%            .A     - major axis
%            .B     - minor axis
%            .Rho   - correlation
% Tested : Matlab 7.0
%     By : Eran O. Ofek                     April 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [MatX,MatY]=meshgrid([1:1:10],[1:1:10]);
%          G=bivar_gauss(MatX,MatY,[4.834 5.12 3 3 0]);
%          N = 1000; B=100;
%          [M]=center2d_cm(MatX,MatY,sqrt(G.*N+B));
% Notes: previously called: enter2d_cm.m
% Reliable: 2
%------------------------------------------------------------------------------

Sum      = sumnd(Data);
SizeData = size(Data);

% barycenter
M.X    = sumnd(Data.*X)./Sum;
M.Y    = sumnd(Data.*Y)./Sum;

% 2nd moments
M.X2   = sumnd(Data.*X.^2)./Sum - M.X.^2;
M.Y2   = sumnd(Data.*Y.^2)./Sum - M.Y.^2;
M.XY   = sumnd(Data.*X.*Y)./Sum - M.X.*M.Y;

Npix = prod(SizeData);
M.ErrX = sqrt(M.X2./Npix);
M.ErrY = sqrt(M.Y2./Npix);

Tan2Theta = 2.*M.XY./(M.X2 - M.Y2);
M.Theta    = 0.5.*atan(Tan2Theta);
R2        = 0.5.*(M.X2 + M.Y2);
M.R       = sqrt(R2);
Diff2     = 0.25.*(M.X2 - M.Y2).^2 - M.XY.^2;
M.A       = sqrt(R2 + sqrt(Diff2));
M.B       = sqrt(R2 - sqrt(Diff2));

%--- The correlation coef. ---
M.Rho = Tan2Theta.*(M.A.^2 - M.B.^2)./(2.*M.A.*M.B);



