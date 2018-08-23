function S=int2d(X,Y,Z,Interp)
% Numerical integration of a 2-D matrix
% Package: Util.math
% Description: Numerically interagte a 2-D matrix.
% Input  : - Vector of X values referes to the columns of the matrix
%            to integrate.
%          - Vector of Y values referes to the rows of the matrix
%            to integrate.
%          - A matrix to integrate.
%          - Interpolation method for integration.
%            See interp2.m for options. Default is 'linear'.
%            Note that ExtraVal is set to 0.
% Output : - Matrix numerical integral.
% Tested : Matlab R2011a
%     By : Eran O. Ofek         Dec 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=[1:0.2:10]'; Y=[1:0.1:10]; [MX,MY]=meshgrid(X,Y);           
%          Z=exp(-((MX-5).^2+(MY-5).^2)./2);
%          S=Util.math.int2d(X',Y,Z);
% Reliable: 2
%------------------------------------------------------------------------------
ExtrapVal    = 0;

Def.Interp   = 'linear';
if (nargin==3)
   Interp   = Def.Interp;
elseif (nargin==4)
   % do nothing
else
   error('Illegal number of input arguments');
end

% make sure X and Y are row vectors
if (size(X,1)>size(X,2))
   X = X.';
end
if (size(Y,1)>size(Y,2))
   Y = Y.';
end



% no oversampling is needed
OverX   = X;
OverY   = Y;
InterpZ = Z;


%--- integrate ---
% area of rectangular base
DOX  = diff(OverX);
DOY  = diff(OverY);
Area = bsxfun(@times,DOX.',DOY);

MidOverX  = (OverX(1:end-1)+OverX(2:end)).*0.5;
MidOverY  = (OverY(1:end-1)+OverY(2:end)).*0.5;


MidVal    = interp2(OverX,OverY,InterpZ,...
                    MidOverX.',MidOverY,Interp,ExtrapVal);

MidVolume = Area.*MidVal.';

S         = nansum(MidVolume(:));



