function YI=interp1_nan(X,Y,varargin)
% Interpolate over NaNs in 1-D vector.
% Package: Util.interp
% Description: Interpolate over NaNs in 1-D vector.
% Input  : - X
%          - Y
%          * Additional parameters to pass to interp1.m.
%            E.g., 'cubic'.
% Output : - Y values at all X values including those at which Y was NaN.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=(1:1:5); Y=[1 2 NaN 4 5]; YI=Util.interp.interp1_nan(X,Y);
% Reliable: 2
%--------------------------------------------------------------------------

YI = Y;
YI(isnan(Y)) = interp1(X(~isnan(Y)), Y(~isnan(Y)), X(isnan(Y)), varargin{:});

