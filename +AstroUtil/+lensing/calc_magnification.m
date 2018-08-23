function [Mag]=calc_magnification(A_11,A_22,A_12)
% Magnification from mapping matrix
% Package: AstroUtil.lensing
% Description: Given the gravitational lensing mapping matrix
%              d\alpha/d\theta, calculate the magnification.
% Input  : - Mapping matrix, A_11, d\alpha/d\theta _11.
%            This can be a scalar, vector or matrix.
%          - Mapping matrix, A_22, d\alpha/d\theta _22.
%            This can be a scalar, vector or matrix.
%          - Mapping matrix, A_12, d\alpha/d\theta _12.
%            This can be a scalar, vector or matrix.
% Output : - Magnification at each \theta.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%------------------------------------------------------------------------------

Mag = 1./((1-A_11).*(1-A_22) - A_12.^2);
