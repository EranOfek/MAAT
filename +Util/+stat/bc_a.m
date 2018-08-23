function [A1,A2]=bc_a(Z0,Ac,Percentile);
%------------------------------------------------------------------------------
% bc_a function                                                      AstroStat
% Description: Bias correction and acceleration for bootstrap and
%              Jackknife estimetors of the confidence interval.
%              The bc_a bias coorection and acceleration can be used
%              to estimate the bias corrected and accelerated confiedence
%              interval (CI).
% Input  : - bc_a bias correction (Z0).
%          - bc_a acceleration (Ac).
% Output : - The new lower (A1) percentile to plug in err_cl.m in order
%            to get the bias corrected and accelerated CI.
%          - (A2) Like A1, but for the upper percentile.
% Reference : Efron B. & Tibshirani, R.J., 1993,
%             in: An Introduction to the Bootstrap
% see also: bootstrap_std.m, bc_a.m, err_cl.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                   October 2002
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: % Looking for the CI for the 0.9545 CL:
%          [A1,A2]=bc_a(Z0,Ac,0.9545);  % estimate "corrected" percentile,
%          CI=err_cl(Theta,0.9545);
% Reliable: 1
%------------------------------------------------------------------------------
Za   = norminv(Percentile,0,1);
Z1a  = norminv(1-Percentile,0,1);

Arg1 = Z0 + (Z0 + Za) ./(1 - Ac.*(Z0 + Za ));
Arg2 = Z0 + (Z0 + Z1a)./(1 - Ac.*(Z0 + Z1a));

A2   = normcdf(Arg1);
A1   = normcdf(Arg2);
