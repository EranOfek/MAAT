function []=stoch_time_delay(varargin)
% SHORT DESCRIPTION HERE
% Package: AstroUtil.microlensing
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------



DefV. = 
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

%%
RAD = 180./pi;
Nsim = 1e5;
XY = rand(Nsim,2).*2 - 1;
XY = XY./RAD.*5;
PM = randn(Nsim,2);
dt = 1;
G = constant.G;
M = constant.SunM;
c = constant.c;
Beta0 = XY;
Beta1 = XY + PM./(RAD.*3600.*100) .*dt;

Beta0 = sqrt(sum(Beta0.^2,2));
Beta1 = sqrt(sum(Beta1.^2,2));

sum( 2.*G.*M./(c.^3) .* log(1 - Beta0) ) - sum( 2.*G.*M./(c.^3) .* log(1 - Beta1) )
