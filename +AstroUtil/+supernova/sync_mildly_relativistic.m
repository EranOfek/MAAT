function nuLnu=sync_mildly_relativistic(t,varargin)
% Synchrotron emission from mildly relativistic ejecta interacting with ISM
% Package: AstroUtil.supernova
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: nuLnu=AstroUtil.supernova.sync_mildly_relativistic
% Reliable: X
%--------------------------------------------------------------------------

if (nargin==0)
    t = [];
end
if (isempty(t))
    t = logspace(0,3,100)'.*86400;
end


DefV.n                   = 1e-2;
DefV.gamma_max           = 1e5;
DefV.Eps_e               = 0.1;
DefV.Eps_B               = 0.01;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


gamma = 1./sqrt(1-beta.^2);
nuLnu = 4.*pi.* InPar.Eps_e .* gamma.^2.*beta.^5.*InPar.n.*constant.mp.*constant.c.^5.*t.^2./( 2.*log(InPar.gamma_max).*(1-beta).^3);

%nu_c = ...