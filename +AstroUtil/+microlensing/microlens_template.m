function [Templates]=microlens_template(T,Pars,varargin)
% SHORT DESCRIPTION HERE
% Package: AstroUtil
% Description: 
% Input  : - Vector of times.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            't0' - Vector of all possible t0 for which to generate
%                   templates.
%            'tE' - Vector of all possible Einstein radius crossing times 
%                   for which to generate templates.
%            'beta' - Vector of all possible beta (minimal impact
%                   parameters) for which to generate templates.
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

if (nargin==0)
    T = (0:10:1000)';
end

DefV.t0                   = (0:1:1000)';
DefV.tE                   = logspace(-3,0,100);   % ER/day
DefV.beta                 = (0.1:0.05:1).^2;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


U2 = In.Beta.^2 + (In.V.*(T - In.T0)).^2;        % Distance^2 [ER units]
BetaT = sqrt(U2);

U2 = BetaT.^2;         % Distance^2 [ER units]
Res.Mu = (U2 + 2)./(sqrt(U2).*sqrt(U2+4));  % Total Magnification

BaseFlux = 1;
Res.Flux = (1-In.Alpha).*BaseFlux  + In.Alpha.*BaseFlux.*Res.Mu;
Mag = In.BaseMag-2.5.*log10(Res.Flux);

Res.Mu1    = 0.5.*Res.Mu + 0.5;
Res.Mu2    = 0.5.*Res.Mu - 0.5;
Res.Theta1 = 0.5.*(Res.BetaT + sqrt(Res.BetaT.^2 + 4));
Res.Theta2 = 0.5.*(Res.BetaT - sqrt(Res.BetaT.^2 + 4));
Res.Theta  = (Res.Theta1.*abs(Res.Mu1) + Res.Theta2.*abs(Res.Mu2))./Res.Mu;
Res.U      = sqrt(U2);
