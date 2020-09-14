function []=fit_flux_powerspec_timedelay(Freq,ObsPS)
% SHORT DESCRIPTION HERE
% Package: AstroUtil
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.A1                  = [];
DefV.A2                  = [];
DefV.Tau                 = (10:0.1:50).';
DefV.gamma               = 3;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


SigmaF = (InPar.A1.^2 + InPar.A2.^2 + 2.*InPar.A1.*InPar.A2.*cos(Omega.*InPar.Tau))./InPar.Freq.^InPar.gamma + InPar.Err.^2;

P_F = -0.5.*log(abs(2.*pi.*SigmaF)) - ObsPS./(2.*SigmaF);


Omega = 2.*pi.*Freq;
Nf = numel(ObsPS);
Ngamma = numel(InPar.gamma);
Ntau   = numel(InPar.Tau);

for Itau=1:1:Ntau
    for Igamma=1:1:Ngamma
        Tau   = InPar.Tau(Itau);
        gamma = InPar.gamma(Igamma);
        
        H = [ones(Nf,2), 2.*cos(Omega.*Tau)]./(Omega.^gamma);
        
        Par = H\ObsPS;
    end
end

        