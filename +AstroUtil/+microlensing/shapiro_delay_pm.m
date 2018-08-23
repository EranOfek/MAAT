function [Dt,Dt_dt1,Dt_dt2,Dt_dt3]=shapiro_delay_pm(t,t0,bm,mu,M)
% The Shaprio time delay given proper motion between two stars.
% Package: AstroUtil.microlensing
% Description: The Shaprio time delay given proper motion between two stars.
% Input  : - time [Julian years].
%          - time of minimum impact parameter [Julian years].
%          - minimum impact parameter [arcsec].
%          - relative proper motion [mas/yr].
%          - Lens mass [solar mass].
% Output : - Time delay [s].
%          - 1st time derivative of time delay.
%          - 2nd time derivative.
%          - 3rd time derivative.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% [Dt,Dt_dt1,Dt_dt2,Dt_dt3]=AstroUtil.microlensing.shapiro_delay_pm(t,t0,bm,mu,M)
% Reliable: 
%--------------------------------------------------------------------------
RAD        = 180./pi;
JYear      = 365.25;
SEC_DAY    = 86400;
ARCSEC_DEG = 3600;

A = -2.*constant.G.*constant.SunM.*M./(constant.c.^3);

T  = t - t0;
T  = T.*JYear.*SEC_DAY;
bm = bm./(RAD.*ARCSEC_DEG); 
mu = mu./(RAD.*ARCSEC_DEG.*1000.*JYear.*SEC_DAY); 

b = sqrt(bm.^2 + (mu.*T).^2);

Dt = A*log(1-cos(b));

% derivatives
% syms b bm mu t t0 A T
% b = sqrt(bm^2 + (mu*T)^2);
% dt = A*log(1-cos(b));
% vectorize(simplify(diff(dt,T,1)))

if (nargout>1)
    Dt_dt1 = -(A.*T.*mu.^2.*sin((bm.^2 + T.^2.*mu.^2).^(1./2)))./((bm.^2 + T.^2.*mu.^2).^(1./2).*(cos((bm.^2 + T.^2.*mu.^2).^(1./2)) - 1));
    if (nargout>2)
        Dt_dt2 = -(A.*mu.^2.*(bm.^2.*sin((bm.^2 + T.^2.*mu.^2).^(1./2)) - T.^2.*mu.^2.*(bm.^2 + T.^2.*mu.^2).^(1./2)))./((bm.^2 + T.^2.*mu.^2).^(3./2).*(cos((bm.^2 + T.^2.*mu.^2).^(1./2)) - 1));
        if (nargout>3)
            Dt_dt3 = (A.*T.*mu.^4.*(3.*bm.^2.*cos((bm.^2 + T.^2.*mu.^2).^(1./2)).*(bm.^2 + T.^2.*mu.^2).^(1./2) - 3.*bm.^2.*sin((bm.^2 + T.^2.*mu.^2).^(1./2)) - 3.*bm.^2.*(bm.^2 + T.^2.*mu.^2).^(1./2) + T.^4.*mu.^4.*sin((bm.^2 + T.^2.*mu.^2).^(1./2)) + 3.*bm.^2.*cos((bm.^2 + T.^2.*mu.^2).^(1./2)).*sin((bm.^2 + T.^2.*mu.^2).^(1./2)) + T.^2.*bm.^2.*mu.^2.*sin((bm.^2 + T.^2.*mu.^2).^(1./2))))./((bm.^2 + T.^2.*mu.^2).^(5./2).*(cos((bm.^2 + T.^2.*mu.^2).^(1./2)) - 1).^2);
        end
    end
end






