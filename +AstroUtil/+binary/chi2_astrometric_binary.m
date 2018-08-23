function [Chi2,Res]=chi2_astrometric_binary(Pars,P,t,X,Y,ErrX,ErrY,FitPM,FitPlx,MeanRA,MeanDec,Coo)
% Find astrometric binary elements from observations
% Package: AstroUtil.binary
% Description: Given the astrometric observations of a binary and its
%              eccentricity, Time of periastron and period, fit linearly
%              the remaining four orbital elements (omega, Omega, i, a).
%              This function is used by fit_astrometric_binary.m
% Input  : - Two elements vector of [e, T].
%          - Astrometric binary trial period.
%          - Vector of time of observations (days).
%            If parallax is fitted than this parameter must be in JD.
%          - Column vector of observed X positions.
%          - Coloum vector of observed Y positions.
%          - Column vector of errors in X positions. Default is empty [].
%          - Column vector of errors in Y positions. Default is empty [].
%          - Fit proper motion {true|false}. Default is false.
%          - Fit parallax {true|false}. Default is false.
%          - Mean J2000 RA [rad]. Required if fitting parallax.
%          - Mean J2000 Dec [rad]. Required if fitting parallax.
%          - Earth barycentric [X, Y, Z] position [au] in ecliptic rectangular
%            coordinates, referred to equinox and ecliptic of J2000.0.
% Output : - Chi2 of best fit.
%          - Structure of best fit parameters.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: fit_astrometric_binary.m; astrometric_binary.m
% Example: [Orbit,Par]=AstroUtil.binary.astrometric_binary;
%          [Chi2,Res]=AstroUtil.binary.chi2_astrometric_binary([365,0.5,1],Orbit.JD,Orbit.X,Orbit.Y,[],[]);
% Reliable: 2
%--------------------------------------------------------------------------
JYear= 365.25;      % Julian year

e = Pars(1);
T = Pars(2);


Def.ErrX   = [];
Def.ErrY   = [];
Def.FitPM  = false;
Def.FitPlx = false;
if (nargin==4)
    ErrX   = Def.ErrX;
    ErrY   = Def.ErrY;
    FitPM  = Def.FitPM;
    FitPlx = Def.FitPlx;
elseif (nargin==5)
    ErrY   = Def.ErrY;
    FitPM  = Def.FitPM;
    FitPlx = Def.FitPlx;
elseif (nargin==6)
    FitPM  = Def.FitPM;
    FitPlx = Def.FitPlx;
elseif (nargin==7)
    FitPlx = Def.FitPlx;       
else
    % do nothing
end

N   = length(t);
if (isempty(ErrX))
   ErrX = ones(N,1);
end
if (isempty(ErrY))
   ErrY = ones(N,1);
end

if (e<1 && e>0)
    
    n = 2.*pi./P;
    M = n.*(t-T);
    a = 1;  % scale a
    [Nu,R] = celestial.Kepler.kepler_elliptic(M,a.*(1-e),e,NaN);

    % fit the linear parameters
    % design matrix
    Vec = [X;Y];
    Err = [ErrX;ErrY];
    Hx   = [R.*cos(Nu), R.*sin(Nu)];
    Hy   = [R.*cos(Nu), R.*sin(Nu)];
    H    = [Hx, zeros(N,2); zeros(N,2), Hy];
    
    Res.RefTime = [];
    if (FitPM)
        Res.RefTime = mean(t);
        Hpmx = [ones(N,1), (t-Res.RefTime)./JYear];
        Hpmy = [ones(N,1), (t-Res.RefTime)./JYear];
        Hpm  = [Hpmx, zeros(N,2); zeros(N,2), Hpmy];
        H    = [H, Hpm];
    end
    if (FitPlx)
        %[Coo,Vel] = calc_vsop87(t, 'Earth', 'e', 'E');
        X    = Coo(1,:).';
        Y    = Coo(2,:).';
        Z    = Coo(3,:).';
        Hplx = [(X.*sin(MeanRA) - Y.*cos(MeanRA)); (X.*cos(MeanRA).*sin(MeanDec) + Y.*sin(MeanRA).*sin(MeanDec) - Z.*cos(MeanDec))];
        H    = [H, Hplx];
    end

    Npar = size(H,2);    
    [Res.Par,Res.ParErr] = lscov(H,Vec,1./(Err.^2));
    
    
    Res.Resid  = Vec - H*Res.Par;
    Res.ResidX = Res.Resid(1:N);
    Res.ResidY = Res.Resid(N+1:end);
    Res.Chi2   = sum((Res.Resid./Err).^2);
    Res.Dof    = 2.*N - Npar;  % 4 is for [a, omega, Omega, i]
    Res.Npar   = Npar;
    Chi2       = Res.Chi2;

    Res.El = celestial.Kepler.thiele_innes2el(Res.Par(1),Res.Par(2),Res.Par(3),Res.Par(4));
else
    % e>=1 return very high chi^2
    Chi2 = 1e9;
    Res = [];
end
