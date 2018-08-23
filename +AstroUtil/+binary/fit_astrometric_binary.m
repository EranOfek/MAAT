function BestFit=fit_astrometric_binary(t,X,Y,ErrX,ErrY,Freq,FitPM,FitPlx,MeanRA,MeanDec)
% Fit an elliptical-orbit binary orbit to astrometric data.
% Package: AstroUtil.binary
% Description: Fit an elliptical-orbit binary orbit to astrometric data.
%              The fit has 7 free parameters, 4 of which (omega, Omega, i,
%              and a) are fitted linearly while two (T, e) are fitted
%              non-linaerly and one (Period) is scanned.        
% Input  : - Column vector of times of the observations (days).
%            If parallax is fitted than this parameter must be in JD.
%          - Column vector of observed X positions.
%          - Coloum vector of observed Y positions.
%          - Column vector of errors in X positions. Default is empty [].
%          - Column vector of errors in Y positions. Default is empty [].
%          - Vector frequencies (1/days) in which to search for
%            binary-period solution. Default is (0.001:0.0001:0.01).'
%            corresponding to search between 100 to 1000 days.
%          - Mean J2000 RA [rad]. Required if fitting parallax.
%          - Mean J2000 Dec [rad]. Required if fitting parallax.
% Output : - Structure containing the best fit parameters. The following
%            fields are available:
%            .Par  - Linear fitted parameters:
%                    [a, omega, Omega, i, MeanPos_x, PM_x, MeanPos_y, PM_y,
%                    Plx].
%            .ParErr - Error in linear fitted parameters.
%            .RefTime - Reference epoch of the proper motion.
%            .El   - Structure of the binary orbital elements.
%            .Chi2 - Chi2 of best fit
%            .Dof  - Number of degrees of freedom.
%            .Npar - Number of free parameters in the model.
%            .ResidX  - Vector of X-axis residuals.
%            .ResidY  - Vector of Y-axis residuals.
%            .VecChi2 - Vector of Chi2 for each trial frequency.
%            .Freq    - Vector of (input) trial frequencies.
%          - Fit proper motion {true|false}. Default is false.
%          - Fit parallax {true|false}. Default is false.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: astrometric_binary.m; chi2_astrometric_binary.m
% Example: [Orbit,Par]=astrometric_binary;
%          Freq  = 1./(330:1:400).';
%          BestFit=AstroUtil.binary.fit_astrometric_binary(Orbit.JD,Orbit.X,Orbit.Y,[],[],Freq);
% Reliable: 2
%--------------------------------------------------------------------------

Def.ErrX = [];
Def.ErrY = [];
Def.Freq = (0.001:0.0001:0.01).';
Def.FitPM  = false;
Def.FitPlx = false;
if (nargin==3)
    ErrX   = Def.ErrX;
    ErrY   = Def.ErrY;
    Freq   = Def.Freq;
    FitPM  = Def.FitPM;
    FitPlx = Def.FitPlx;
elseif (nargin==4)
    ErrY   = Def.ErrY;
    Freq   = Def.Freq;
    FitPM  = Def.FitPM;
    FitPlx = Def.FitPlx;
elseif (nargin==5)
    Freq   = Def.Freq;
    FitPM  = Def.FitPM;
    FitPlx = Def.FitPlx;
elseif (nargin==6)
    FitPM  = Def.FitPM;
    FitPlx = Def.FitPlx;
elseif (nargin==7)
    FitPlx = Def.FitPlx;
end

if (FitPlx)
    if (nargin~=10)
        error('Mean RA/Dec is required while fitting parallax');
    end
else
    MeanRA  = [];
    MeanDec = [];
end
    

if (isempty(Freq)),
    Freq = Def.Freq;
end
Nfreq = numel(Freq);

if (FitPlx)
    Coo = celestial.SolarSys.calc_vsop87(t, 'Earth', 'e', 'E');
else
    Coo = [];
end

Njd = length(t);
GuessPars = [0.0, 0];   % e T
P = 365;
for Ifreq=1:1:Nfreq
    %[Ifreq, Nfreq]
    P = 1./Freq(Ifreq);
    [Fit(Ifreq).X,Fit(Ifreq).Chi2]=Util.fit.fminsearch_my({@chi2_astrometric_binary, P, t, X, Y, ErrX, ErrY, FitPM, FitPlx, MeanRA, MeanDec, Coo},GuessPars);
end

% best fit
[MinChi2,MinInd]=min([Fit.Chi2]);
El.e = Fit(MinInd).X(1);
El.T = Fit(MinInd).X(2);
El.P = 1./Freq(MinInd);

[Chi2,Res] = AstroUtil.binary.chi2_astrometric_binary([El.e, El.T], El.P, t, X, Y, ErrX, ErrY,FitPM,FitPlx,MeanRA,MeanDec,Coo);
Dof        = Njd.*2 - 3 - Res.Npar;
El.a       = Res.El.a;
El.I       = Res.El.I;
El.omega   = Res.El.omega;
El.Omega   = Res.El.Omega;

BestFit.Par     = Res.Par;
BestFit.ParErr  = Res.ParErr;
BestFit.RefTime = Res.RefTime;
BestFit.El      = El;
BestFit.Chi2    = Chi2;
BestFit.Dof     = Dof;
BestFit.Npar    = 3 + Res.Npar;
BestFit.ResidX  = Res.ResidX;
BestFit.ResidY  = Res.ResidY;

BestFit.VecChi2 = [Fit.Chi2].';
BestFit.Freq    = Freq;




