function [Mag,Res]=microlens_psb(Pars,BetaT)
% Microlensing lightcurve for a point source as a function of angular dist
% Package: AstroUtil.microlensing
% Description: Calculate the microlensing light curve for a point source,
%              in a list of a given source-lens distances.
% Input  : - Either vector of free parameters
%            [Alpha,BaseMag]
%            or a matrix of the same free parameters,
%            or a structure array with these fields and parameters.
%            Alpha is the blending parameter, 0<Alpha<=1, Alpha=1 if no blending.
%            and BaseMag is the base line magnitude [mag].
%            Default is [1, 18].
%          - Vector of lens-source distances in units of the Einstein
%            radius. Default is (-2:0.01:2).'.
% Output : - Magnitude for each time.
%          - Structure array of additional parameters for each time.
%            .BetaT  - The source-lens distances for which the parameters
%                      were calculated (Einstein radius units).
%            .Mu     - Total magnification.
%            .Mu1    - 1st image magnification.
%            .Mu2    - 2nd image magnification.
%            .Theta1 - The position of the 1st image in units of the
%                      Einstein radius, as measured relative to the
%                      lens.
%            .Theta2 - The position of the 2nd image in units of the
%                      Einstein radius, as measured relative to the
%                      lens.
%            .Theta  - The flux weighted combined image position.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Mag,Res]=microlens_psb; plot(Res.BetaT,Mag);
% Reliable: 2
%--------------------------------------------------------------------------


Def.BetaT = abs(-2:0.01:2).';
%         Alpha,BaseMag
Def.Pars  = [1, 18];

if (nargin==0)
    Pars  = Def.Pars;
    BetaT = Def.BetaT;
elseif (nargin==1)
    BetaT = Def.BetaT;
elseif (nargin==2)
    % do nothing
else
    error('Illegal number of input arguments');
end

if (isstruct(Pars))
    In = Pars;
elseif (isnumeric(Pars))
    In.Alpha   = Pars(:,1);
    In.BaseMag = Pars(:,2);
else
    error('Unknown Pars input');
end
    
Res.BetaT  = BetaT;

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
