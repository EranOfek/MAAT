function [Mag,Res]=microlens_ps(Pars,T)
% Microlening light curve for a point source as a function of time
% Package: AstroUtil.microlensing
% Description: Calculate the microlensing light curve for a point source,
%              as a function of time (assuming linear motion).
% Input  : - Either vector of free parameters
%            [T0,Beta,V,Alpha,BaseMag]
%            or a matrix of the same free parameters,
%            or a structure array with these fields and parameters.
%            T0 is time of min. impact parameters [days].
%            Beta is the impact parameter in units of the Einstein radius.
%            V is the relative velocity (e.g., Einstein radius per day).
%            Alpha is the blending parameter, 0<Alpha<=1, Alpha=1 if no blending.
%            and BaseMag is the base line magnitude [mag].
%            Default is [0, 0.05, 0.1, 1, 18].
%          - Vector of times. Default is (-20:0.1:20).'.
% Output : - Magnitude for each time.
%          - Structure array of additional parameters for each time.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Feb 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Mag,Res]=AstroUtil.microlensing.microlens_ps; plot(Res.T,Mag);
% Reliable: 2
%--------------------------------------------------------------------------

Def.T    = (-20:0.1:20).';
%         T0,Beta,V,Alpha,BaseMag
Def.Pars = [0, 0.05, 0.1, 1, 18];

if (nargin==0)
    Pars = Def.Pars;
    T    = Def.T;
elseif (nargin==1)
    T    = Def.T;
elseif (nargin==2)
    % do nothing
else
    error('Illegal number of input arguments');
end

if (isstruct(Pars))
    In = Pars;
elseif (isnumeric(Pars))
    In.T0      = Pars(:,1);
    In.Beta    = Pars(:,2);
    In.V       = Pars(:,3);
    In.Alpha   = Pars(:,4);
    In.BaseMag = Pars(:,5);
else
    error('Unknown Pars input');
end
    
U2 = In.Beta.^2 + (In.V.*(T - In.T0)).^2;        % Distance^2 [ER units]
BetaT = sqrt(U2);
[Mag,Res]=AstroUtil.microlensing.microlens_psb(In,BetaT);

% Res.Mu = (U2 + 2)./(sqrt(U2).*sqrt(U2+4));  % Total Magnification
% 
% BaseFlux = 1;
% Res.Flux = (1-In.Alpha).*BaseFlux  + In.Alpha.*BaseFlux.*Res.Mu;
% Mag = In.BaseMag-2.5.*log10(Res.Flux);
% 
% Res.Mu1  = 0.5.*Res.Mu + 0.5;
% Res.Mu2  = 0.5.*Res.Mu - 0.5;
% Res.U    = sqrt(U2);
Res.T    = T;
