function [Mag,Res]=microlens_pspar(Par,JD,varargin)
% Microlening light curve including annual parallax
% Package: AstroUtil.microlensing
% Description: Calculate the microlensing light curve for a point source,
%              as a function of time include annual parallax effects.
% Input  : - A vector of free parameters
%            [T0,ER,RA,Dec,Beta_ra,Beta_dec,Mu_ra,Mu_dec,Par_ls,Alpha,BaseMag].
%            T0 - the JD of minimum Barycentric impact parameter.
%            ER - Einstein radius in consistent units (e.g., mas).
%            RA - J2000 RA of source [rad].
%            Dec - J2000 Dec of source [rad].
%            Beta_ra - Minimum impact parameter in the RA axis in
%                  consistent units (e.g., mas).
%            Beta_dec - Minimum impact parameter in the Dec axis in
%                  consistent units (e.g., mas).
%            Mu_ra - Proper motion in ra in consistent units
%                  (e.g., mas/day).
%            Mu_dec - Proper motion in dec in consistent units
%                  (e.g., mas/day).
%            Par_ls - parallax of lens - parallax of source in consistent
%                  units (e.g., mas).
%            Alpha - the blending parameter, 0<Alpha<=1, Alpha=1 if no blending.
%            BaseMag - the base line magnitude [mag].
%          - Vector of JDs in which to calculate the light curve.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'EpemType' - Which formulae to use in order to calculate
%                         the observer Barycentric position.
%                         'ple' - use ple_earth.m (default).
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
%            .Theta1 - The position of the 2nd image in units of the
%                      Einstein radius, as measured relative to the
%                      lens.
%            .Theta  - The flux weighted combined image position.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: JD = 2451545+[-10:0.01:10].';
%          [Mag,Res]=AstroUtil.microlensing.microlens_pspar([2451545 0.001 4.71 -0.35 0.0001 0.0001 0.0001 0.0001 0.0001 1 18],JD); 
%          [Mag1,Res1]=AstroUtil.microlensing.microlens_pspar([2451545 0.001 4.71 -0.35 0.0001 0.0001 0.0001 0.0001 0.000001 1 18],JD); 
%          plot(JD,Mag); hold on; plot(JD,Mag1,'r-'); plot.invy
% Reliable: 2
%--------------------------------------------------------------------------

DefV.EpemType = 'ple';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

T0       = Par(1);
ER       = Par(2);
RA       = Par(3);
Dec      = Par(4);
Beta_ra  = Par(5);
Beta_dec = Par(6);
Mu_ra    = Par(7);
Mu_dec   = Par(8);
Par_ls   = Par(9);
Alpha    = Par(10);
BaseMag  = Par(11);

Obl = celestial.coo.obliquity(JD);

switch lower(InPar.EpemType)
    case 'ple'
        [Long,~,Rad] = celestial.SolarSys.ple_earth(JD);
        X = Rad.*cos(Long);
        Y = Rad.*sin(Long);
    otherwise
        error('Unknown EpemType option');
end

Alpha_ls = Beta_ra  + (JD - T0).*Mu_ra./cos(Dec)  + Par_ls.*(X.*sin(RA).*sec(Dec) - Y.*cos(RA).*sec(Dec));
Delta_ls = Beta_dec + (JD - T0).*Mu_dec           + Par_ls.*(X.*cos(RA).*sin(Dec) - Y.*(tan(Obl).*cos(Dec) - sin(RA).*sin(Dec)));

BetaT = sqrt((Alpha_ls.*cos(Dec)).^2 + Delta_ls.^2)./ER;
[Mag,Res]=AstroUtil.microlensing.microlens_psb([Alpha, BaseMag],BetaT);
