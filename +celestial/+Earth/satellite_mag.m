function [Res]=satellite_mag(varargin)
% Satellite apparent magnitude
% Package: celestial
% Description: Satellite apparent magnitude
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Dist' - Satellite distance [km]. Default is 500.
%            'Area' - Satellite area [cm^2]. Default is 1.
%            'Albedo' - Albedo. Default is 0.1.
%            'SunAbsMag' - Sun abs. mag. Default is -26.74 (V Vega).
% Output : - Structure containing the following fields:
%            'Mag' - Satellite apparent magnitude.
%            'AngV' - Satellite angular velicity assuming Keplerian orbit
%                     with no projection [arcsec/s].
%            'MagResEl' - Magnitude per resolution element (FWHM-length).
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Dec 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=celestial.Earth.satellite_mag
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.Dist                 = 500;     % km
DefV.Area                 = 1;       % [cm^2]
DefV.Albedo               = 0.1;
DefV.FWHM                 = 3;
DefV.SunAbsMag            = -26.74;  % V
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Res.Mag = InPar.SunAbsMag - 2.5.*log10(InPar.Area ./ (4.*pi.*(InPar.Dist.* 1e5)^2)) - 2.5.*log10(InPar.Albedo);

K = celestial.Kepler.kepler3law(constant.EarthM,'a',constant.EarthR+InPar.Dist.*1e5);
Res.AngV = K.v./(InPar.Dist.*1e5).*RAD.*3600;  % ''/s

Res.MagResEl = Res.Mag + 2.5.*log10(Res.AngV);
