function [Rate,zLim,LimMag,SN]=event_rate(varargin)
% Calculate event rate for ULTRASAT given event and host abs. mag.
% Package: ultrasat
% Description: Calculate event rate for ULTRASAT given event and host
%              abs. mag.
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            ''
% Output : - Event rate per year in FoV.
%          - Redshift detection limit.
%          - Limiting magnitude as a function of z [z, LimMag].
%          - S/N structure for zLim/event rate.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Nov 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Rate,zLim,LimMag,SN]=ultrasat.event_rate
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.EventAbsMag          = -14; %-14.8 + 0.2;
DefV.HostAbsMag           = -14; %-19;
DefV.Rate                 = 3e-5;   % event /Mpc^3/yr

DefV.SN                   = 5;
DefV.FWHM                 = 15;     % FWHM [arcsec]
DefV.PSFeff               = 0.8;    % PSF efficiency
DefV.Aper                 = 33;     % [cm]
DefV.FL                   = 36;     % [cm]
DefV.PixSize              = 9.5;    % [micron]
DefV.RN                   = 3.5;    % [e-]
DefV.DC                   = 1e-2;   % [e-/pix/s]
DefV.Gain                 = 2;      % [e-/ADU]
DefV.ExpTime              = 300;    % [s]
DefV.Nim                  = 3;
DefV.ClearAper            = 0.75;
DefV.Trans                = 1; %0.99.^8 .* 0.99.^4;
DefV.Reflection           = 1; %0.96;
DefV.QE                   = 1; %0.7;
DefV.TargetSpec           = 2e4;    % if given override Mag
DefV.BackSpec             = 23;    % per arcsec^2

DefV.ZodiMagV             = 23.3;
DefV.CerenkovFlux         = 'DailyMax_75flux';
DefV.CerenkovSupp         = 6;
DefV.HostType             = 'Gal_Sc.txt';
DefV.HostMag              = 25;
DefV.HostFilterFamily     = 'SDSS';
DefV.HostFilter           = 'r';
DefV.HostMagSys           = 'AB';

%DefV.HostMag              = 16.5;
DefV.HostFilter           = 'r';
DefV.HostFilterFamily     = 'SDSS';
DefV.CernekovSupp         = 6;

DefV.DetSize              = 9.5;  % [cm]

DefV.Filter               = 'Req4m3';
DefV.FilterFamily         = 'ULTRASAT';
DefV.MagSys               = 'AB';
DefV.Wave                 = (1000:10:25000).';  % [ang]

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


PixScale  = InPar.PixSize.*1e-4./InPar.FL .*RAD.*3600;

PsfEffAreaAS  = 4.*pi.*(InPar.FWHM./2.35).^2;
PsfEffAreaPix = PsfEffAreaAS./(PixScale.^2);

SkyFOV = (InPar.DetSize./InPar.FL).^2./(4.*pi);

SkyFOV = 180./(4.*pi.*RAD.^2);

%        -14.8 (10)       -17.3 (236)
%        -13.5 (1.7)      -15.4 (21)
%        -12.2 (0.3)      -13.6 (2)

AbsMag     = InPar.EventAbsMag;
HostAbsMag = InPar.HostAbsMag;



Vecz = logspace(-3,0,30).';
Nz   = numel(Vecz);
Diff = zeros(Nz,1);
LimMag = [Vecz, zeros(Nz,1)];
for Iz=1:1:Nz
    [LumDist,dM]=AstroUtil.cosmo.lum_dist(Vecz(Iz));
    
    %[~,Vc]=AstroUtil.cosmo.comoving_volume(z);
    %Rate = Vc./1e18 .*3e-5 .* SkyFOV;
    
    AppMag     = AbsMag+dM;  % app mag of event
    HostAppMag = HostAbsMag+dM; % app mag of host


    BackCompFunPar = {'ZodiMagV',InPar.ZodiMagV,'CerenkovFlux',InPar.CerenkovFlux,...
                      'CerenkovSupp',InPar.CerenkovSupp,'HostType',InPar.HostType,...
                      'HostMag',HostAppMag,...
                      'HostFilterFamily',InPar.HostFilterFamily,'HostFilter',InPar.HostFilter,...
                      'HostMagSys',InPar.HostMagSys,'SkyMag',35};

    BackSpec = telescope.sn.back_comp('PsfEffAreaAS',PsfEffAreaAS,BackCompFunPar{:});

    [SN]=telescope.sn.snr('SN',InPar.SN,'FWHM',InPar.FWHM,'PSFeff',InPar.PSFeff,...
                          'Aper',InPar.Aper,'FL',InPar.FL,'PixSize',InPar.PixSize,...
                          'RN',InPar.RN,'DC',InPar.DC,'Gain',InPar.Gain,...
                          'ExpTime',InPar.ExpTime,'Nim',InPar.Nim,...
                          'ClearAper',InPar.ClearAper,...
                          'Trans',InPar.Trans,'Reflection',InPar.Reflection,'QE',InPar.QE,...
                          'TargetSpec',InPar.TargetSpec,...
                          'BackSpec',BackSpec,...
                          'Filter',InPar.Filter,'FilterFamily',InPar.FilterFamily,...
                          'MagSys',InPar.MagSys,'Wave',InPar.Wave);
    Diff(Iz) = SN.LimMag - AppMag;
    
    LimMag(Iz,2) = SN.LimMag;
end

zLim = interp1(Diff,Vecz,0,'pchip');

[LumDist,dM]=AstroUtil.cosmo.lum_dist(zLim);
    
[~,Vc]=AstroUtil.cosmo.comoving_volume(zLim);
Rate = Vc./1e18 .*InPar.Rate .* SkyFOV;

AppMag     = AbsMag+dM;  % app mag of event
HostAppMag = HostAbsMag+dM; % app mag of host


BackCompFunPar = {'ZodiMagV',InPar.ZodiMagV,'CerenkovFlux',InPar.CerenkovFlux,...
                  'CerenkovSupp',InPar.CerenkovSupp,'HostType',InPar.HostType,...
                  'HostMag',HostAppMag,...
                  'HostFilterFamily',InPar.HostFilterFamily,'HostFilter',InPar.HostFilter,...
                  'HostMagSys',InPar.HostMagSys,'SkyMag',35};

BackSpec = telescope.sn.back_comp('PsfEffAreaAS',PsfEffAreaAS,BackCompFunPar{:});

[SN]=telescope.sn.snr('SN',InPar.SN,'FWHM',InPar.FWHM,'PSFeff',InPar.PSFeff,...
                      'Aper',InPar.Aper,'FL',InPar.FL,'PixSize',InPar.PixSize,...
                      'RN',InPar.RN,'DC',InPar.DC,'Gain',InPar.Gain,...
                      'ExpTime',InPar.ExpTime,'Nim',InPar.Nim,...
                      'ClearAper',InPar.ClearAper,...
                      'Trans',InPar.Trans,'Reflection',InPar.Reflection,'QE',InPar.QE,...
                      'TargetSpec',InPar.TargetSpec,...
                      'BackSpec',BackSpec,...
                      'Filter',InPar.Filter,'FilterFamily',InPar.FilterFamily,...
                      'MagSys',InPar.MagSys,'Wave',InPar.Wave);


