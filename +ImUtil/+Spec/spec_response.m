function Tran=spec_response(StdObsSpec,Std,varargin)
%--------------------------------------------------------------------------
% spec_response function                                            ImSpec
% Description: Given the observed spectrum of a flux standard, and the
%              flux standard calibrated spectrum outside the atmosphere,
%              attempt to find the transmission function of the
%              spectrograph and telescope.
% Input  : - Observed spectrum of a standard star [Wave(Ang), Intensity].
%          - The standard star flux calibrated spectrum outside the
%            atmosphere. This can be either:
%            (1) A matrix of [Wavelength(Ang), Flux(erg/cm^2/s/Ang)].
%            (2) A structure which was returned by search_specphot_stand.m
%            (3) A string which contain standard star name
%                (e.g., 'BD+33 2642'),
%            (4) A cell array of coordinates, in which the first element
%                is the std J2000 RA in radians, or a sexagesimal string,
%                the second element is the std J2000 Dec in radians, or a
%                sexagesimal string, and the third element is the search
%                radius in arcsec. Default for the third element is 60.
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'AM'    - Airmass at which the standard star was observed.
%                      Default is []. If empty then use the Altitude
%                      in order to calculate the airmass.
%            'Alt'   - Altitude (radians) at which the standard star was
%                      observed. Default is pi./2.
%                      If both Altitude and Airmass are provided then
%                      the Airmass will be used.
%            'ExpTime'- Exposure time (s) of standard star observations.
%                      This is required if the transmission is needed
%                      in non arbitrary units.
%                      Default is 1.
%            'Ext'   - Atmospheric extinction type. See atmospheric_ext.m
%                      for options. Default is
%                      'KPNO_atmospheric_extinction.dat'.
%                      If empty then will not apply atmospheric extinction
%                      correction.
%            'InterpMethod' - Interpolation method. See interp1.m for
%                      options. Default is 'linear'.
%            'TranInterpMethod' - Interpolation method for the transmission
%                      calculation from WaveGrid. See interp1.m for
%                      options. Default is 'cubic'.
%            'WaveGrid' - Grid points of wavelengths [Ang] which will be
%                      used as ancor points in the interpolation of
%                      the spectra over emission/absorbtion lines.
%                      This points are returned by search_specphot_stand.m
%                      If empty matrix, and WaveGrid is not available
%                      in the Standard star structure returned by
%                      search_specphot_stand.m, then the user will be
%                      prompted to select this points manually.
%                      If 'manual' then will prompt the user to select the
%                      points manually.
%                      This can also be used to override the grid points
%                      returned by search_specphot_stand.m
%                      Default is empty.
%            'SmoothObs' - Smooth observed spectrum before using it.
%                      options are {'none','medfilt1',runmean1'}.
%                      Default is 'none'.
%            'R'     - The resolution of the smoothed spectrum.
%                      Default is 500. If empty then do not expand
%                      Telluric lines according to spectral resolution.
%            'Filter'- Filter to use in the spectral degradation by
%                      the resolution. Default is 'mean'.
%                      See spec_convolve_resolution.m for options.
%            'TranMethod' - Calculating the response tranformation in
%                      linear flux or log flux {'linear'|'log'}.
%                      Default is 'log'.
%            'Plot'  - Plot response fit:
%                      'y' - plot
%                      'n' - do not plot (default).
%                      Otherwise this will be a file name and an eps plot
%                      will be saved as this file.
% Output : - Structure containing information about the atmosphere
%            transmission function. The following fields are available:
%            .VecWave - Wavelength at which the inverse transmission
%                       is evaluated.
%            .InvTran - inverse transmission required to convert a 1 s
%                       integration
%                       spectrum from which the atmospheric extinction 
%                       was removed (i.e., outside the atmosphere) to
%                       flux in erg^-1 cm^-2 s^-1 Ang^-1. This has units
%                       of [ADU erg^-1 cm^2 Ang].
%            .InvTranS - A smoothed version of .InvTran.
%                       This is the recomended response function to use.
%            .AM      - Airmass of the standard star observation.
%            .StdAirlessSpec - The observed Std spectra converted
%                       to airless extinction.
%            .StdCalSpec - Stand star calibrated spectrum.
%            .GridWave - The GridWave used.
%            .Resid    - Residuals vector between the observed and
%                        calibrated spectra in units of the calibrated
%                        spectra.
%            .Ratio    - Ratio vector between the observed and
%                        calibrated spectra in units of the calibrated
%                        spectra.
%            .FlagTelluric - Flag indicating if each one of the observed
%                        spectra wavelengths is within a Telluric
%                        absorption. 0 for no, >0 for yes.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Sep 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 2
%--------------------------------------------------------------------------

Def.SearchRad = 60;  % arcsec


DefV.AM               = [];
DefV.Alt              = pi./2;
DefV.ExpTime          = 1;
DefV.Ext              = 'KPNO_atmospheric_extinction.dat';
DefV.R                = 500;        % approximate resultion of observed spectrum
DefV.Filter           = 'mean';     % degrading to resultion method
DefV.InterpMethod     = 'linear';
DefV.GridWave         = [];
DefV.OverSample       = -4;
DefV.SmoothTran       = 'median';   % {'median','mean','gauss','poly','no'}
DefV.SmoothDeg        = 5;             
DefV.SmoothR          = 50;         % smoothing resolution
DefV.SmoothLog        = true;       % smooth in log space
% 
% DefV.TranInterpMethod = 'linear'; %'cubic';
% 
% DefV.SmoothObs        = 'none';
% DefV.R                = 500;
% DefV.TranMethod       = 'log';
% DefV.Plot             = 'y';
% DefV.FitMethod        = 'interp';  % interp | poly
% DefV.ResponseDeg      = 7;
% %DefV.SmoothTran       = 'no';   % no | poly
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

% set AirMass
if (isempty(InPar.AM)),
    InPar.AM = hardie(pi./2-InPar.Alt);
end

% correct spectrum for exposure time
% spectrum in 1s exposure
StdObsSpec(:,2) = StdObsSpec(:,2)./InPar.ExpTime;  % ADU/s

%--- Atmospheric extinction ---
% If required, correct the observed std spectrum from
% extincted to airles spectrum
if (isempty(InPar.Ext)),
    % no atmospheric extinction correction - already airless
    StdObsAirlessSpec = StdObsSpec;
else
    StdObsAirlessSpec = atmospheric_ext(StdObsSpec,InPar.AM,InPar.Ext,InPar.InterpMethod);
end



%--- vector of wavelength ---
VecWave = StdObsAirlessSpec(:,1);
Nwave   = numel(VecWave);

%--- oversample ---
%   either a negaive scalar, posative scalar, vector or empty
if (~isempty(InPar.OverSample)),
     if (numel(InPar.OverSample)==1),
         if (InPar.OverSample>0),
             % linear sampling
             VecWave = linspace(min(VecWave),max(VecWave),ceil(Nwave.*InPar.OverSample)).';
         else
             % logarithmic sampling
             VecWave = logspace(log10(min(VecWave)),log10(max(VecWave)),ceil(Nwave.*abs(InPar.OverSample))).';
         end
     else
         % assume OverSample is a vector of grid points
         VecWave = InPar.OverSample;
     end
     StdObsAirlessSpec = [VecWave, interp1(StdObsAirlessSpec(:,1),StdObsAirlessSpec(:,2),VecWave,InPar.InterpMethod)];          
end


%--- Extend Telluric region to match spectral resolution ---
Telluric = load2('Telluric.mat');
FlagTelluric   = Util.array.find_ranges_flag(VecWave,Telluric.Cat(:,1:2));
OutSpec = spec_convolve_resolution([VecWave, FlagTelluric],InPar.R,'mean');
FlagTelluric   = OutSpec(:,2)>0;

%--- read Std spectrum ---
% Deal with std formats:
%            (1) A matrix of [Wavelength(Ang), Flux(erg/cm^2/s/Ang)].
%            (2) A structure which was returned by search_specphot_stand.m
%            (3) A string which contain standard star name
%                (e.g., 'BD+33 2642'),
%            (4) A cell array of coordinates, in which the first element
%                is the std J2000 RA in radians, or a sexagesimal string,
%                the second element is the std J2000 Dec in radians, or a
%                sexagesimal string, and the third element is the search
%                radius in arcsec. Default for the third element is 60.
if (isnumeric(Std)),
    StdCalSpec = Std;
    GridWave   = [];
elseif (isstruct(Std)),
    StdCalSpec = Std.Spec;
    GridWave   = Std.GridWave;
elseif (ischar(Std)),
    StdStruct  = search_specphot_stand(Std);
    StdCalSpec = StdStruct.Spec;
    GridWave   = StdStruct.GridWave;
elseif (iscell(Std)),
    % search by coordinates
    RA  = Std{1};
    Dec = Std{2};
    if (numel(Std)<3),
        SearchRad = Def.SearchRad;
    else
        SearchRad = Std{3};
    end
    StdStruct  = search_specphot_stand(RA,Dec,SearchRad);
    StdCalSpec = StdStruct.Spec;
    GridWave   = StdStruct.GridWave;
else
    error('Unknown Std type');
end

% define GridWave
if (isempty(InPar.GridWave) && isempty(GridWave)),
    InPar.GridWave = 'manual';
end

if (~isempty(InPar.GridWave)),
    % override GridWave from stand star DB with user supplied GridWave
    GridWave       = InPar.GridWave;
end

%--- Interpolate over Telluric regions ---
% observed spectrum only
StdObsAirlessSpecI = StdObsAirlessSpec;
StdObsAirlessSpecI(:,2) = interp1(VecWave(FlagTelluric==0),StdObsAirlessSpec(FlagTelluric==0,2),VecWave,InPar.InterpMethod,'extrap');

%--- Degrade observed and stdcal spectra ---
DegradedObs = spec_convolve_resolution(StdObsAirlessSpecI,InPar.R,InPar.Filter);
DegradedCal = spec_convolve_resolution(StdCalSpec,InPar.R,InPar.Filter);

% equalize the sampling of the two spectra
[DegradedObs,DegradedCal]=eq_sampling(DegradedObs,DegradedCal,DegradedObs(:,1),InPar.InterpMethod);


Tran.VecWave           = VecWave;
Tran.InvTran           = DegradedCal(:,2)./DegradedObs(:,2);
Tran.FlagTelluric      = FlagTelluric;
Tran.StdCalSpec        = StdCalSpec;
Tran.StdCalSpecInt     = interp1(StdCalSpec(:,1),StdCalSpec(:,2),Tran.VecWave,InPar.InterpMethod);
Tran.StdObsAirlessSpec = StdObsAirlessSpec;

%--- Smooth transmission ---
switch lower(InPar.SmoothTran)
    case {'no','none'}
        Tran.InvTranS = Tran.InvTran;
    case {'median','mean','gauss'}
        OutSpec = spec_convolve_resolution([VecWave,Tran.InvTran],InPar.SmoothR,InPar.SmoothTran);
        Tran.InvTranS = OutSpec(:,2);
    case 'poly'
        if (InPar.SmoothLog),
            % fit polynomial in log space
            Par = polyfit(VecWave(FlagTelluric==0),log10(Tran.InvTran(FlagTelluric==0)),InPar.SmoothDeg);
            Tran.InvTranS = 10.^(polyval(Par,VecWave));
       
        else
            % fit polynomial in linear space
            Par = polyfit(VecWave(FlagTelluric==0),Tran.InvTran(FlagTelluric==0),InPar.SmoothDeg);
            Tran.InvTranS = polyval(Par,VecWave);
            
        end
    otherwise
        error('Unknown SmoothTran option');
end

Tran.Diff  = Tran.StdObsAirlessSpec(:,2).*Tran.InvTranS - Tran.StdCalSpecInt;
Tran.Ratio = Tran.StdObsAirlessSpec(:,2).*Tran.InvTranS./Tran.StdCalSpecInt;


