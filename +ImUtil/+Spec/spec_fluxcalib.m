function Tran=spec_fluxcalib(StdObsSpec,Std,varargin)
%--------------------------------------------------------------------------
% spec_fluxcalib function                                           ImSpec
% Description: Given the observed spectrum of a flux standard, and the
%              flux standard calibrated spectrum outside the atmosphere,
%              attempt to find the transmission function of the
%              spectrograph and telescope.
% Input  : - Observed spectrum of a standard star [Wave(Ang), Intensity].
%          - The standard star flux calibrated spectrum outside the
%            atmosphere. This can be either: a matrix of
%            [Wavelength(Ang), Flux(erg/cm^2/s/Ang)].
%            or a structure which was returned by search_specphot_stand.m
%            or a string which contain standard star name
%            (e.g., 'BD+33 2642').
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
%                      This points are returned (for at least some stars)
%                      by search_specphot_stand.m
%                      If empty matrix, and WaveGrid is not available
%                      in the Standard star structure returned by
%                      search_specphot_stand.m, then the user will be
%                      prompted to select this points manually.
%                      If 'manual' then will prompt the user to select the
%                      points manually.
%                      This can also be used to overrid the grid points
%                      returned by search_specphot_stand.m
%                      Default is empty.
%            'SmoothObs' - Smooth observed spectrum before using it.
%                      options are {'none','medfilt1',runmean1'}.
%                      Default is 'none'.
%            'R'     - The resolution of the smoothed spectrum.
%                      Default is 500.
% Output : - Structure containing information about the transmission
%            function. The following fields are available:
%            .Wave    - Wavelength at which the inverse transmission
%                       is evaluated.
%            .InvTran - inverse transmission required to convert a
%                       spectrum from which the atmospheric extinction 
%                       was removed (i.e., outside the atmosphere) to
%                       flux in erg^-1 cm^-2 s^-1 Ang^-1. This has units
%                       of [ADU erg^-1 cm^2 Ang].
%            .AM      - Airmass.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Sep 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 2
%--------------------------------------------------------------------------



DefV.AM               = [];
DefV.Alt              = pi./2;
DefV.ExpTime          = 1;
DefV.Ext              = 'KPNO_atmospheric_extinction.dat';
DefV.InterpMethod     = 'linear';
DefV.TranInterpMethod = 'cubic';
DefV.GridWave         = [];
DefV.SmoothObs        = 'none';
DefV.R                = 500;

InPar = set_varargin_keyval(DefV,'y','use',varargin{:});


if (isempty(InPar.AM)),
    InPar.AM = hardie(pi./2-InPar.Alt);
end

% correct spectrum for exposure time
% spectrum in 1s exposure
StdObsSpec(:,2) = StdObsSpec(:,2)./InPar.ExpTime;  % ADU/s

% FFU: smooth spectra
switch lower(InPar.SmoothObs)
    case 'none'
        % do nothing
    otherwise
        % resample in log space to have a uniform resultion sampling
        SampleWave = logspace(log10(min(StdObsSpec(:,1))),...
                              log10(max(StdObsSpec(:,1))),...
                              length(StdObsSpec(:,1)).*2);
        SampleWaveInt = interp1(StdObsSpec(:,1),StdObsSpec(:,2),SampleWave,InPar.InterpMethod);
        FilterSize = ceil(SampleWave(1)./(SampleWave(2)-SampleWave(1))./InPar.R);
        
        switch lower(InPar.SmoothObs)
            case 'medfilt1'              
                SampleWaveInt = medfilt1(SampleWaveInt,FilterSize);
            case 'runmean1'
                SampleWaveInt = runmean1(SampleWaveInt,floor(FilterSize./2));
            otherwise
                error('Unknown SmoothObs option');
        end
        % resample to original sampling
        StdObsSpec(:,2) = interp1(SampleWave,SampleWaveInt,StdObsSpec(:,1),InPar.InterpMethod);
end

        
        



if (isempty(InPar.Ext)),
    % no atmospheric extinction correction - already airless
    StdAirlessSpec = StdObsSpec;
else
    StdAirlessSpec = atmospheric_ext(StdObsSpec,InPar.AM,InPar.Ext,InPar.InterpMethod);
end

if (isnumeric(Std)),
    StdCalSpec = Std;
    GridWave   = [];
elseif (isstruct(Std)),
    StdCalSpec = Std.Spec;
    GridWave   = Std.GridWave;
elseif (ischar(Std)),
    StdStruct = search_specphot_stand(Std);
    StdCalSpec = StdStruct.Spec;
    GridWave   = StdStruct.GridWave;
else
    error('Unknown Std type');
end

if (isempty(InPar.GridWave) & isempty(GridWave)),
    InPar.GridWave = 'manual';
end
     
if (ischar(InPar.GridWave)),
    if (strcmpi(InPar.GridWave,'manual')),
        % select grid points manually
        Res=plot_int({StdCalSpec(:,1),StdCalSpec(:,2),'-',{},'X','Y'});
        waitfor(gcf,'KeyPressFcn','');

        GridWave = StdCalSpec(Res.IndRM,1);
    end
end

% resample the spectra at WaveGrid
StdAirlessSpec_GW = interp1(StdAirlessSpec(:,1),StdAirlessSpec(:,2),GridWave,InPar.InterpMethod);
StdCalSpec_GW     = interp1(StdCalSpec(:,1),StdCalSpec(:,2),GridWave,InPar.InterpMethod);

Tran.Wave         = StdObsSpec(:,1);

Tran.InvTran      = interp1(GridWave,StdCalSpec_GW./StdAirlessSpec_GW,StdObsSpec(:,1),InPar.TranInterpMethod);
Tran.AM           = InPar.AM;

Tran.StdAirlessSpec = StdAirlessSpec(:,2);
        

    
    
