function [FluxCalibTran,Par,WaveCalib]=spec_fluxcal(StdImage,StdName,ArcImage,ArcType,varargin);
%--------------------------------------------------------------------------
% spec_fluxcal function                                             ImSpec
% Description: Given a standar star spectrum calculate the sensitivity
%              function of the CCD needed in order to convert
%              an observed spectrum (at airmass zero) to flux.
% Input  : - FITS image of standard star (preferably after bias
%            subtraction and flat field correction).
%          - String containing standard star name, or alternatively,
%            a two column matrix containing the standard star calibrated
%            spectrum [Wavelength(A), Flux(units)]
%            If empty, then attempt to get standard
%            star name from the coordinates in the FITS header (RA and Dec
%            header keywords). See get_specstand.m.
%          - FITS image of a wavelength calibration arc.
%            If empty (default), then attempt to use sky lines for
%            wavelength calibration.
%          - Type of arc image (see get_arclines.m for more details).
%            Options are:
%            {'Ar'|'Cd'|'Hg'|'Ne'|'Zn'|'SA'|'SS'}, default is 'SS'.
%          * Arbitrary number of pair of arguments: ...,keyword,value,...,
%            where keywod can be one of the followings:
%            'ExpKey'    - Exposure time image header keyword, default
%                          is 'EXPOSURE'.
%            'DispAxis'  - Dispersion axis {'x'|'y'}, default is 'x'.
%            'AutoSel'   - Automatically select the position in which to
%                          start the spectrum tracing {'y'|'n'},
%                          default is 'n'. NOT YET AVAILABLE.
%            'Int'       - Interactive mode {'y'|'n'}. If 'n', then try
%                          to use non-interactive mode whenever possible.
%                          Default is 'y'.
%            'InterpMethod' - Interpolation method (see interp1.m),
%                          default is 'linear'.
%            'KeyAirMass'- String containing the header keyword name
%                          for airmass, or alternatively the airmass
%                          of the Standard star spectrum.
%                          Default is 'AIRMASS'.
%            'Ext'       - Extinction type. This can be a file name or a
%                          matrix containing the atmospheric extinction
%                          curve. The format of the file or the matrix is: 
%                          [wavelength(Ang), Extinction_in_unit_airmass(mag)].
%                          Examples: 'VLT_atmospheric_extinction.dat'
%                                    'KPNO_atmospheric_extinction.dat'
%                          default is: 'VLT_atmospheric_extinction.dat'.
% Output : - The transformation needed in order to flux calibrate the
%            observation [Wavelength(A), Factor] from the DN units 
%            as observed at airmass=0, to the units of the flux calibrated
%            standard star. Where facrtor is the number that you need to
%            multiply the observed (after transforming to airmass=0)
%            signal to flux. The transformation is sampled the same as
%            the input spectrum.
%          - Polynomial parmeters (see polyfit_sc.m) of the best fit
%            transformation.
%          - Structure containing the wavelength calibration for the
%            standard star spectrum. The structure contains the
%            following fields (see wavecalib_1d_int.m):
%            .WaveCalibOut
%            .RMS
%            .ParWaveCalib
% Tested : Matlab 7.0
%     By : Eran O. Ofek                   January 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------------
DefArcImage   = [];
DefArcType    = [];
if (nargin==2),
   ArcImage    = DefArcImage;
   ArcType     = DefArcType;
elseif (nargin==3),
   ArcType     = DefArcType;
else
   % do nothing
end

% set default values:
ExpKey     = 'EXPOSURE';
DispAxis   = 'x';
AutoSel    = 'n';
Int        = 'y';
InterpMethod = 'linear';
KeyAirMass = 'AIRMASS';
Ext        = 'VLT_atmospheric_extinction.dat';

% get keywords value:
Narg = length(varargin);
for Iarg=1:2:Narg-1,
   switch lower(varargin{Iarg})
    case 'expkey'
       ExpKey         = varargin{Iarg+1};
    case 'dispaxis'
       DispAxis       = varargin{Iarg+1};
    case 'autosel'
       AutoSel        = varargin{Iarg+1};
    case 'interpmethod'
       InterpMethod   = varargin{Iarg+1};
    case 'keyairmass'
       KeyAirMass     = varargin{Iarg+1};
    case 'ext'
       Ext            = varargin{Iarg+1};
    otherwise
       error('Unknown keyword option');
   end
end


disp('ds9 should be open')

% trace spectrum
switch lower(AutoSel)
 case 'y'
    error('AutoSel=y not availble yet');
 case 'n'
    [Std.Trace,Std.ExtractedSpec]=spec_trace(StdImage,[],'DispAxis',DispAxis,'Int',Int,'PeakMethod','wmean','Deg',[]);
    close;
 otherwise
    error('Unknwon AutoSel option');
end

% subtract sky from StD spectrum
[Std.SkySubSpec,Std.SkyVal,Std.SkyRMS]=spec_skysub(Std.ExtractedSpec,[],'SkyMet','poly1');
close;

% extract StD spectrum
[Std.Fit] = spec_extract1(Std.SkySubSpec,Std.SkyVal,'Method','Num');

%[Std.Fit] = spec_extract1(Std.SkySubSpec,Std.SkyVal,'Method','Aper','Aper',[54 59]);


% wavelength calibration 
SizeES  = size(Std.ExtractedSpec,1);
SpecPos = round((SizeES+1).*0.5); 

if (isempty(ArcImage)==1),
   % attempt to use sky lines for wavelength calibration
   [SelWaveCalib,WaveCalibOut,RMS,ParWaveCalib]=wavecalib_1d_int(Std.SkyVal(SpecPos,:).');

else
   % use arc for wavelength calibration
   TraceMat = [Std.Trace.Disp, Std.Trace.SpatP];
   ExtSemiW = 50;
   ExtractedArc = extract_trace(ArcImage,TraceMat,ExtSemiW,InterpMethod,DispAxis);

   SizeEA  = size(ExtractedArc,1);
   ArcPos  = round((SizeEA+1).*0.5); 
   [SelWaveCalib,WaveCalibOut,RMS,ParWaveCalib]=wavecalib_1d_int(ExtractedArc(ArcPos,:).');
end
close;
close;

disp(sprintf('Wavelength calibration RMS = %f',RMS))

Std.ObservedSpec = [WaveCalibOut(:,2), Std.Fit.Sum];

% get StD spectrum
[SpecSTD] = get_specstand(StdName);

% sample StD spectrum as Std.ObservedSpec
[Std.ActualSpec,Junk] = eq_sampling(SpecSTD,Std.ObservedSpec,Std.ObservedSpec(:,1));

% get airmass of StD spectrum
if (isstr(KeyAirMass)==1),
   KeyVal  = get_fits_keyword(StdImage,{KeyAirMass});
   AirMass = KeyVal{1};
else
   AirMass = KeyAirMass;
end


if (isstr(AirMass)==1),
   AirMass = str2num(AirMass);
end
Std.ExtCorrSpec = atmospheric_ext(Std.ObservedSpec,AirMass,Ext,InterpMethod);

PrelimCal = [Std.ExtCorrSpec(:,1),Std.ActualSpec(:,2)./Std.ExtCorrSpec(:,2)];
I = find(isnan(PrelimCal(:,2))==0);
PrelimCal = PrelimCal(I,:);
PlotOption = 1;
SigClip    = [-Inf Inf];
MaxIter    = 1;
[CleanPrelimCal,J1,J2,Par,Res,I1,I2,GH] = plot_rm(PrelimCal(:,1),PrelimCal(:,2),'b.','polyfit_sc',5,[],SigClip,MaxIter,PlotOption);

FluxCalibTran = polyval(Par,Std.ObservedSpec(:,1));


WaveCalib.WaveCalibOut  = WaveCalibOut;
WaveCalib.RMS           = RMS;
WaveCalib.ParWaveCalib  = ParWaveCalib;
