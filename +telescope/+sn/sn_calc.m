function [SN,OutPars,StProp]=sn_calc(varargin)
% A signal-to-noise calculator for astronomical telescopes.
% Package: telescope.sn
% Description: A signal-to-noise calculator for astronomical telescopes.
%              Calculate S/N or limiting magnitude and field of view
%              properties.
% Input  : * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'Mag'   - Magnitude column vector at which to calculate S/N.
%                      If the 'SN' field is provided, then this input
%                      is ignored.
%            'SN'    - Column vector of S/N at which to calculate
%                      limiting magnitude.
%                      Default is empty []. If empty, then calculate the
%                      S/N at the magnitude given in 'Mag'.
%            'AM'    - Airmass of observation. If empty than assume
%                      there is no atmospheric extinction. Default is 1.3.
%            'AtmExtFile' - Atmospheric extinction file.
%                      See atmospheric_ext.m for details.
%                      Default is 'KPNO_atmospheric_extinction.dat'.
%            'ScatL' - Fraction of scattered light from the entire
%                      star light in the field of view. Default is 0.0.
%            'Tel'   - Name of telescope system with predefined
%                      parameters {'ULTRASAT','PTF','ZTF','Wise100cm'}.
%                      If empty than don't use telescope.
%                      NOT YET AVAILABLE.
%            'Aper'  - Telescope diameter [cm]. Default is 60 cm.
%            'FocalLength' - Telescope effective focal length [cm].
%                      Default is 95 cm.
%            'Family'- Observation band family (e.g., 'SDSS').
%                      See get_astfilter.m for options.
%            'Band'  - Observation band name (e.g., 'r').
%                      See get_astfilter.m for options.
%            'MagType'- Magnitude type {'AB','Vega'}. Default is 'AB'.
%            'QOE'   - Total quantum and optical efficiency of the
%                      telescope including obscurations. Default is 0.5.
%            'FWHM'  - Effective FWHM [arcsec]. Default is 2.5.
%            'Bin'   - Binning size. Default is 1.
%            'PixSize'- Pixel size [microns]. Default is 6.5.
%            'SizeCCD'- CCD size in pixels. Default is [2048 2048].
%            'Nccd'  - Number of CCDs. Default is 1.
%            'ReadTime'- CCD readout time [s]. Default is 1.
%            'RN'    - Readout noise [e-]. Default is 2.
%            'DC'    - Dark current [e-/s/pix]. Default is 0.016.
%            'WD'    - CCD well depth. Default is 1e5.
%            'Gain'  - Detector gain [e-/ADU]. Default is 1.5.
%            'CTE'   - Detector CTE. Default is 0.99999.
%            'PosCCD'- Position of source on CCD to calculate CTE losses.
%                      [pix pix]. Default is [1024 1024].
%            'ExpTime'- Exposure time [s]. Default is 1.
%            'CRsplit'- Number of images. Default is 1.
%            'Back'  - Background [mag/arcsec^2]. Default is 20.3.
%            'BackUnits' - Background units. Default is 'mag/arcsec^2'.
%            'PhotRad'- Row vector of photometric aperture radius, to use
%                      in the search for the optimal photometry radius.
%                      Defaults is (1.5:0.24).
%            'MagVecSN'- Vector of magnitudes in which to look for
%                      limiting magnitude given S/N.
%                      Default is (5:0.5:30).'.
%            'Temp'  - Effective black body temperature of source [K].
%                      Default is 5770 K.
%            'InterpMethod' - Interpolation method. Default is 'cubic'.
%                      See interp1.m for options.
%            'RA'    - J2000.0 R.A. in which to calculate the properties
%                      of stars brighter than the magnitude.
%                      Default is '10:00:00.0'.
%                      See convertdms.m for options.
%            'Dec'   - J2000.0 Dec. Default is '+20:00:00.0'.
%                      See convertdms.m for options.
%            'FOV'   - Field of view in arcsec in which to query the USNO-B
%                      catalog. Default is 1000.
% Output : - Structure containing the S/N properties. The following fields
%            are available:
%            .Mag    - Magnitude or limiting mag.
%            .SN     - S/N.
%            .VarRatio - Contribution of the various noise to the total
%                      variance.
%            .BestPhotRad - Optimal photometric aperture radius [arcsec].
%            .FracInCentralPix - Fraction of flux in central pixel.
%            .SaturationMag - Saturation magnitude.
%            .SNrad - S/N as a function of photometric radius.
%            .PhotRad - photometric radius [arcsec].
%          - Structure containing additional parameters.
%          - Structure containing the stellar density in the field.
%            This is calculated using the USNO-B1 R2 mag, and the colors
%            are B2-R2.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [SN,OP,S]=telescope.sn.sn_calc('ExpTime',1./30,'SN',10); SN.Mag, OP.PixScale, OP.totFOVarea, OP.totFOVarea.*S.StarDensity
%          [SN,OP]=telescope.sn.sn_calc('Mag',[16;17]);
%          [SN,OP]=telescope.sn.sn_calc('SN',[5;10]);
% Reliable: 2
%--------------------------------------------------------------------------


RAD      = 180./pi;

%--- Main input ---
DefV.Mag          = 0;
DefV.SN           = [];

%--- Telescope ----
%DefV.Tel          = [];
DefV.Aper         = 27.9; %3.5./1.4; %57; %48*2.54; %30; %55; %1.8; %55;        % [cm] Telescope aperture
DefV.FocalLength  = 62; %3.5; %108; %48.*2.54.*2.8; %30.*1.44; %108; %1.7; %108;        % [cm] Telescope focal length
DefV.Family       = 'SDSS'; %'Cousins'; %'Cousins'; %'Cousins'; %'Johnson'; %'SDSS';  %
DefV.Band         = 'g'; %'R'; %'R'; %'R'; %'r';        %
DefV.MagType      = 'AB'; %'Vega'; %'AB';       % {'AB','Vega'}
%--- Optics ---
DefV.QOE          = 0.75; %0.3 %0.45; %[0.6 0.75 0.85];        % []   quantum-optical effiency [including CCD, optics and filter]
DefV.FWHM         = 2.8; % 2.5; %2.5; %3.7;          % ["]  Point source image FWHM (FWHM=2.3548\sigma) [including pixel size, optics, diffraction, jitter, charge leakage]
%--- CCD ---
DefV.Bin          = 1;           % [pix]     pixel binning
DefV.PixSize      = 3.6; %11.5; %6.5; %6.5;         % [micron]  CCD pixel size
DefV.SizeCCD      = [6000 6000]; %[2048 2500]; % [pix pix] size of a single CCD
DefV.Nccd         = 7;           % []        number of CCDs or telescopes.
%DefV.ReadRate     = 2e5;         % [pix/s]   CCD readout rate
DefV.ReadTime     = 0;           % [s]       CCD readout time - if given, then used instead of ReadRate
DefV.RN           = 1.5; %1.5; %1.5; %1.5; %1.5; %4.*sqrt(DefV.ReadRate./250);;           % [e-]      CCD Read noise
DefV.DC           = 0.1; %16;       % [e-/s/pix]CCD dark current
DefV.WD           = 1e5;         % [e-]      CCD wall depth/saturation limit per pixel
DefV.Gain         = 1.5;         % [e-/DN]   CCD Gain
DefV.CTE          = 0.99999;     % []        CCD CTE
DefV.PosCCD       = [1024 1024]; % [pix pix] position of source on CCD
%--- Exposure ---
DefV.ExpTime      = 15; %30; %1./30;           % expoure time per a single CRsplit
%DefV.CycleExpTime = 1;          % [s]       Exposure time including overheads and CRsplit
DefV.CRsplit      = 1; %1;           % []        Number of images in Exposure time
%     ReadoutTime = SizeCCD(1)./ReadRate./Namp/(Bin.^2)
%     DeadTime    = ReadoutTime.*CRsplit
%     ExpTime     = CycleExpTime - DeadTime;
%     ReadNoise   = RN.*sqrt(CRsplit);      note that in practice this maybe higher

%--- additional noise ---
DefV.ScintNoise   = 0; %0.1;
DefV.ScintFreq    = 30;

%--- Background ---
DefV.Back         = 21.0; %19.9; %15.6; %19.9; %26.78; %19.9; %20.3;  % 21.3 % mag/arcsec^2 in filter
DefV.BackUnits    = 'mag/arcsec^2';  % {'mag/arcsec^2','photons/cm^2/s/arcsec^2'}
DefV.ScatL        = 0.0;  %0.01;

%--- Atmosphere ---
DefV.AM           = 1.3;
DefV.AtmExtFile   = 'KPNO_atmospheric_extinction.dat';  % atmospheric extinction file
%--- Photometry ---
DefV.PhotRad      = (0.2:0.2:100);  % ["] radius of photometry aperture
%DefV.Annulus      = [20 50];      % ["] Inner and outer Background annulus radius
%DefV.FFerr        = 0.005;        % Flat Field relative error.

%--- Flux input ---
DefV.MagVecSN     = (0:0.5:30).';
%DefV.InputFlux    = 'radtemp';     % {[radtemp] | bollumtemp}
%DefV.Radius       = 696000e5.*500; % [cm]    source radius
DefV.Temp         = 1e5; %5770;           % [K]     source temperature
%DefV.BolLum       = 1e43;           % [erg/s] source bolumetric luminosity
%DefV.Dist         = 61e6;          % [pc]    source luminosity distance


%--- general ---
DefV.InterpMethod = 'pchip';

%--- stellar properties ---
DefV.RA  = '18:51:00.0';
DefV.Dec = '-06:00:00.0';

%DefV.RA  = '20:23:12';
%DefV.Dec = '+40:46:42';

DefV.FOV = 1000;  % [arcsec]

% hidden pars
DefV.ScatLback = 0;

%InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

InPar.QOE = prod(InPar.QOE);

%if (~isempty(InPar.Tel)),
%   switch lower(InPar.Tel)
%       case 'ultrasat'

if (InPar.ScatL>0)
    % take into account scattered light
    ScatLprop = stellar_prop(InPar,25);
    
    
    [SN,OutPars]=signal2noise_calc(InPar);
    SumMag = sum_mag(ScatLprop.Mag);
    
    %.*OutPars.EffArea.*InPar.ExpTime
    % scattered light background [counts/arcsec^2/s/cm^2 in filter]
    InPar.ScatLback = InPar.ScatL.*OutPars.Counts0.*10.^(-0.4.*SumMag)./(ScatLprop.FOVareaDeg2.*3600.^2);

else
    InPar.ScatLback = 0;
end
 

if (isempty(InPar.SN))
    % calculate S/N for given magnitude
    [SN,OutPars]=signal2noise_calc(InPar);
else
    % calculate limiting magnitude for S/N
    InPar.Mag = InPar.MagVecSN;
    [SN,OutPars] = signal2noise_calc(InPar);
    MagSN = interp1(SN.SN,InPar.MagVecSN,InPar.SN,InPar.InterpMethod);
    InPar.SN = [];
    InPar.Mag = MagSN;
    [SN,OutPars] = signal2noise_calc(InPar);
    
end

if (nargout>2)
    StProp = stellar_prop(InPar,SN.Mag);
end
    






function [SN,OutPars]=signal2noise_calc(InPar)
RAD     = 180./pi;
SigmaB  = constant.sigma;
h       = constant.h;
SolR    = constant.SunR;

% Telescope/focal plane parameters
OutPars.PixScale    = InPar.Bin.*InPar.PixSize.*1e-4./InPar.FocalLength .*RAD.*3600;  % [arcsec/pix]
OutPars.ccdFOVdim   = InPar.SizeCCD.*OutPars.PixScale;      % [arcsec]
OutPars.ccdFOVarea  = prod(OutPars.ccdFOVdim)./(3600.^2);   % [deg^2]
OutPars.totFOVarea  = OutPars.ccdFOVarea.*InPar.Nccd;       % [deg^2]

% Exposure time
OutPars.ExpTime     = InPar.ExpTime;
OutPars.totExpTime  = (InPar.ExpTime+InPar.ReadTime).*InPar.CRsplit;

% get band
if (ischar(InPar.Family))
    Filter = AstFilter.get(InPar.Family,InPar.Band);
    % normalize band maximum to 1
else
    %Filter.nT{1} = InPar.Family;
    Filter.nT = InPar.Family;
end
% normalize band maximum to 1
Filter.nT(:,2) = Filter.nT(:,2)./max(Filter.nT(:,2));
%Filter.nT{1} = [4999 0; 5000 1; 6000 1; 6800 1; 6801 0];

%------------------------------
%--- Atmospheric extinction ---
%------------------------------
if (~isempty(InPar.AM))
    [~,ExtFactor]=ImUtil.Spec.atmospheric_ext(Filter.nT,InPar.AM);
    % calculate the weighted mean of extinction factor within filter
    % assuming flat f_lambda spectrum
    OutPars.ExtinFactor = sum(ExtFactor.*Filter.nT(:,2))./sum(Filter.nT(:,2));
else
    OutPars.ExtinFactor = 1;
end

%------------
% Input flux
%------------
if (isempty(InPar.Mag))
%    % Input flux is specified by BolLum and Temp
%    Temp   = InPar.Temp;
%    Dist   = InPar.Dist;
%    Radius = sqrt(InPar.BolLum./(4.*pi.*SigmaB.*Temp.^4));
% 
%    Ner = length(Radius);
%    Net = length(Temp);
%    Ned = length(Dist);
%    MNe = max([Ner;Net;Ned]);
% 
%    if (MNe>Ner),
%       Radius = Radius.*ones(MNe,1);
%    end
%    if (MNe>Net),
%        Temp   = Temp.*ones(MNe,1);
%    end
%    if (MNe>Ned),
%        Dist   = Dist.*ones(MNe,1);
%    end
% 
%    Flux   = zeros(MNe,1);
%    Counts = zeros(MNe,1);
%    for Imn=1:1:MNe,
%       [Flux(Imn),Counts(Imn)] = spec_photon_counts(Temp(Imn),Filter.nT{1},[],Radius(Imn),Dist(Imn));
%    end
% 
%    Mag = blackbody_mag_c(Temp,InPar.Family,InPar.Band,InPar.MagType,Radius,Dist);
%  
else
   % Input flux is specified by Magnitude and eff. Temperature
   % convert Magnitude to counts
   Net      = length(InPar.Temp);
   Mag      = InPar.Mag;

   Flux     = zeros(Net,1);
   Counts   = zeros(Net,1);
   %DeltaMag = zeros(Net,1);
   for Imn=1:1:Net
      [Flux(Imn),Counts(Imn)] = telescope.sn.spec2photons(InPar.Temp(Imn),Filter.nT,[],SolR,10);
   end
   
   % for Mag 0 star
   DeltaMag0  = AstroUtil.spec.blackbody_mag_c(InPar.Temp,InPar.Family,InPar.Band,InPar.MagType,SolR,10);
   Flux0      = Flux.*10.^(+0.4.*DeltaMag0);
   Counts0    = Counts.*10.^(+0.4.*DeltaMag0);
   
   
   DeltaMag  = AstroUtil.spec.blackbody_mag_c(InPar.Temp,InPar.Family,InPar.Band,InPar.MagType,SolR,10) - Mag;
   Flux      = Flux.*10.^(+0.4.*DeltaMag);
   Counts    = Counts.*10.^(+0.4.*DeltaMag);

   % correct Flux and Counts to atmospheric extinction
   Flux      = Flux./OutPars.ExtinFactor;
   Counts    = Counts./OutPars.ExtinFactor;
   
   %Mag1 = blackbody_mag_c(InPar.Temp,InPar.Family,InPar.Band,InPar.MagType,1,10);

end

OutPars.Flux0   = Flux0;      % flux of 0 mag star (in band) per cm^2 per second
OutPars.Counts0 = Counts0;    % number of photons of 0 mag star (in band) per cm^2 per second


% System aperture and total throuput
OutPars.Area    = pi.*(InPar.Aper./2).^2;   % Telescope aperturea area [cm^2]
OutPars.EffArea = OutPars.Area.*InPar.QOE;

% Calculate background
switch lower(InPar.BackUnits)
    case 'mag/arcsec^2'
        OutPars.BckCountsAS  = Counts0.*10.^(-0.4.*InPar.Back);    % backgrounds counts/s/cm^2/arcsec^2 in filter
        OutPars.BckCountsAS  = OutPars.BckCountsAS.*OutPars.EffArea; % backgrounds counts/s/arcsec^2 in filter and telescope
        OutPars.BckCountsPix = OutPars.BckCountsAS.*OutPars.PixScale.^2;  % backgrounds counts/s/pix in filter and telescope

    case 'photons/cm^2/s/arcsec^2'
        OutPars.BckCountsAS  = InPar.Back;    % backgrounds counts/s/cm^2/arcsec^2 in filter
        OutPars.BckCountsAS  = OutPars.BckCountsAS.*OutPars.EffArea; % backgrounds counts/s/arcsec^2 in filter and telescope
        OutPars.BckCountsPix = OutPars.BckCountsAS.*OutPars.PixScale.^2;  % backgrounds counts/s/pix in filter and telescope

    otherwise
        error('Unknown BackUnits option');
end


% scattered light background
OutPars.ScBckCountsAS  = InPar.ScatLback;   % scat l. backgrounds counts/s/cm^2/arcsec^2 in filter
OutPars.ScBckCountsAS  = OutPars.ScBckCountsAS.*OutPars.EffArea; % scat l. backgrounds counts/s/arcsec^2 in filter and telescope
OutPars.ScBckCountsPix = OutPars.ScBckCountsAS.*OutPars.PixScale.^2;  % scat l. backgrounds counts/s/pix in filter and telescope


% total counts in background annulus
% pi.*(max(InPar.Annulus).^2 - min(InPar.Annulus).^2).*OutPars.BckCountsPix 

% fraction of light within photometric aperture
PSF_Sigma = InPar.FWHM./(2.*sqrt(2.*log(2)));  % [arcsec]
FracLight = 1 - exp(-InPar.PhotRad.^2./(2.*PSF_Sigma.^2));
NumberPixPhotRad = ceil(pi.*(InPar.PhotRad./OutPars.PixScale).^2);

% Calculate S/N:
%Signal = Counts.*OutPars.EffArea.*FracLight.*InPar.ExpTime;
CTE_Factor = InPar.CTE.^(sum(InPar.PosCCD));
Signal = bsxfun(@times,Counts,FracLight).*OutPars.EffArea.*InPar.ExpTime.*InPar.CRsplit;
Var.Source = Signal.*CTE_Factor;
Var.Bck    = OutPars.BckCountsPix.*NumberPixPhotRad.*InPar.ExpTime.*InPar.CRsplit;
Var.ScBck  = OutPars.ScBckCountsPix.*NumberPixPhotRad.*InPar.ExpTime.*InPar.CRsplit;
Var.DC     = InPar.DC.*NumberPixPhotRad.*InPar.ExpTime.*InPar.CRsplit;
Var.RN     = NumberPixPhotRad.*InPar.RN.^2.*InPar.CRsplit;
Var.Digi   = NumberPixPhotRad.*(0.289.*InPar.Gain).^2.*InPar.CRsplit;

% scintilation noise (relayive)
ScintRelNoise  = InPar.ScintNoise./sqrt(InPar.ExpTime.*InPar.ScintFreq);

Noise2 = bsxfun(@plus, Var.Source + Var.Source.*ScintRelNoise, ...
                       Var.Bck + Var.ScBck + Var.DC + Var.RN + Var.Digi);
                       
SN.Mag    = Mag;


SN.VarRatio.Source = bsxfun(@rdivide,Var.Source,Noise2);
SN.VarRatio.Bck    = bsxfun(@rdivide,Var.Bck,Noise2);
SN.VarRatio.ScBck  = bsxfun(@rdivide,Var.ScBck,Noise2);
SN.VarRatio.DC     = bsxfun(@rdivide,Var.DC,Noise2);
SN.VarRatio.RN     = bsxfun(@rdivide,Var.RN,Noise2);
SN.VarRatio.Digi   = bsxfun(@rdivide,Var.Digi,Noise2);

SNrad = Var.Source.*CTE_Factor./sqrt(Noise2);

[SN.SN,MaxInd] = max(SNrad,[],2);
SN.BestPhotRad = InPar.PhotRad(MaxInd);
SN.FracInCentralPix = 1 - exp(-(OutPars.PixScale.*0.6).^2./(2.*PSF_Sigma.^2));  % fraction of light in central pixel
SN.SaturationMag = SN.Mag - 2.5.*log10(InPar.WD./(Signal(sub2ind(size(Signal),(1:length(InPar.Mag)).',MaxInd)).*SN.FracInCentralPix));



SN.SNrad = SNrad;
SN.PhotRad = InPar.PhotRad;
SN.Signal    = Signal(sub2ind(size(SN.VarRatio.Source),(1:length(InPar.Mag)).',MaxInd));
SN.Var       = Var;
SN.Noise2    = Noise2;

SN.VarRatio.Source = SN.VarRatio.Source(sub2ind(size(SN.VarRatio.Source),(1:length(InPar.Mag)).',MaxInd));
SN.VarRatio.Bck    = SN.VarRatio.Bck(sub2ind(size(SN.VarRatio.Bck),(1:length(InPar.Mag)).',MaxInd));
SN.VarRatio.ScBck  = SN.VarRatio.ScBck(sub2ind(size(SN.VarRatio.ScBck),(1:length(InPar.Mag)).',MaxInd));
SN.VarRatio.DC     = SN.VarRatio.DC(sub2ind(size(SN.VarRatio.DC),(1:length(InPar.Mag)).',MaxInd));
SN.VarRatio.RN     = SN.VarRatio.RN(sub2ind(size(SN.VarRatio.RN),(1:length(InPar.Mag)).',MaxInd));
SN.VarRatio.Digi   = SN.VarRatio.Digi(sub2ind(size(SN.VarRatio.Digi),(1:length(InPar.Mag)).',MaxInd));

SN.NpixBestPhotAper = NumberPixPhotRad(MaxInd);

function StProp=stellar_prop(InPar,LimMag)

InvRAD = pi./180;

RA  = convertdms(InPar.RA,'gH','r');
Dec = convertdms(InPar.Dec,'gD','R');

[Cat,ColCell,Col,UnitsCell] = wget_usnob1(RA,Dec,InPar.FOV.*InvRAD./3600);
StProp.FOVareaDeg2 = pi.*(InPar.FOV./3600).^2;

Nlm = length(LimMag);
for Ilm=1:1:Nlm
   Ist = find(Cat(:,Col.R2mag)<LimMag(Ilm));
   StProp.Nstar = length(Ist);
   StProp.StarDensity = StProp.Nstar./(pi.*(InPar.FOV./3600).^2);
   StProp.MedianMag = nanmedian(Cat(Ist,Col.R2mag));
   StProp.Mag   = Cat(Ist,Col.R2mag);
   StProp.Color = Cat(Ist,Col.B2mag) - Cat(Ist,Col.R2mag);
   
   % estimate stars angular radius
   Par = [-0.05084      0.29303     -0.70919       4.1467];
   StProp.Color(StProp.Color<-0.17) = -0.17;
   Tstar = 10.^(polyval(Par,StProp.Color));
   StProp.Theta = 0.47.*sqrt((Tstar./5700).^(-4).*10.^(-0.4.*(StProp.Mag-4.7)));
   
   %[StProp.Theta,StProp.FlagTheta]=star_ang_rad(StProp.Mag,StProp.Color.*1.188+0.15,'g','g-r');
   FT = StProp.Theta<0.005./2; %0.044./2;
   StProp.FlaggedStarDensity = sum(FT)./StProp.FOVareaDeg2;  % per deg^2
   
   
end




