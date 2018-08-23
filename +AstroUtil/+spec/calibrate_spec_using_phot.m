function [ObsRestSpec,CorrSpec,Par]=calibrate_spec_using_phot(Spec,SpecDate,LC,FiltFam,FiltName,FiltSys,varargin)
%--------------------------------------------------------------------------
% calibrate_spec_using_phot function                             AstroSpec
% Description: Given a spectrum and a light curve in one or more bands,
%              applay synthetic photometry to the spectrum to calibrate
%              it against the light curves. Optionally, if RA, Dec are
%              provided, dereddned the spectrum for Galactic extinction.
%              Furthermore, if redshift is provided then the spectrum
%              will be converted to the rest frame, and if additional
%              host extinction is given then another deredenning at the
%              rest frame wavelength will be applied.
%              If more than one filter is provided than the spectrum
%              calibration will be performed using a polynomial fitting
%              when the deg of the polynom is N-1, where N is the good
%              filters.
% Input  : - Spectrum [Wave(Ang), Flux(per unit wavelength)]
%            Alternatively this can be a file name containing the
%            spectrum name.
%          - Date [Day Month Year Hour Min Sec] or JD at which the
%            spectrum was obtained.
%          - Calibrated light curve of the source [JD Mag].
%            Optionally if a cell array then each cell contains the 
%            light curve in a different band.
%            If a string or a cell array of strings is provided then
%            the program will attempt to use the read_mark_lc.m program
%            to read the file into a light curve.
%          - String of filter family (see AstFilter.get.m for options).
%            e.g., 'PTF', 'SDSS'.
%            If the light curve input argument is a cell array then
%            this a cell array of filter families should be provided.
%            Alternatively, this can be an AstFilter class object.
%          - String or cell array of filter names (e.g., 'R').
%          - String or cell array of the magnitude system for the
%            light curve (i.e., 'AB','Vega').
%          * Arbitrrary number of pairs of input arguments: ...,key,val,...
%            The following keywords are available:
%            'RA'    - J2000.0 Right Ascension [radians] of the source.
%                      Default is empty matrix.
%                      If empty matrix then Galactic extinction is not
%                      applied, unless the 'Ebv' keyword is provided.
%            'Dec'   - J2000.0 Declination [radians] of the source.
%                      Default is empty matrix.
%                      If empty matrix then Galactic extinction is not
%                      applied, unless the 'Ebv' keyword is provided.
%            'Ebv'   - E_{B-V} [mag] at the source position.
%                      Default is 0, unless 'RA' and 'Dec' are provided.
%            'z'     - Redshift of the source. Default is 0.
%            'Ebv_z' - Additional E_{B-V} extinction at the redshift of
%                      the source to applay to the rest frame spectrum.
%                      Default is 0.
%            'Rv'    - Galactic R_{V}. Default is 3.08.
%            'Rv_z'  - Host galaxy R_{V}. Default is 3.08.
%            'MaxFlag'- maximum fraction of filter outside spectrum range.
%                      Default is 0.05.
%            'PolyDeg'- Degree of polynomial to fit to the syn. phot.
%                      minus interpolated phot. Default is number of
%                      good filters - 1.
%            'InterpMethod' - Light curve interpolation method.
%                      See interp1.m for options. Default is 'linear'.
% Output : - Output calibrated observed frame spectrum
%            [Wave(Ang), Flux(erg/cm^2/s/Ang)]
%            i.e., only syn. phot. calibration is applied.
%          - Output calibrated spectrum with all requested corrections
%            applied (rest frame).
%          - Structure of various parameters.
% Tested : Matlab R2013A
%     By : Eran O. Ofek                   July 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [ObsRestSpec,CorrSpec,Par]=calibrate_spec_using_phot(Spec,...
%                                     [5 6 2013 0.5],LC,'PTF','R','AB');
%          [ObsRestSpec,CorrSpec,Par]=calibrate_spec_using_phot(Spec,...
%                                     2461545.3,'PTF10cwx.out_PTF48R','PTF','R','AB');
%          [ObsRestSpec,CorrSpec,Par]=calibrate_spec_using_phot(Spec,...
%                                     2461545.3,'PTF10cwx.out_PTF48R','PTF','R','AB',...
%                                     'z',0.8,'RA',1,'Dec',0.2,'Rv_z',2.6);
% Reliable: 2
%--------------------------------------------------------------------------

DefV.RA      = [];
DefV.Dec     = [];
DefV.Ebv     = 0;
DefV.z       = 0;
DefV.Ebv_z   = 0;
DefV.Rv      = 3.08;  % R_{V} of Galactic extinction
DefV.Rv_z    = 3.08;  % R_{V} of host extinction
DefV.MaxFlag = 0.05;  % maximum fraction of filter outside spectrum range.
DefV.PolyDeg = -1;
DefV.InterpMethod = 'linear';

%InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% spectrum date into JD
if (length(SpecDate)>1)
    SpecJD = convert.date2jd(SpecDate);
else
    SpecJD = SpecDate;
end

% read spectrum
if (ischar(Spec))
    SpecName = Spec;
    clear Spec;
    Spec = Util.IO.load2(SpecName);
end

% move LC and filters into cell array
if (~iscell(LC))
    LC = {LC};
end
if (AstFilter.isAstFilter(FiltFam))
    FiltName = {FiltFam.band};
    FiltFam  = {FiltFam.family};
end

if (~iscell(FiltFam))
    FiltFam = {FiltFam};
end
if (~iscell(FiltName))
    FiltName = {FiltName};
end
if (~iscell(FiltSys))
    FiltSys = {FiltSys};
end

% if only one FiltSys is provided then duplicate for all filters
if (iscell(FiltSys))
    if (length(FiltSys)==1 && length(FiltName)>1)
        [FiltSys{1:1:length(FiltName)}]=deal(FiltSys{1});
    end
end

% for eachg LC
Nlc         = length(LC);
SynPhot     = zeros(Nlc,1);
FlagSynPhot = zeros(Nlc,1);
FilterWave  = zeros(Nlc,1);
InterpLC    = zeros(Nlc,1);
for Ilc=1:1:Nlc
    % read LC
    if (ischar(LC{Ilc}))
        % attempt to read LC using read_mark_lc
        LCdata{Ilc} = VO.PTF.read_mark_lc(LC{Ilc});
        LCdata{Ilc}(:,1) = LCdata{Ilc}(:,1) + 2400000.5;  % convert MJD to JD
    else
        LCdata{Ilc} = LC{Ilc};
    end
    % interpolate LCs to JD
    InterpLC(Ilc) = interp1(LCdata{Ilc}(:,1),LCdata{Ilc}(:,2),SpecJD,InPar.InterpMethod);
    
    % run syn phot. on spectrum
    [SynPhot(Ilc), FlagSynPhot(Ilc),FilterWave(Ilc)] = AstroUtil.spec.synphot(Spec,FiltFam{Ilc},FiltName{Ilc},FiltSys{Ilc});
end
%[SynPhot, FlagSynPhot,FilterWave, InterpLC];

Igood = find(FlagSynPhot<InPar.MaxFlag & ~isnan(SynPhot));
Ngood = length(Igood);
if (Ngood==0)
    error('Filters doesnot cover spectrum completly or spectrum JD is not covered by LC');
end
if (InPar.PolyDeg<0)
   InPar.PolyDeg = Ngood - 1;
end

PolyPar = polyfit(FilterWave, InterpLC(Igood) - SynPhot(Igood), InPar.PolyDeg);
FluxCorr = 10.^(-0.4.*polyval(PolyPar,Spec(:,1)));
ObsRestSpec = Spec;
ObsRestSpec(:,2) = ObsRestSpec(:,2).*FluxCorr;

% Get extinction from Schlegal maps
if (~isempty(InPar.RA) && ~isempty(InPar.Dec))
    InPar.Ebv = sky_ebv(InPar.RA,InPar.Dec);
else
    % use default Ebv
end

GalExtinCorr = 10.^(+0.4.*(AstroUtil.spec.extinction(InPar.Ebv,ObsRestSpec(:,1)./10000,[],InPar.Rv)));
CorrSpec = ObsRestSpec;
CorrSpec(:,2) = CorrSpec(:,2).*GalExtinCorr;

% cosmology
if (InPar.z~=0)
   CorrSpec = shift_spec(CorrSpec,InPar.z,'w','f');
end

% rest frame (e.g., host galaxy) extinction
RestExtinCorr = 10.^(+0.4.*(AstroUtil.spec.extinction(InPar.Ebv_z,CorrSpec(:,1)./10000,[],InPar.Rv_z)));
CorrSpec(:,2) = CorrSpec(:,2).*RestExtinCorr;

Par.PolyDeg       = InPar.PolyDeg;
Par.GalExtinCorr  = GalExtinCorr;
Par.RestExtinCorr = RestExtinCorr;
Par.SynPhot       = SynPhot;
Par.FlagSynPhot   = FlagSynPhot;
Par.FilterWave    = FilterWave;
Par.Ebv           = InPar.Ebv;