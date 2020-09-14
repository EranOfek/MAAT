function [Spec]=zodiac_bck(RA,Dec,Date,varargin)
% Calculate the zodi calibrated spectrum and magnitude at a given position.
% Package: +ultrasat
% Description: Calculate the zodiacal calibrated spectrum and magnitude
%              at a given sky position.
% Input  : - Vector of J2000.0 R.A., or longitude [radians].
%          - Vector of J2000.0 Dec, or latitude [radians].
%          - Date [D M Y F], or JD (column matrix). Default is empty.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Wave' - Vector of wavelength [Ang] in which to calculate
%                   spectrum. If empty, then use original data.
%                   Default is empty.
%            'OutType' - 'mat'|'astspec'. Default is 'mat'.
%            'CooSys'  - Input coordinate system. Default is 'j2000.0'.
%                        Options include 'hec' - heliocentric ecliptic,
%                        'e','g', etc.
%            'FilterFamily' - Default is 'SDSS'.
%            'FilterName'   - Default is 'g'.
%            'FilterSys'  - Default is 'AB'.
%            'InterpMethodSpec' - Default is 'linear'.
%            'InterpMethodCoo'  - Default is 'linear'.
% Output : - A structure with the following fields:
%            'RA'  - Vector of RA/long corresponding to spectra.
%            'Dec' - Vector of Dec/latitude.
%            'JD'  - Vectoe of JD.
%            'Wave'- Vector of wavelength.
%            'Spec'- A matrix in which each column is the zodi spectrum
%                    in units of erg/cm^2/s/Ang/arcsec^2.
%                    Or a vector of AstSpec objects with the spectra.
%            'Mag' - Vector of magnitude in requested band.
%            'MagV'- Vector of Vega V-band magnitudes.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Spec]=ultrasat.zodiac_bck(1,1,2451545);
 %         [Spec]=ultrasat.zodiac_bck(1,1,2451545,'CooSys','g')
 %         [Spec]=ultrasat.zodiac_bck(1,1,[],'CooSys','hec')
% Reliable: 2
%--------------------------------------------------------------------------



DefV.Wave                 = [];
DefV.OutType              = 'mat';
DefV.CooSys               = 'j2000.0';  % 'g','ec'
DefV.FilterFamily         = 'SDSS';
DefV.FilterName           = 'g';
DefV.FilterSys            = 'AB';

DefV.InterpMethodSpec     = 'linear'; 
DefV.InterpMethodCoo      = 'linear';
%
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

RA  = RA(:);
Dec = Dec(:);

Ncol = max(numel(RA),numel(Dec));
RA   = RA.*ones(Ncol,1);
Dec  = Dec.*ones(Ncol,1);




switch lower(InPar.CooSys)
    case 'hec'
        % Helio-ecliptic coordinates
        % doesn't requires time
        
        if ~isempty(Date)
            error('When using the Helio-Ecliptic coordinate system - Date/time must be empty');
        end
          
        Long = RA;
        Lat  = Dec;
        
    otherwise
        % non hec coordinates
        % time must be provided
        if isempty(Date)
            error('When using the non-Helio-Ecliptic coordinate system - Date/time must be provided');
        else
            if (size(Date,2)==1)
                %Date = Date;
            else
                Date = celestial.time.julday(Date);
            end
            
        end
        
        Coo = celestial.coo.coco([RA, Dec],InPar.CooSys,'e');
        Long = Coo(:,1);
        Lat  = Coo(:,2);
end

SpecZ    = ultrasat.zodiac_spectrum('InterpMethod',InPar.InterpMethodSpec,'BackType','zodi','Wave',InPar.Wave);

ZodiVmag = ultrasat.zodiac_bck_V(Long,Lat,Date,'InterpMethod',InPar.InterpMethodCoo);


Mag = [];
%AstroUtil.spec.synthetic_phot(SpecZ,'Johnson','V','Vega')
%MagVnorm = AstroUtil.spec.synphot(SpecZ,'Johnson','V','Vega')
MagVnorm = 22.1; % ultrasat.zodiac_spectrum should return spectrum normalized to V=22.1 (Vega)
CorrFactor = 10.^(-0.4.*(ZodiVmag-MagVnorm));

Ncol = max(Ncol,numel(Date(:)));

Spec.RA   = RA(:).'.*ones(1,Ncol);
Spec.Dec  = Dec(:).'.*ones(1,Ncol);
if isempty(Date)
    Spec.JD = [];
else
    Spec.JD   = Date(:).'.*ones(1,Ncol);
end
Spec.Wave = SpecZ(:,1);
Spec.Spec = SpecZ(:,2).*[CorrFactor(:).'];
MagBand   = AstroUtil.spec.synphot(SpecZ,InPar.FilterFamily,InPar.FilterName,InPar.FilterSys);
Spec.Mag  = MagBand + (ZodiVmag-MagVnorm);
Spec.MagV = ZodiVmag;

switch lower(InPar.OutType)
    case 'mat'
        % do nothing
    case 'astspec'
        N = numel(Spec.RA);
        AS = AstSpec(N,1);
        for I=1:1:N
            AS(I).Wave = Spec.Wave;
            AS(I).Int  = Spec.Spec(:,I);
        end
        Spec.Spec = AS;
    otherwise
        error('Unknown OutType option');
end


%AstroUtil.spec.synthetic_phot(Spec,InPar.FilterFamily,InPar.FilterName,InPar.FilterSys)
%AstroUtil.spec.synphot(Spec,InPar.FilterFamily,InPar.FilterName,InPar.FilterSys)
