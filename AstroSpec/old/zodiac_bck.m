function [ZodiMag,ZodiVmag,Flux,Counts]=zodiac_bck(Long,Lat,Date,Filter_family,Filter_name,FilterSys)
%--------------------------------------------------------------------------
% zodiac_bck function                                            AstroSpec
% Description: Calculate the zodiac magnitude and flux in a given filter
%              and sky position. The zodiac spectrum and position
%              dependent flux are adopted from the HST WFC3 handbook.
%              OBSOLETE: Use AstSpec.zodiac_bck instead.
% Input  : - Either Ecliptic longitude or Helio-ecliptic longitude
%            (i.e., L - L_sun) in radians
%            (see convertdm.m for additional options).
%            By default this should be the Helio-ecliptic longitude.
%          - Ecliptic latitude [radians].
%            (see convertdm.m for additional options).
%          - Date [D M Y] or JD at which to calculate the Helio-ecliptic
%            longitude. Defaut is empty. If empty than treat the first
%            input argument as Helio-ecliptic longitude.
%          - Filter family (e.g., 'GALEX'). See get_filter.m for
%            options. Default is 'SDSS'.
%          - Filter band name (e.g., 'NUV'). See get_filter.m for
%            options. Default is 'r'.
%          - Mag system {'AB'|'Vega'}. Default is 'AB'.
% Output : - Zodiacal light magnitude in filter and position.
%          - Zodiacal light V-band magnitude in position.
%          - Zodiacal light flux in filter [erg cm^-2 s^-1 arcsec^2]
%          - Zodiacal light photon rate [photons cm^-2 s^-1 arcsec^2]
% Tested : Matlab R2014a
%     By : Ilan Sagiv                      Sep 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [F,V,Flux,Counts]=zodiac_bck(45./RAD,80./RAD,[],'LIM','NUV')
% Reliable: NEED TO VERIFY zodi spec
%--------------------------------------------------------------------------
ReplaceNaNVal  = 21.5;

InterpMethod   = 'linear';
RAD            = 180./pi;
h              = get_constant('h','cgs');
c              = get_constant('c','cgs');
hc             = h*c*1e8; %[erg]*[Ang]

Def.Date          = [];
Def.Filter_family = 'SDSS';
Def.Filter_name   = 'r';
Def.FilterSys     = 'AB';
if (nargin==2),
    Date          = Def.Date;
    Filter_family = Def.Filter_family;
    Filter_name   = Def.Filter_name;
    FilterSys     = Def.FilterSys;
elseif (nargin==3),
    Filter_family = Def.Filter_family;
    Filter_name   = Def.Filter_name;  
    FilterSys     = Def.FilterSys;
elseif (nargin==4),
    Filter_name   = Def.Filter_name;   
    FilterSys     = Def.FilterSys;
elseif (nargin==5),
    FilterSys     = Def.FilterSys;    
else
    % do nothing
end
    

if (~isempty(Date)),
    % assume ecliptic longitude as input
    % Transform Ecliptic coordinates to Helio-ecliptic coordinates
    HelioEcLong = ecliptic2helioecliptic(Long,Date);
else
    HelioEcLong = Long;
end

% convert to 0-pi range
HelioEcLong(HelioEcLong>pi) = 2.*pi - HelioEcLong(HelioEcLong>pi);


% get sodi spectrum
Spec=zodiac_spectrum;


%wavelength=data(:,1);
%Zodiac    =data(:,3);
%Spectrum=[wavelength,Zodiac];

%if(0)% plot Spectrum
%	figure()
%	semilogy(wavelength,Zodiac,'r');
%	xlabel('Wavelength [A]')
%	ylabel('Flux [erg / sec\cdotcm^2 \cdot A]')
%    axis([1000 12e3 1e-25 1e-16])
%end


% approximate zodiacal sky background as a function of
% Heliocentric ecliptic longitude and ecliptic latitude
% (in Vmag / arcsec^2)
% from HST WFC3 Handbook; (table 9.4)
% V-band magnitude of zodi
% adopted from Table 9.4 in:
% http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c09_exposuretime08.html#389841
EcLon=(0:15:180)./RAD;
EcLat=(0:15:90)./RAD;
%lat-->
%v lon
ZodiVmagTable=[...
nan ,nan ,nan ,nan ,22.6,23.0,23.3;
nan ,nan ,nan ,nan ,22.6,23.1,23.3;
nan ,nan ,nan ,22.3,22.7,23.1,23.3;
nan ,nan ,22.1,22.5,22.9,23.1,23.3;
21.3,21.9,22.4,22.7,23.0,23.2,23.3;
21.7,22.2,22.6,23.9,23.1,23.2,23.3;
22.0,22.3,22.7,23.0,23.2,23.3,23.3;
22.2,22.5,22.9,23.1,23.3,23.3,23.3;
22.4,22.6,22.9,23.2,23.3,23.3,23.3;
22.4,22.6,22.9,23.2,23.3,23.3,23.3;
22.4,22.6,22.9,23.1,23.3,23.3,23.3;
22.3,22.5,22.8,23.0,23.2,23.4,23.3;
22.1,22.4,22.7,23.0,23.2,23.4,23.3];

% Interpolate zodi Vmag to requested positions
ZodiVmag = interp2(EcLat,EcLon,ZodiVmagTable,abs(Lat),HelioEcLong,InterpMethod);
ZodiVmag(isnan(ZodiVmag)) = ReplaceNaNVal;

Vmag0 = 22.1;   % synphot(Spec,'Johnson','V','Vega')
DeltaMag = ZodiVmag-Vmag0;

ZodiMag = DeltaMag + synphot(Spec,Filter_family,Filter_name,FilterSys);

Filter = get_filter(Filter_family,Filter_name);
FilterTr = [Filter.nT{1}(:,1), Filter.nT{1}(:,2)./max(Filter.nT{1}(:,2))];
[Flux,Counts] = spec_photon_counts([Spec(:,1),Spec(:,2).*10.^(-0.4.*DeltaMag)],FilterTr,[],1,1./3.08e18);

