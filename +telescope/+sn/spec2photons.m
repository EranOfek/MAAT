function [Flux,Counts]=spec2photons(Spectrum,Filter,Extin,Radius,Dist)
% Spectrum to photon counts in band.
% Package: telescope.sn
% Description: Given a spectrum and the effective area of an instrument as
%              a function of wavelength, calculate the the total recieved
%              flux and the photons count rate in the instrument.
% Input  : - Spectra of the source (emitted from a cm^2 on the source):
%            [Wavelengh(Ang), Emmitence(erg cm^-2 s^-1 A^-1)].
%            Alternatively, a scalar representing a black-body
%            temperature [K].
%            Alternatively this can be an AstSpec object.
%          - Filter [Wavelength(Ang), EffectiveAreaOfInstrument(cm^2)].
%            If two element vector than assumes it to be the wavelength
%            range [Ang] of a top-hat filter with efective area of 1 cm^2.
%          - Transmission law [Wavelength(Ang), Transmission(frac)].
%            Where the extinction is measured in fraction
%            (e.g., 1 means no extinction; 0.5 means half the photons
%             are absorbed).
%            If empty matrix (i.e., []), then assumes no extinction.
%          - Source radius [cm], default is 1cm.
%          - Source distance [pc], default is 10pc.
% Output : - Observed flux [erg s^-1] on instrument.
%          - Observed photon count rate on instrument.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Nov 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See Also: blackbody_mag_c.m, blackbody_flux.m
% Example: [Flux,Counts]=telescope.sn.spec2photons(5770,[100 300000],[],696000e5,10);
%          Flux.*4.*pi.*(10.*3.08e18).^2;   % Bolometric luminosity of the Sun
% Reliable: 2
%--------------------------------------------------------------------------
Nstep        = 1000;
InterpMethod = 'linear';

Pc  = constant.pc;


DefExtin  = [];
DefRadius = 1;
DefDist   = 10;
if (nargin==2)
   Extin  = DefExtin;
   Radius = DefRadius;
   Dist   = DefDist;
elseif (nargin==3)
   Radius = DefRadius;
   Dist   = DefDist;
elseif (nargin==4)
   Dist   = DefDist;
elseif (nargin==5)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (numel(Filter)==2)
   Step = abs(range(Filter))./Nstep;
   Filter = [min(Filter):Step:max(Filter)]';
   Filter = [Filter, ones(size(Filter))];
end

if (AstSpec.isastspec(Spectrum))
    Spectrum = astspec2mat(Spectrum,'mat');
end

if (numel(Spectrum)==1)
   % Assumes black-body spectrum
   %[Il,In,IlA,ImJy,Ip] = black_body(Spectrum,Filter(:,1));
   Spectrum = AstSpec.blackbody(Spectrum,Filter(:,1),'cgs/A','ang','mat');
   %Spectrum = [Filter(:,1), IlA];
end

% Equalize Filter and Spectrum:
[NewSpec,NewFilt]=AstroUtil.spec.eq_sampling(Spectrum,Filter,[],InterpMethod);

% applay extinction:
if (isempty(Extin)==1)
   Extin = [NewSpec(:,1), ones(size(NewSpec,1),1)];
end
[NewSpec,NewExtin] = AstroUtil.spec.eq_sampling(NewSpec,Extin,NewSpec(:,1),InterpMethod);

% Observd Spectrum:
ObsSpec = [NewSpec(:,1), NewSpec(:,2) .* 4.*pi.* Radius.^2 ./(4.*pi.*(Dist.*Pc).^2)];

% Applay filter + extinction:  
ObsSpec = [ObsSpec(:,1), ObsSpec(:,2).*NewFilt(:,2).*NewExtin(:,2)];   % [erg s^-1 A^-1]

% Flux:
Flux    = trapz(ObsSpec(:,1), ObsSpec(:,2));


% Count rate:
EnConv  = convert.energy('A','erg',ObsSpec(:,1));     % conversion between erg to Ang.
Counts  = trapz(ObsSpec(:,1), ObsSpec(:,2)./EnConv);
