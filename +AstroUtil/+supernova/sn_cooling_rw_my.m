function [L,Tc,R,Time,Mag]=sn_cooling_rw(Time,varargin)
% Shock cooling light curve (Rabinak & Waxman 2011)
% Package: AstroUtil.supernova
% Description: Calculate the shock cooling light curve (following the
%              shock breakout) of a supernova, based on the
%              Rabinak & Waxman (2011) model.
% Input  : - Vector of times [days] in which to calculate the shock
%            cooling luminosity.
%            Default is logspace(log10(1./24),log10(10),30).';
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Type' - {'rsg','bsg','wr'}. Default is 'rsg'.
%            'E51'  - Energy in units of 10^51 erg. Defaukt is 1.
%            'Ms'   - Ejecta mass in solar masses. Default is 10.
%            'Rs'   - Progenitor radius in units of the solar radii.
%                     Default is 500.
%            'f_rho'- f_{\rho} parameter. If empty use 0.1 for 'rsg',
%                     0.05 for 'bsg', and 0.1 for 'wr'.
%            'kappa'- Opacity. Default is 0.34 cm^2 gr^-1.
%            'Tcorr'- Temperature correction factor. Default is 1.2.
%            'Z'    - He abundabnce for WR models (Z=0 He; Z=1 no He).
%                     Default is 0.
%            'redshift' - SN redshift. Default is 0.023061.
%            'Dist' - SN distance [Mpc]. Default is 100.
%                     If empty will use the redshift field, otherwise
%                     will override the redshift field. Default is empty.
%            'Ebv'  - Extinction E(B-V). Default is 0.
%            'Rv'   - Selective extinction (R_V). Default is 3.08.
%            'FiltFam' - Filter family (see get_filter.m for options).
%                     Default is 'LIM'.
%            'FiltName'- Filter name (see get_filter.m for options).
%                     Default is 'JPL'.
%            'FiltSys' - Filter system {'AB','Vega'}. Default is 'AB'.
%            'Wave'    - Wavelength [A] of AB monochromatic magnitude to
%                     calculate. If provided the returned magnitude will
%                     corresponds to this wavelength (instead of the filter
%                     family and name). Default is empty.
% Output : - Intrinsic luminosity [erg/s] as a function of time.
%          - Intrinsic effective temperature [K] as a function of time.
%          - photospheric radius [cm] as a function of time.
%          - Time [days].
%          - Observed apparnt Magnitude as a function of time.
%            Corrected for redshift.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Dec 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [L,Tc,R,Time,Mag]=AstroUtil.supernova.sn_cooling_rw_my;
% Reliable: 2
%--------------------------------------------------------------------------
SolRad = 696000e5;  % [cm]

Def.Time = logspace(log10(1./24),log10(10),30).';
if (nargin==0)
    Time = Def.Time;
end
if (isempty(Time))
    Time = Def.Time;
end


DefV.Type      = 'rsg';
DefV.E51       = 1;
DefV.Ms        = 10;
DefV.Rs        = 500;
DefV.Z         = 0;
DefV.f_rho     = [];  %0.1;
DefV.kappa     = 0.34;
DefV.Tcorr     = 1.2;
DefV.redshift  = 0.023061;
DefV.Dist      = [];
DefV.Ebv       = 0;
DefV.Rv        = 3.08;
DefV.FiltFam   = 'LIM';
DefV.FiltName  = 'JPL';
DefV.FiltSys   = 'AB';
DefV.Wave      = [];

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

switch lower(InPar.Type)
    case 'rsg'
        % Hydrogen/convective envelope
        % n=3/2: RW Eqs. 13a+14
        if (isempty(InPar.f_rho))
            InPar.f_rho = 0.1;
        end
       
        L = 8.5e42.*InPar.E51.^0.92.*((InPar.Rs.*SolRad)./1e13).*(Time.*86400./1e5).^(-0.16)./...
               (InPar.f_rho.^0.27 .*InPar.Ms.^0.84.*(InPar.kappa./0.34).^0.92);
        Tc = 18567.*InPar.f_rho.^-0.037.*InPar.E51.^0.027.*((InPar.Rs.*SolRad)./1e13).^0.25.*(Time.*86400./1e5).^-0.45 ./...
               (InPar.Ms.^0.054.*(InPar.kappa./0.34).^0.28);
        Tc = Tc.*InPar.Tcorr;

    case 'bsg'
        % Hydrogen/radiative envelope
        % RW Eqs. 13b+15
        if (isempty(InPar.f_rho))
            InPar.f_rho = 0.05;
        end
        Tc = 18567.*InPar.f_rho.^-0.022 .* InPar.E51.^0.016 .* ((InPar.Rs.*SolRad)./1e13).^(1./4) .* ...
            (Time.*86400./1e5).^-0.47 .* InPar.Ms.^-0.033 .* (InPar.kappa./0.34).^-0.27;
        L  = 9.9e42.* InPar.E51.^0.85 .* ((InPar.Rs.*SolRad)./1e13) .*(Time.*86400./1e5).^-0.35 .* ...
            InPar.f_rho.^-0.16 .* InPar.Ms.^-0.69 .* (InPar.kappa./0.34).^-0.85;
        Tc = Tc.*InPar.Tcorr;
    case 'wr'
         % C/O envelope

         if (isempty(InPar.f_rho))
            InPar.f_rho = 0.1;
         end
       
        R12 = ((InPar.Rs.*SolRad)./1e12);
        t5  = (Time.*86400./1e5);

        if (InPar.Z==0)
            % RW Eq. 29 (R13->R12) + 27 (He+CO)

            InPar.Z = 0;  % He

            Tc  = convert.energy('eV','T',1.33).* InPar.f_rho.^-0.02 .* R12.^0.20.* t5.^-0.38;   % T_ph >= 1.07
            t5_b = (1.07./(1.33.*InPar.f_rho.^-0.02.*R12.^0.20)).^(-1./0.38);
            I   = find(Tc<convert.energy('eV','T',1.07));
            %t5(I(1))
            Tc(I) = convert.energy('eV','T',1.07).* (t5(I)./t5_b).^-0.12;   % T_ph <  1.07
   
             L = 3.3e42 .* InPar.E51.^0.84 .*R12.^0.85 .* InPar.f_rho.^-0.15 .* InPar.Ms.^-0.67 .* t5.^-0.03;

        else
            InPar.Z = 1;  % C/O (no He)

            Tc  = convert.energy('eV','T',1.33).* InPar.f_rho.^-0.02 .* R12.^0.20.* t5.^-0.38;   % T_ph >= 1.07
            t5_b = (1.07./(1.33.*InPar.f_rho.^-0.02.*R12.^0.20)).^(-1./0.38);
            I   = find(Tc<convert.energy('eV','T',1.07));
            % t5(I(1))
            Tc(I) = convert.energy('eV','T',1.07).* (t5(I)./t5_b).^-0.12;   % T_ph <  1.07

            L = 3.3e42 .* InPar.E51.^0.84 .*R12.^0.85 .* InPar.f_rho.^-0.15 .* InPar.Ms.^-0.67 .* t5.^-0.03;

            kappa_He    = 0.085 .* (1-InPar.Z).*(Tc./convert.energy('eV','T',1.07)).^0.88;   % T>= 1.07 eV
            kappa_He(I) = 0.085 .* (1-InPar.Z).*(Tc(I)./convert.energy('eV','T',1.07)).^10;   % T<  1.07 eV
   
            kappa_Z     = 0.043.*InPar.Z.*(Tc./convert.energy('eV','T',1.0)).^1.27;

            Ikappa = find(kappa_Z>kappa_He);

            Tc(Ikappa) = convert.energy('eV','T',1.5) .* InPar.f_rho.^-0.017 .* InPar.Z.^-0.2 .*...
                R12.^0.19 .*t5(Ikappa).^-0.35;
            % Note that according to RW erratum: t^+0.07 (and not -0.07)
            L(Ikappa)  = 4.7e42 .* InPar.E51.^0.83 .*R12.^0.8 .* InPar.f_rho.^-0.14.* InPar.Ms.^-0.67 .* ...
                t5(Ikappa).^0.07.*InPar.Z.^-0.63;
    end
    Tc = Tc.*InPar.Tcorr;

 otherwise
    error('Unknown Progenitor option');
end

SigmaB = constant.sigma;
R      = sqrt(L./(4.*pi.*SigmaB.*Tc.^4));

if (nargout>3)
    if (isempty(InPar.Dist))
        Dist = AstroUtil.cosmo.lum_dist(InPar.redshift);
    else
        Dist = InPar.Dist.*1e6;  % [pc]
        InPar.redshift = AstroUtil.cosmo.inv_lum_dist(Dist);
    end
    
    
    
    if (isempty(InPar.Wave))
        % calculate full black body magnitude
        Mag = AstroUtil.spec.blackbody_mag_c(Tc./(1+InPar.redshift),InPar.FiltFam,InPar.FiltName,InPar.FiltSys,R,Dist);
        Mag = Mag - 2.5.*log10((1+InPar.redshift)^4);
    else
        Pc = get_constant('pc');
        [~,In]=AstroUtil.spec.black_body(Tc./(1+InPar.redshift),InPar.Wave);
        In = In.*(1+InPar.redshift).^4;
        Mag = convert.flux(In.*R.^2./(Dist.*Pc).^2,'cgs/Hz','AB',InPar.Wave,'A');
    end
    if (InPar.Ebv>0)
        A = AstroUtil.spec.extinction(InPar.Ebv,InPar.FiltFam,InPar.FiltName,InPar.Rv);
        Mag = Mag + A;
    end
end
        



