function [L,Tc,R,t_d,Mag,t_min,t_max,t_opac,t_tr,t_delta]=sn_cooling_sw(Time,varargin)
% Supernova shock cooling (Sapir & Waxman 2016)
% Package: AstroUtil.supernova
% Description: Calculate the shock cooling light curve (following the
%              shock breakout) of a supernova, based on the
%              Sapir & Waxman (2016) model (Draft version June 22,2016).
% Input  : - Vector of times [days] in which to calculate the shock
%            cooling luminosity in the obeserver frame.
%            Default is logspace(log10(1./24),log10(10),30).';
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Type' - {'rsg','bsg','wr'}. Default is 'rsg'.
%            'E51'  - Energy in units of 10^51 erg. Defaukt is 1.
%            'Ms'   - Ejecta mass in solar masses. Default is 10.
%            'Rs'   - Progenitor radius in units of the solar radii.
%                     Default is 500.
%            'f_rho'- f_{\rho} parameter. If empty use 1 for 'rsg',
%                     0.05 for 'bsg', and 0.1 for 'wr'.
%            'kappa'- Opacity. Default is 0.34 cm^2 gr^-1.
%            'Tcorr'- Temperature correction factor. 
%                     Default is 1.1 for RSG and 1.0 for BSG.
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
%          - photospheric effective radius [cm] as a function of time
%            (using T_col and not T_ph). 
%          - t_d [days] in the rest frame.
%          - Observed apparnt Magnitude as a function of time.
%            Corrected for redshift.
%          - t_min [days]  - the earliest time the model is valid given in
%                            rest frame.  
%          - t_max [days]  - the latest time the dynamical model is valid
%                            given in rest frame. 
%          - t_opac [days] - the time T_ph = 0.7eV and the constant opacity
%                            approximation breaks. The time is given in
%                            rest frame.
%          - t_delta [days]- the time where the small delta approximation
%                            breaks (delta ~ 0.1) and the emmision becomes
%                            dependent in the envelope density profile.
%                            The time is given in rest frame.
% Tested : Matlab R2016a
%     By : Noam Ganot                     Oct 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [L,Tc,R,t_d,Mag]=AstroUtil.supernova.sn_cooling_sw;
% Reliable: 1
%--------------------------------------------------------------------------

SolR = constant.SunR;  % [cm]
SolM = constant.SunM;  % [g]
SigmaB = constant.sigma;

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
DefV.f_rho     = 1;
DefV.kappa     = 0.34;
DefV.Tcorr     = [];
DefV.redshift  = 0.023061;
DefV.Dist      = [];
DefV.Ebv       = 0;
DefV.Rv        = 3.08;
DefV.FiltFam   = 'LIM';
DefV.FiltName  = 'JPL';
DefV.FiltSys   = 'AB';
DefV.Wave      = [];
DefV.Model     = 'extended';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% convert input Time to rest frame Time
if (isempty(InPar.Dist))
    Dist = AstroUtil.cosmo.lum_dist(InPar.redshift);
else
    Dist = InPar.Dist.*1e6;  % [pc]
    InPar.redshift = AstroUtil.cosmo.inv_lum_dist(Dist);
end

t_d = Time/(1+InPar.redshift);

switch lower(InPar.Type)
    case 'rsg'
        % Hydrogen/convective envelope
        % n=3/2: SW Eqs. 3, 4     
        beta = 0.191;
        v_sstar = 1.05 .* InPar.f_rho^-beta .* sqrt(InPar.E51.*1e51./(InPar.Ms.*SolM));
        
        eps1 = 0.027; eps2 = 0.086; 
        
        T_phRW = 1.61 .* ( (v_sstar/10^8.5).^2 .* t_d.^2 ./ (InPar.f_rho .* InPar.Ms .* (InPar.kappa./0.34)) ).^ eps1 .* ...
            (InPar.Rs.*SolR./1e13).^0.25 ./ (InPar.kappa./0.34).^0.25 .* t_d.^-0.5;  % [eV]             
        if isempty(InPar.Tcorr)
            Tcorr = 1.1;
        else
            Tcorr = InPar.Tcorr;
        end        
        Tc    = convert.energy('eV','T',T_phRW.*Tcorr);
        
        L_RW   = 2e42 .* ( (v_sstar/10^8.5) .* t_d.^2 ./ (InPar.f_rho .* InPar.Ms .* (InPar.kappa./0.34)) ).^ -eps2 .* ...
            (v_sstar/10^8.5).^2 .* (InPar.Rs.*SolR./1e13) ./ (InPar.kappa./0.34);  % [erg/s]
        t_tr   = 19.5 .* sqrt((InPar.kappa./0.34).*InPar.Ms./(v_sstar/10^8.5));  % [d]
        
        A = 0.94; a = 1.67; alpha = 0.8;
        LtoL_RW = A.*exp(-(a.*t_d./t_tr).^alpha);
        L = L_RW .* LtoL_RW;
        
        % model validity limits:
        % SW Eq. 5 - RW validity.
        t_min   = 0.2.*(InPar.Rs.*SolR./1e13)./(v_sstar/10^8.5).* ...
                        max(0.5, (InPar.Rs.*SolR./1e13).^0.4./ ...
                            ((InPar.f_rho.*(InPar.kappa./0.34).*InPar.Ms).^0.2 .* (v_sstar/10^8.5).^0.7) );  % [d]
        t_delta = 3*InPar.f_rho.^-0.1.*sqrt((InPar.kappa./0.34).*InPar.Ms)./(v_sstar/10^8.5);  % [d]
        % Opacity condition - Tph,RW>0.7eV: SW Eqs. 23
        t_opac  = 7.4.*((InPar.Rs.*SolR./1e13)./(InPar.kappa./0.34)).^0.55;  % [d]
        t_max   = t_tr/a;  % [d] - text below eq. 18

    case 'bsg'
        % Hydrogen/radiative envelope
        % RW Eqs. 13b+15
        beta = 0.186;
        v_sstar = 1.05 .* InPar.f_rho^-beta .* sqrt(InPar.E51.*1e51./(InPar.Ms.*SolM));
        
        eps1 = 0.027; eps2 = 0.086; 
        
        T_phRW = 1.61 .* ( (v_sstar/10^8.5).^2 .* t_d.^2 ./ (InPar.f_rho .* InPar.Ms .* (InPar.kappa./0.34)) ).^ eps1 .* ...
            (InPar.Rs.*SolR./1e13).^0.25 ./ (InPar.kappa./0.34).^0.25 .* t_d.^-0.5;  % [eV]             
        if isempty(InPar.Tcorr)
            Tcorr = 1.0;
        else
            Tcorr = InPar.Tcorr;
        end        
        Tc    = convert.energy('eV','T',T_phRW.*Tcorr);
        
        L_RW   = 2e42 .* ( (v_sstar/10^8.5) .* t_d.^2 ./ (InPar.f_rho .* InPar.Ms .* (InPar.kappa./0.34)) ).^ -eps2 .* ...
            (v_sstar/10^8.5).^2 .* (InPar.Rs.*SolR./1e13) ./ (InPar.kappa./0.34);  % [erg/s]
        t_tr   = 19.5 .* sqrt((InPar.kappa./0.34).*InPar.Ms./(v_sstar/10^8.5));  % [d]
        
        A = 0.94; a = 1.67; alpha = 0.8;
        LtoL_RW = ones(size(L_RW));
        LtoL_RW = A.*exp(-(a.*t_d./t_tr).^alpha);
        L = L_RW .* LtoL_RW;
        
        % model validity limits:
        % SW Eq. 5 - RW validity.
        t_min   = 0.2.*(InPar.Rs.*SolR./1e13)./(v_sstar/10^8.5).* ...
                        max(0.5, (InPar.Rs.*SolR./1e13).^0.4./ ...
                            ((InPar.f_rho.*(InPar.kappa./0.34).*InPar.Ms).^0.2 .* (v_sstar/10^8.5).^0.7) );  % [d]
        t_delta = 3*InPar.f_rho.^-0.1.*sqrt((InPar.kappa./0.34).*InPar.Ms)./(v_sstar/10^8.5);  % [d]
        % Opacity condition - Tph,RW>0.7eV: SW Eqs. 23
        t_opac  = 7.4.*((InPar.Rs.*SolR./1e13)./(InPar.kappa./0.34)).^0.55;  % [d]
        t_max   = t_tr/a;  % [d] - text below eq. 18
               
otherwise
    error('Unknown Progenitor option');
end

R      = sqrt(L./(4.*pi.*SigmaB.*Tc.^4));

if (nargout>3)
    if (isempty(InPar.Wave))
        % calculate full black body magnitude
        Mag = AstroUtil.spec.blackbody_mag_c(Tc(:)./(1+InPar.redshift),InPar.FiltFam,InPar.FiltName,InPar.FiltSys,R(:),Dist);
        Mag = Mag - 2.5.*log10((1+InPar.redshift).^4);
        Mag = reshape(Mag, size(Tc));
    else
        Pc = get_constant('pc');
        [~,In]=AstroUtil.spec.black_body(Tc./(1+InPar.redshift),InPar.Wave).*(1+InPar.redshift).^4;
        Mag = convert.flux(In.*R.^2./(Dist.*Pc).^2,'cgs/Hz','AB',InPar.Wave,'A');
    end
    if (InPar.Ebv>0)
        A = AstroUtil.spec.extinction(InPar.Ebv,InPar.FiltFam,InPar.FiltName,InPar.Rv);
        Mag = Mag + A;
    end
end
        



