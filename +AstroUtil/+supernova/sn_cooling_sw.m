function [L,Tc,R,t_d,Mag,t_min,t_max,t_opac,t_tr,t_delta]=sn_cooling_sw(Time,varargin)
%--------------------------------------------------------------------------
% sn_cooling_sw function                                           General
% Description: Calculate the shock cooling light curve (following the
%              shock breakout) of a supernova, based on the
%              Sapir & Waxman (2017) model (ApJ 838:130).
% Input  : - Vector of times [days] in which to calculate the shock
%            cooling luminosity in the obeserver frame.
%            Default is logspace(log10(1./24),log10(10),30).';
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Type' - {'rsg','bsg','wr'}. Default is 'rsg'.
%            'Vs' -   The shock velocity paramaeter in units of cm/s.
%                     Default is 10^8.5 cm/s (see eq. 2,3).
%            'E51'  - Energy in units of 10^51 erg. Default is 1.
%                     If empty will use Vs field, otherwise will override
%                     the Vs field. Default is empty.
%            'Rs'   - Progenitor radius in units of the solar radii.
%                     Default is 500.
%            'Ms'   - Ejecta mass in solar masses. Default is 10.
%            'f_rho'- f_{\rho} parameter. If empty use 1 for 'rsg',
%                     0.05 for 'bsg', and 0.1 for 'wr'.
%            'kappa'- Opacity. Default is 0.34 cm^2 gr^-1.
%            'Tcorr'- Temperature correction factor. 
%                     Default is 1.1 for RSG and 1.0 for BSG.
%            'redshift' - SN redshift. Default is 0.023061.
%            'Dist' - SN distance [Mpc]. 
%                     If empty will use the redshift field, otherwise
%                     will override the redshift field. Default is empty.
%            'Ebv'  - Extinction E(B-V). Default is 0.
%            'Rv'   - Selective extinction (R_V). Default is 3.08.
%            'FiltFam' - Filter family (see get_filter.m for options).
%                     Default is 'LIM'.
%                     If FiltFam is a struct in the format returned by
%                     get_filter.m, the input filter is used instead of
%                     calling get_filter.m
%            'FiltName'- Filter name (see get_filter.m for options).
%                     Default is 'JPL'.
%            'FiltSys' - Filter system {'AB','Vega'}. Default is 'AB'.
%            'Wave'    - Wavelength [A] of AB monochromatic magnitude to
%                     calculate. If provided the returned magnitude will
%                     corresponds to this wavelength (instead of the filter
%                     family and name). Default is empty.
%            'Model'- Determine weither to use SK17 extension to the
%                     original RW11 model by supressing the bolometric
%                     luminocity or not. Extending the model beyond the
%                     small delta approximation introduces dependencies on
%                     the progenitor's structure. Default is 'RW'.
%                     Use 'SW' to extend.
% Output : - Intrinsic luminosity [erg/s] as a function of time.
%          - Intrinsic color temperature [K] as a function of time.
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
%          - t_tr [days]   - the transperancy time which is the timescale
%                            for the luminosity suppression.
%          - t_delta [days]- the time where the small delta approximation
%                            breaks (delta ~ 0.1) and the emmision becomes
%                            dependent in the envelope density profile.
%                            The time is given in rest frame.
%
% Comment: The funtion support calculation for multiple cases in one call. 
%          Input parameters on the same dimension are considered as a set
%          for a single case (e.g. different pairs of progenitor radius and
%          and redshifts. Input parameters on different dimensions are
%          not considered as a part of a set and all the cases results from
%          the combinations of the different dimensions are calculated.
%          see the second example.
%
% Tested : Matlab R2017a
%     By : Noam Ganot                    Jan 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:  [L,Tc,R] = AstroUtil.supernova.sn_cooling_sw((3000:3000:90000)/86400);
%           [L,Tc,R,t_d,Mag,t_min,t_max,t_opac,t_tr,t_delta]=...
%                AstroUtil.supernova.sn_cooling_sw((3000:3000:90000)/86400,'Rs',...
%                      [500;1000],'Vs',shiftdim([0.5;1;2]*10^8.5,-2),...),...
%                      'redshift',[0.01;0.07],'FiltFam','Galex','FiltName','NUV')
% Reliable: 1
%--------------------------------------------------------------------------

SolR = constant.SunR;  % [cm]
SolM = constant.SunM;  % [g]
SigmaB = constant.sigma; % [erg cm-2 K-4]

Def.Time = logspace(log10(1./24),log10(10),30).'; % [days]
if (nargin==0)
    Time = Def.Time;
end
if (isempty(Time))
    Time = Def.Time;
end


DefV.Type      = 'rsg';
DefV.E51       = [];
DefV.Vs        = 10^8.5;
DefV.Ms        = 10;
DefV.Rs        = 500;
DefV.f_rho     = 1;
DefV.kappa     = 0.34;
DefV.Tcorr     = [];
DefV.redshift  = 0.023061;
DefV.Dist      = [];
DefV.Ebv       = 0;
DefV.Rv        = 3.08;
DefV.FiltFam   = [];
DefV.FiltName  = [];
DefV.FiltSys   = 'AB';
DefV.Wave      = [];
DefV.Model     = 'Original';

%InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

[S] = Util.array.bsx_nsize(size(Time),size(InPar.Vs),size(InPar.E51),size(InPar.Rs),size(InPar.f_rho),...
                           size(InPar.Ms),size(InPar.kappa),size(InPar.redshift),size(InPar.Dist),...
                           size(InPar.Tcorr),size(InPar.Ebv),size(InPar.Rv)                              );
                       
if (isempty(InPar.Dist))
    Dist = AstroUtil.cosmo.lum_dist(InPar.redshift);
else
    Dist = InPar.Dist.*1e6;  % [pc]
    InPar.redshift = AstroUtil.cosmo.inv_lum_dist(Dist);
end

% convert input Time to rest frame Time
t_d = Time./(1+InPar.redshift); 

switch lower(InPar.Type)
    case 'rsg'
        % Hydrogen/convective envelope
        % n=3/2:  
        %Eqs. 3, 4     
        beta = 0.191;
        if isempty(InPar.E51)
            v_sstar = InPar.Vs;
        else
            v_sstar = 1.05 .* InPar.f_rho.^-beta .* sqrt(InPar.E51.*1e51./(InPar.Ms.*SolM));
        end
        
        eps1 = 0.027; eps2 = 0.086; 
        T_phRW = 1.61 .* ( (v_sstar/10^8.5).^2 .* t_d.^2 ./ (InPar.f_rho .* InPar.Ms .* (InPar.kappa./0.34)) ).^ eps1 .* ...
            (InPar.Rs.*SolR./1e13).^0.25 ./ (InPar.kappa./0.34).^0.25 .* t_d.^-0.5;  % [eV]
%         T_phRW = 1.61 .* ( (v_sstar(:).'/10^8.5).^2 .* t_d.^2 ./ (InPar.f_rho(:).' .* InPar.Ms(:).' .* (InPar.kappa(:).'./0.34)) ).^ eps1 .* ...
%             (InPar.Rs(:).'.*SolR./1e13).^0.25 ./ (InPar.kappa(:).'./0.34).^0.25 .* t_d.^-0.5;  % [eV]
         
        if isempty(InPar.Tcorr)
            Tcorr = 1.1;
        else
            Tcorr = InPar.Tcorr;
        end        
        Tc    = convert.energy('eV','T',T_phRW.*Tcorr);
        
        L_RW   = 2e42 .* ( (v_sstar/10^8.5) .* t_d.^2 ./ (InPar.f_rho .* InPar.Ms .* (InPar.kappa./0.34)) ).^ -eps2 .* ...
            (v_sstar/10^8.5).^2 .* (InPar.Rs.*SolR./1e13) ./ (InPar.kappa./0.34);  % [erg/s]
%         L_RW   = 2e42 .* ( (v_sstar(:).'/10^8.5) .* t_d.^2 ./ (InPar.f_rho(:).' .* InPar.Ms(:).' .* (InPar.kappa./0.34)) ).^ -eps2 .* ...
%             (v_sstar(:).'/10^8.5).^2 .* (InPar.Rs(:).'.*SolR./1e13) ./ (InPar.kappa(:).'./0.34);  % [erg/s]

        % f_rho = (Menv/Mc)^0.5 for R = 3/2; Mej=M=Menv+Mc:
        Menv = (InPar.f_rho).^2./(1+(InPar.f_rho).^2).*InPar.Ms;
        %eq. 13,14:
        t_tr   = 19.5 .* sqrt((InPar.kappa./0.34).*Menv./(v_sstar/10^8.5));  % [days]       
        A = 0.94; a = 1.67; alpha = 0.8;
        if strcmpi(InPar.Model,'SW')
            LtoL_RW = A.*exp(-(a.*t_d./t_tr).^alpha);
            L = L_RW .* LtoL_RW;
        else
            L = L_RW;
        end
        
        % model validity limits:
        % SW Eq. 5 - RW validity.
        t_min   = 0.2.*(InPar.Rs.*SolR./1e13)./(v_sstar/10^8.5).* ...
                        max(0.5, (InPar.Rs.*SolR./1e13).^0.4./ ...
                            ((InPar.f_rho.*(InPar.kappa./0.34).*InPar.Ms).^0.2 .* (v_sstar/10^8.5).^0.7) );  % [days]
        t_delta = 3*InPar.f_rho.^-0.1.*sqrt((InPar.kappa./0.34).*InPar.Ms)./(v_sstar/10^8.5);  % [days]
        % Opacity condition - Tph,RW>0.7eV: SW Eqs. 24
        t_opac  = 7.4.*((InPar.Rs.*SolR./1e13)./(InPar.kappa./0.34)).^0.55;  % [days]
        t_max   = t_tr/a;  % [days] - text below eq. 18

    case 'bsg'
        % Hydrogen/radiative envelope
        % RW Eqs. 13b+15
        beta = 0.186;
        if isempty(InPar.E51)
            v_sstar = InPar.Vs;
        else
            v_sstar = 1.05 .* InPar.f_rho^-beta .* sqrt(InPar.E51.*1e51./(InPar.Ms.*SolM));
        end
        
        eps1 = 0.016; eps2 = 0.175; 
        
        T_phRW = 1.69 .* ( (v_sstar/10^8.5).^2 .* t_d.^2 ./ (InPar.f_rho .* InPar.Ms .* (InPar.kappa./0.34)) ).^ eps1 .* ...
            (InPar.Rs.*SolR./1e13).^0.25 ./ (InPar.kappa./0.34).^0.25 .* t_d.^-0.5;  % [eV]             
%         T_phRW = 1.61 .* ( (v_sstar(:)/10^8.5).^2 .* t_d.^2 ./ (InPar.f_rho(:) .* InPar.Ms(:) .* (InPar.kappa(:)./0.34)) ).^ eps1 .* ...
%             (InPar.Rs(:).*SolR./1e13).^0.25 ./ (InPar.kappa(:)./0.34).^0.25 .* t_d.^-0.5;  % [eV]             
        if isempty(InPar.Tcorr)
            Tcorr = 1.0;
        else
            Tcorr = InPar.Tcorr;
        end        
        Tc    = convert.energy('eV','T',T_phRW.*Tcorr);
        
        L_RW   = 2.1e42 .* ( (v_sstar/10^8.5) .* t_d.^2 ./ (InPar.f_rho .* InPar.Ms .* (InPar.kappa./0.34)) ).^ -eps2 .* ...
            (v_sstar/10^8.5).^2 .* (InPar.Rs.*SolR./1e13) ./ (InPar.kappa./0.34);  % [erg/s]
%         L_RW   = 2e42 .* ( (v_sstar(:)/10^8.5) .* t_d.^2 ./ (InPar.f_rho(:) .* InPar.Ms(:) .* (InPar.kappa(:)./0.34)) ).^ -eps2 .* ...
%             (v_sstar(:)/10^8.5).^2 .* (InPar.Rs(:).*SolR./1e13) ./ (InPar.kappa(:)./0.34);  % [erg/s]
        
        if true % Rc/R<<1
            Menv = InPar.f_rho.*(InPar.Ms)./(0.08+InPar.f_rho);
        else    % Rc/R~0.01
            Menv = InPar.f_rho.*(InPar.Ms)./(0.24+InPar.f_rho);
        end
        t_tr   = 19.5 .* sqrt((InPar.kappa./0.34).*Menv./(v_sstar/10^8.5));  % [d]
        
        A = 0.79; a = 4.57; alpha = 0.73;
        if strcmpi(InPar.Model,'SW')
            LtoL_RW = A.*exp(-(a.*t_d./t_tr).^alpha);
            L = L_RW .* LtoL_RW;
        else
            L = L_RW;
        end
        
        % model validity limits:
        % SW Eq. 5 - RW validity.
        t_min   = 0.2.*(InPar.Rs.*SolR./1e13)./(v_sstar/10^8.5).* ...
                        max( 0.5, (InPar.Rs.*SolR./1e13).^0.4./ ...
                                  ( (InPar.f_rho.*(InPar.kappa./0.34).*InPar.Ms).^0.2 .* (v_sstar/10^8.5).^0.7 ) );  % [d]
                        
        t_delta = 3*InPar.f_rho.^-0.1.*sqrt((InPar.kappa./0.34).*InPar.Ms)./(v_sstar/10^8.5);  % [d]
        % Opacity condition - Tph,RW>0.7eV: SW Eqs. 24
        t_opac  = 7.4.*((InPar.Rs.*SolR./1e13)./(InPar.kappa./0.34)).^0.55;  % [d]
        t_max   = t_tr/a;  % [d] - text below eq. 18
               
otherwise
    error('Unknown Progenitor option');
end

R      = sqrt(L./(4.*pi.*SigmaB.*Tc.^4));

if (nargout>3) && (~isempty(InPar.FiltFam) || ~isempty(InPar.Wave))
    Pc = constant.pc;
    if (isempty(InPar.Wave))
        if (ischar(InPar.FiltFam))
            Filter = AstFilter.get(InPar.FiltFam,InPar.FiltName);
            WaveRange = (Filter.min_wl: (Filter.max_wl - Filter.min_wl)./100 :Filter.max_wl).';
        elseif (AstFilter.isAstFilter(InPar.FiltFam))
            Filter = InPar.FiltFam;
            WaveRange = (Filter.min_wl: (Filter.max_wl - Filter.min_wl)./100 :Filter.max_wl).';        
        else
            Filter = InPar.FiltFam;
            Norm = trapz(Filter(:,1),Filter(:,2));
            Filter(:,2) = Filter(:,2)./Norm;
            WaveRange = (min(Filter(:,1)):range(Filter(:,1))./100:max(Filter(:,1))).';
        end
        
        redshiftedTc = Tc./(1+InPar.redshift);
        Spec = AstSpec.blackbody(redshiftedTc(:),WaveRange');
%        Spec = AstSpec.blackbody(redshiftedTc,WaveRange');
        F  = SigmaB .* (redshiftedTc).^4;

        SpecFactor = L./(4.*pi.*(Dist.*Pc).^2)./F;
        Spec=Spec.*SpecFactor(:);
        clear F

        Mag = synphot(Spec,Filter,[],InPar.FiltSys);
    else 
        [~,In]=AstroUtil.spec.black_body(Tc./(1+InPar.redshift),InPar.Wave).*(1+InPar.redshift).^4;
        Mag = convert.flux(In.*R.^2./(Dist.*Pc).^2,'cgs/Hz','AB',InPar.Wave,'A');
    end
    if (InPar.Ebv>0)
        A = AstroUtil.spec.extinction(InPar.Ebv,InPar.FiltFam,InPar.FiltName,InPar.Rv);
        Mag = Mag + A;
    end
    Mag = reshape(Mag,S);
else
    Mag = [];
end

% L  = reshape(L,shape);
% Tc = reshape(Tc,shape);
% R  = reshape(R,shape);
        



