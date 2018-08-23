function [L,Res]=macronova_multizone(t,varargin)
% Calculate the Waxman et al. analytic multizone kilonova model 
% Package: AstroUtil.supernova
% Description: Calculate the Waxman et al. analytic multizone kilonova
%              model.
% Input  : - time [seconds]. If empty or not provided, then default is
%            logspace(2,6,10000).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Case' - Select predefined private cases:
%                     '1a' - XM<0 & s>0
%                     '1b' - XM<0 & s<0
%                     '2a' - XM>0 & s>0
%                     '2b' - XM>0 & s<0
%            'M'    - Ejecta mass [solar mass].
%            'vM'   - Ejecta v_M velocity [speed of light].
%            'Alpha'- Ejecta velocity distribution: v(m)=vM m^-Alpha.
%            'Beta' - Heating power-law decay.
%            'Gamma'- Evolution in kappa.
%            't0'   - t0.
%            'kM'   - kappa(t=tM) [cm^2/gr]
%            'kG'   - kappa_\gamma [cm^2/gr]
%            'epsdot0' - Energy deposition rate at t0 [erg/s/gr].
% Output : - Luminosity [erg/s].
%          - Structure array with additional parameters.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: t=logspace(1,7,30)'; [L,Res]=AstroUtil.supernova.kilonova_multizone([]);
% Reliable: 2
%--------------------------------------------------------------------------

DefV.Case                 = [];
DefV.M                    = 5e-2;  % solar mass
DefV.vM                   = 0.15;   % fraction of speed of light
DefV.Alpha                = 0.7;
DefV.Beta                 = 1.0;
DefV.Gamma                = 0.6;
DefV.t0                   = 300;     % [s]
DefV.kM                   = 0.3;     % [cm^2/gr]
DefV.ke                   = 0.4;       % [cm^2/gr]
DefV.epsdotM              = 0.6e10; %convert.energy('eV','erg',1e6)./constant.mp./(86400.*10)./300;
DefV.LM                   = [];
DefV.XM                   = [];
InPar = InArg.populate_keyval(DefV,varargin,mfilename);
        
M      = InPar.M.*constant.SunM;
vM     = InPar.vM.*constant.c;   % convert velocity from [c] to [cm/s]
Alpha  = InPar.Alpha;
Beta   = InPar.Beta;
Gamma  = InPar.Gamma;
t0     = InPar.t0;
kM     = InPar.kM;
ke     = InPar.ke;
epsdotM= InPar.epsdotM;


% model demonstration / example
if (nargin==0)
    t = [];
end
if (isempty(t))
    t = logspace(2,6,10000)';
end

if (~isempty(InPar.Case))
    switch InPar.Case
        case '2a'
            % XM>1 & s>0
            kM    = 5;
            ke    = 0.1;
            Alpha = 0.4;
            Beta  = 1.1;
            Gamma = 0.0;
            %Text  = 'case 1a: $X_{M}>1$ \& $s>0$';
        case '2b'
            % XM>1 & s<0
            kM    = 5;
            ke    = 0.1;
            Alpha = 0.4;
            Beta  = 1.1;
            Gamma = 1;
            %Text  = 'case 1b: $X_{M}>1$ \& $s<0$';
        case '1a'
            % XM<1 & s>0
            kM    = 0.3;
            ke    = 1;
            Alpha = 0.4;
            Beta  = 1.1;
            Gamma = 0;

        case '1b'
            % XM<1 & s<0
            kM    = 0.1;
            ke    = 1;
            Alpha = 0.4;
            Beta  = 1.1;
            Gamma = 1;        

        otherwise
    end
end
        

tM    = (kM.*M./(4.*pi.*(1+2.*Alpha).*constant.c.*vM)).^(1./2);

%epsdot = epsdot0.*exp(-t./t0);
%It0     = t>t0;
%epsdot(It0) = epsdot0.*(t(It0)./t0).^(-Beta);
epsdot = epsdotM.*(t./tM).^(-Beta);

%teps  = sqrt(kG.*M./(4.*pi.*(1+2.*Alpha).*vM.^2 ));
teps  = (constant.c.*ke.*M./(4.*pi.*Alpha.*(2+3.*Alpha).*vM.^3)).^(1./2);

%epsdotM     = epsdot0.*(tM./t0).^(-Beta);
if (isempty(InPar.LM))
    LM    = M.*epsdotM;
else
    LM    = InPar.LM;
end
rph   = constant.c.*tM.*(vM./constant.c).^((1+Alpha)./(1+2.*Alpha)) .*(t./tM).^((1+Alpha.*Gamma)./(1+2.*Alpha));
rphtM = constant.c.*tM.*(vM./constant.c).^((1+Alpha)./(1+2.*Alpha)) .*(tM./tM).^((1+Alpha.*Gamma)./(1+2.*Alpha));
TeffM = (LM./(4.*pi.*constant.sigma.*rphtM.^2)).^(1./4);
%s     = (2.*Alpha - (1+2.*Alpha).*Gamma)./(1+Alpha);
s     = (4.*Alpha - (1+3.*Alpha).*Gamma)./(1+Alpha);
Eta   = 1+ (2-Gamma)./((1+Alpha).*(2-Beta));
%XM    = kM.*vM./(kG.*constant.c);
if (isempty(InPar.XM))
    XM    = (tM./teps).^2;
else
    XM    = InPar.XM;
end
%tGD   = tM.*XM.^(1./s);
teD   = tM.*XM.^(1./s);

%[XM,s]
%[t0,tGD, tM, teps]./86400

L = nan(size(t));
if (tM<teps)
    I = t<teD & s>0;
    %L(I) = LM.*Eta.*XM.^(-1+Beta./2) .*(t(I)./tM).^(2.*(1-Gamma)-Beta-s.*Beta./2);
    if (s>0)
        I = t<teD;
        L(I) = LM.*Eta.*XM.^(-1+Beta./2) .* (t(I)./tM).^( (2-Gamma)./(1+Alpha) - Beta+s.*(1-Beta./2) );
        I = t<tM & t>teD;
        L(I) = LM.*Eta.*(t(I)./tM).^( (2-Gamma)./(1+Alpha)-Beta );
    else
        I = t<tM;
        L(I) = LM.*Eta.*(t(I)./tM).^( (2-Gamma)./(1+Alpha)-Beta );
    end
    I = t>tM & t<teps;
    L(I) = LM.*(t(I)./tM).^(-Beta);
    I = t>teps;
    L(I) = LM.*XM.^(Beta./2) .*(t(I)./teps).^(-Beta-2);
else
    % tM>teps
    if (s<0)
        I = t<teD;
        L(I) = LM.*Eta.*(t(I)./tM).^((2-Gamma)./(1+Alpha) - Beta);
        I = t<tM & t>teD;
        L(I) = LM.*Eta.*XM.^(-1+Beta./2) .* (t(I)./tM).^ ( (2-Gamma)./(1+Alpha) - Beta+s.*(1-Beta./2) );
    else
        I = t<tM;
        L(I) = LM.*Eta.*XM.^(-1+Beta./2) .* (t(I)./tM).^ ( (2-Gamma)./(1+Alpha) - Beta+s.*(1-Beta./2) );
    end
    I = t>tM;
    L(I) = (LM./XM).*(t(I)./tM).^(-Beta-2);
    
end


% photospheric radius before t_eps/t_ph
Res.rph  = constant.c.*tM.*(vM./constant.c).^((1+Alpha)./(1+2.*Alpha)) .*(t/tM).^((1+Gamma.*Alpha)./(1+2.*Alpha));
Res.T    = (L./(4.*pi.*constant.sigma.*Res.rph.^2)).^(1./4);

%loglog(t./86400,L)

Res.t     = t;
Res.XM    = XM;
Res.s     = s;
Res.t0    = t0;
Res.teD   = teD;
Res.tM    = tM;
Res.teps  = teps;
Res.Alpha = Alpha;
Res.Beta  = Beta;
Res.Gamma = Gamma;
Res.epsdotM = epsdotM;
Res.rphtM = rphtM;
Res.TeffM = TeffM;
Res.LM    = LM;


