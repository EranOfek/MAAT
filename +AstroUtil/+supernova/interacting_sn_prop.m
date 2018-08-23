function Data=interacting_sn_prop(L0,Alpha,m,varargin)
% Calculate the properties of CSM-ejecta interacting SN
% Package: AstroUtil.supernova
% Description: Calculate the properties of CSM-ejecta interacting SN based
%              on the Ofek et al. (2014) formulae.
% Input  : - L0
%          - Alpha - Bolometric light curve power-law index slope.
%          - m - stellar envelope power law index
%          * Arbitary number of pairs of arguments: ...,key,val,...
%            The following keywords are available:
%            'tbo' - Wind shock breakout time scale [s].
%            'Eps' - Efficiency.
%            'kappa'
%            'mu_p'
%            'Te_shock'
%            'Te_csm'
%            'Z'
%            'Nu'
%            'tI'
% Output : - Structure pf output parameters.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Data=AstroUtil.supernova.interacting_sn_prop(1e44,-1./3,8)

%DefV.m        = 10;
%DefV.alpha    = -0.3;
DefV.tbo      = 15.*86400;  % [s]
DefV.Eps      = 0.5;
DefV.kappa    = 0.34;
DefV.mu_p     = 0.6;
DefV.Te_shock = 1e8;
DefV.Te_csm   = 1e4;
DefV.Z        = 1;   % atomic number
DefV.Nu       = 10;  % [GHz]
DefV.tI       = 300.*86400; % [s]  time to which to integrate mass
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

tbo      = InPar.tbo; 
Eps      = InPar.Eps;
kappa    = InPar.kappa;
mu_p     = InPar.mu_p;
Te_shock = InPar.Te_shock;
Te_csm   = InPar.Te_csm;
Z        = InPar.Z;
Nu       = InPar.Nu;
tI       = InPar.tI;


c  = constant.c;
mp = constant.mp;

% solve('alpha=((2-w)*(m-3)+3*(w-3))/(m-w)','w')
w   = (Alpha.*m - 2.*m + 15)./(Alpha - m + 6);

vbo = tbo.^((Alpha-1)./3).*(2.*pi.*Eps.*(m-w).*(w-1)./(m-3).*c./(kappa.*L0)).^(-1./3);

rbo = vbo.*tbo.*(m-w)./(m-3);

K = c./kappa .* (w-1).*((m-w)./(m-3)).^(w-1) .*vbo.^(w-2) .*tbo.^(w-1);
% or
K = L0./(2.*pi.*Eps) .* ((m-w)./(m-3)).^(w-2) .* vbo.^(w-5) .*tbo.^(Alpha+w-2);


PI1 = ((m-w).*(3.*w-4)+(w-3).*(3.*w-6)+(2-w).*(m-3))./(3.*(m-w));
PI2 = (3-m).*(w-3)./(m-w);

rI        = rbo.*(tI./tbo).^((m-3)./(m-w));
vI        = vbo.*(tI./tbo).^((w-3)./(m-w));
tauI      = kappa.*K.*rI.^(1-w)./(w-1);
nI        = K./(mu_p.*mp) .* rI.^(-w);
tffcoolI  = 1.8e15 .* (Te_shock./1e8).^(1./2) .* nI.^(-1) .* Z.^(-2);
NI        = K.*rI.^(1-w)./(mu_p.*mp.*(w-1));
tauffI    = 8.5e-28.*(Te_csm./1e4).^(-1.35) .* (Nu./10).^(-2.1) .* K.^2.*rI.^(1-2.*w)./(mu_p.^2 .* mp.^2 .* (2.*w-1));
zetaI     = L0.*tI.^(Alpha)./(nI.*rI.^2) .* max(tauI, ones(size(tauI))); % ionization parameter

MI = 4.*pi.*c.* ((m-w)./(m-3)).^(5./3) .* (w-1).^(2./3)./(3-w) .* (2.*pi.*c.*Eps./L0).^(-1./3) .* kappa.^(-2./3) .* tbo.^PI1 .* tI.^PI2;

Data.In.m     = m;
Data.In.Alpha = Alpha;
Data.In.L0    = L0;
Data.Par.w    = w;
Data.Par.vbo = vbo;
Data.Par.rbo = rbo;
Data.Par.K   = K;
Data.Ev.rI   = rI;
Data.Ev.vI   = vI;
Data.Ev.tauI = tauI;
Data.Ev.nI   = nI;
Data.Ev.tffcoolI   = tffcoolI;
Data.Ev.NI   = NI;
Data.Ev.tauffI   = tauffI;
Data.Ev.zetaI   = zetaI;
Data.Ev.MI   = MI;