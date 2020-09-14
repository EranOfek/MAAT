function [Res,ResEl]=Cerenkov(varargin) 
% Calculate the Cerenkov spectrum generated in a lens
% Package: telescope.Optics
% Description: Calculate the Cerenkov spectrum generated in a lens.
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'T' - Temperature [C]. Default is 20.
%            'Material' - Options are:
%                         'si02_suprasil_2a' (fused silica) default.
%            'FluxOption' - Options:
%                   DailyMin_MeanFlux: 
%                   DailyMax_MeanFlux: 
%                   DailyMin_95flux: 
%                   DailyMax_95flux:  - Default.
%                   DailyMax_50flux: 
%                   DailyMax_75flux: 
%                   DailyMin_50flux: 
%                   DailyMin_75flux: 
%            'Plot' - Default is false.
% Output : - A structure containing the following fields:
%            'Lam' - Wavelength [Ang].
%            'Int' - Intensity of Cerenkov radaition generated in the lens.
%                    Units: [count/cm^2/s/sr/micron]
%            'n'   - Refraction index at wavelength.
%          - A structure containing the following fields:
%            'E' - Electrons energy.
%            'F' - Electrons integrated flux with energy (>E)
%                  [counts(>E)/cm^2/s].
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res,ResEl]=ultrasat.Cerenkov
%%

RAD = 180./pi;


DefV.Material             = 'si02_suprasil_2a';
DefV.FluxOption           = 'DailyMax_95flux';
DefV.Plot                 = false;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Material   = InPar.Material;
FluxOption = InPar.FluxOption;
Plot       = InPar.Plot;

%see X-UV-updated note for derivations

switch lower(Material)
    case {'silica','sio2','si02_suprasil_2a'}
        
        %if Silica %SiO2
        n=1.5;
        rho=2.2;
    case 'sapphire'
        
        %else %Sapphire Al2O3
        n=1.75; %at 1 mu, n(0.25 mu)=1.85
        rho=4;
        
     
    otherwise
        error('Unknown Material option');
end
  
%My read off fig 6 of Kruk et al
F=10 .^[9 8 7 6 4 1]; 
Ee=[0.001 0.05 0.25 0.7 2.5 8];   % Electron energy [MeV]
%Yossi's data
%AE9_daily_electron_flux
%AE9_dailymax_electron_flux

[Fmat,ColCell] = ultrasat.geostat_electrons_spec_flux;
Col = cell2struct(num2cell(1:1:numel(ColCell)),ColCell,2);

Ee1=Fmat(:,1);
F1=Fmat(:,Col.(FluxOption));

ResEl.E = Ee1;
ResEl.F = F1;

%Compare Yossi's data with my read-off
if Plot
    figure;
    loglog(Ee,F,Ee1,Fmat(:,2:5),'k')
    axis([1e-2 10 10 1e9]);
    grid
end
ee=10 .^[log10(0.04):.01:log10(8)];  % energies in which to interpolate
Fe=spline(Ee1,F1,ee);
%pause
%Calculate dEdX
[Ek,dEdX]=ultrasat.dEdX_calc(Material,Plot);
%ultrasat.dEdX_calc
%pause

%Cerenkov L per unit area and wavelength
bM=1./n;
gM=1./sqrt(1-bM.^2);
EM=(gM-1).*0.511;
ij=1:length(ee)-1; 
em=(ee(ij)+ee(ij+1))./2;
gm=1+em./.511;
bm=sqrt(1-1 ./gm.^2); 
ge=1+ee./.511;
be=sqrt(1-1 ./ge.^2); 
Fm=spline(ee,Fe,em);
%j=-diff(Fe)./diff(ee)/pi; 
%loglog(em,j.*em,ee,Fe); grid
%pause
%jn=j/spline(em,j,1);
fC=max(0,1-1./n.^2 ./bm.^2);

% Ek and dEdX are from dEdX_calc

gEE=spline(Ek,1 ./dEdX,em);
%lE=(2/3)*(em.^(3/2)-Em^(3/2)).*(em<=1) + (em-1+2/3*(1-Em^(3/2))).*(em>1);
intg=gEE.*Fm.*fC; 
Lnorm=2*pi/137/rho/1e-4*spline(em,intg,1);
int=intg*diff(ee')/spline(em,intg,1);
Cint=[0 cumsum(intg.*diff(ee))]/(intg*diff(ee'));
L1mu=int*Lnorm; % at lambda=1mu
IC1mu = L1mu/2/pi/n^2;   % intensity at 1 micron ?

if Plot
    figure;
    the=acos(1/n ./be);

    FlagReal = imag(the)==0;

    loglog(ee(FlagReal),Cint(FlagReal));
    hold on;
    loglog(ee(FlagReal),the(FlagReal));
    grid
    axis([0.1 3 1e-2 1])
    xlabel('E [MeV]'); ylabel('I(<E)/I, \theta')
    %qC=(int_qC)*diff(ee')
    %pause
end

%Cerenkov I in -z direction, assuming all 
thM=acos(1/n);
th=thM*10 .^[-2:.01:0];
ik=1:length(th)-1; thm=(th(ik)+th(ik+1))/2;
bth=1/n ./cos(thm); gth=1 ./sqrt(1-bth.^2); eth=0.511*(gth-1);
%gEE=0.5*(eth.^(1/2).*(eth<=1) + 1 .*(eth>1));
gEE=spline(Ek,1 ./dEdX,eth);
%figure(3)
%loglog(Ek,1 ./dEdX,eth,gEE,eth,0.5*min(sqrt(eth),1));
%grid
%pause
int_qC=spline(ee,Fe,eth)/spline(ee,Fe,1).*(gth.*bth).^3.*sin(thm).^3 .*gEE/spline(eth,gEE,1);

qCth=(int_qC)*diff(th');
ICnorm=0.511*n/137/1e-4*spline(ee,Fe,1)*1/rho*spline(eth,gEE,1);
IC1mu_v2 = qCth*ICnorm/n^2;
ICth=[0 cumsum(int_qC.*diff(th))]/qCth;

if Plot
    figure;
    loglog(th,ICth,thm,eth,'k',th,(th./thM).^4,'b--'); grid
    axis([0.2 1 1e-2 2])
    xlabel('\theta'); ylabel('I(<\theta)/I, E_e [MeV]')
end

%[IC1mu, IC1mu_v2



% as a function of lambda
[Lam,n]=telescope.Optics.refraction_index('Material',Material);
Nn = numel(n);
IC1mu = nan(Nn,1);
L1mu = nan(Nn,1);

for In=1:1:Nn

    %Cerenkov L per unit area and wavelength
    bM=1./n(In);
    gM=1./sqrt(1-bM.^2);
    EM=(gM-1).*0.511;
    ij=1:length(ee)-1; 
    em=(ee(ij)+ee(ij+1))./2;
    gm=1+em./.511;
    bm=sqrt(1-1 ./gm.^2); 
    ge=1+ee./.511;
    be=sqrt(1-1 ./ge.^2); 
    Fm=spline(ee,Fe,em);
    %j=-diff(Fe)./diff(ee)/pi; 
    %loglog(em,j.*em,ee,Fe); grid
    %pause
    %jn=j/spline(em,j,1);
    fC=max(0,1-1./n(In).^2 ./bm.^2);

    % Ek and dEdX are from dEdX_calc

    gEE=spline(Ek,1 ./dEdX,em);
    %lE=(2/3)*(em.^(3/2)-Em^(3/2)).*(em<=1) + (em-1+2/3*(1-Em^(3/2))).*(em>1);
    intg=gEE.*Fm.*fC; 
   % Lnorm=2*pi/137/rho/(Lam(In).*1e-8)*spline(em,intg,1);  % per micron
   % Lnorm= Lnorm./(Lam(In)./1e4);   % wavelength is normalized to 1 micron
    
    Lnorm=2*pi/137/rho/( (1e-8.*Lam(In)).^2)*spline(em,intg,1); % count/cm^2/s/cm(wave)
    Lnorm=Lnorm.*1e-4;  % count/cm^2/s/micron(wave)
    
    int=intg*diff(ee')/spline(em,intg,1);
    Cint=[0 cumsum(intg.*diff(ee))]/(intg*diff(ee'));
    L1mu(In)  = int*Lnorm; % at lambda=1mu
    IC1mu(In) = L1mu(In)/2/pi/n(In)^2;   % intensity at 1 micron ?
    
end

Res.Lam = Lam;
Res.Int = IC1mu;
Res.Int_Units   = 'count/cm^2/s/sr/micron';
Res.IntAA       = Res.Int./1e4./((RAD.*3600).^2);  % count/cm^2/s/arcsec^2/Ang
Res.IntAA_Units = 'count/cm^2/s/arcsec^2/Ang';
Res.IntFA       = convert.flux(Res.IntAA,'ph/A','cgs/A',Res.Lam,'A');
Res.IntFA_Units = 'erg/cm^2/s/arcsec^2/Ang';
Res.n   = n;
Res.Lum = L1mu;

%Spec = 
