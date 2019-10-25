%see X-UV-updated note for derivations
global Silica
Silica=1
if Silica %SiO2
    n=1.5;
    rho=2.2;
else %Sapphire Al2O3
    n=1.75; %at 1 mu, n(0.25 mu)=1.85
    rho=4;
end

%My read off fig 6 of Kruk et al
F=10 .^[9 8 7 6 4 1]; 
Ee=[0.001 0.05 0.25 0.7 2.5 8];
%Yossi's data
%AE9_daily_electron_flux
AE9_dailymax_electron_flux
Ee1=Fmat(:,1); F1=Fmat(:,5);
%Compare Yossi's data with my read-off
figure(1)
loglog(Ee,F,Ee1,Fmat(:,2:5),'k')
axis([1e-2 10 10 1e9]); grid
ee=10 .^[log10(0.04):.01:log10(8)]; Fe=spline(Ee1,F1,ee);
%pause
%Calculate dEdX
dEdX_calc
%pause

%Cerenkov L per unit area and wavelength
bM=1/n; gM=1/sqrt(1-bM^2); EM=(gM-1)*0.511
ij=1:length(ee)-1; 
em=(ee(ij)+ee(ij+1))/2; gm=1+em/.511; bm=sqrt(1-1 ./gm.^2); 
ge=1+ee/.511; be=sqrt(1-1 ./ge.^2); 
Fm=spline(ee,Fe,em);
%j=-diff(Fe)./diff(ee)/pi; 
%loglog(em,j.*em,ee,Fe); grid
%pause
%jn=j/spline(em,j,1);
fC=max(0,1-1/n^2 ./bm.^2);
gEE=spline(Ek,1 ./dEdX,em);
%lE=(2/3)*(em.^(3/2)-Em^(3/2)).*(em<=1) + (em-1+2/3*(1-Em^(3/2))).*(em>1);
intg=gEE.*Fm.*fC; 
Lnorm=2*pi/137/rho/1e-4*spline(em,intg,1)
int=intg*diff(ee')/spline(em,intg,1)
Cint=[0 cumsum(intg.*diff(ee))]/(intg*diff(ee'));
L1mu=int*Lnorm % at lambda=1mu
L1mu/2/pi/n^2
figure(2)
the=acos(1/n ./be);
loglog(ee,Cint,ee,the); grid
axis([0.1 3 1e-2 1])
xlabel('E [MeV]'); ylabel('I(<E)/I, \theta')
%qC=(int_qC)*diff(ee')
pause

%Cerenkov I in -z direction, assuming all 
thM=acos(1/n)
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
figure(3)
qCth=(int_qC)*diff(th')
ICnorm=0.511*n/137/1e-4*spline(ee,Fe,1)*1/rho*spline(eth,gEE,1)
qCth*ICnorm/n^2
ICth=[0 cumsum(int_qC.*diff(th))]/qCth;
loglog(th,ICth,thm,eth,'k',th,(th./thM).^4,'b--'); grid
axis([0.2 1 1e-2 2])
xlabel('\theta'); ylabel('I(<\theta)/I, E_e [MeV]')


