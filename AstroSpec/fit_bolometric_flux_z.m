function [TotFlux,MinWL,MaxWL,F_lambda,BestFitTemp,BestFitRadius]=fit_bolometric_flux_z(Data,z,n,Ebv,dist,TempVec,plt)
%--------------------------------------------------------------------------
% fit_bolometric_flux function                                   AstroSpec
% Description: Estimate the Total bolometric flux according to given
%              magnitudes of different filters.
% Input  : - Data to fit in the format: (same as the data supplied
%            to fit_blackbody):
%            {FILTER_FAMILY, FILTER_NAME, MAG_TYPE, MAG, MAG_ERR}.
%          - z is the redshift of the source, input data is observed and
%            the output is in the source's rest frame
%          - n is a parameter governing the SED fit method:
%            n=0: interpolation (cubic)
%            n>0: polynomial fit of order n
%            n<0: blackbody fit, by calling fit_blackbody. 
%                 For this option TempVec and d are used.
%          - Known extinction E_{B-V} value with which to correct
%            the magnitude before the calculation.
%            If two element vector is given, then the second value is
%            R_{V}. Default is [0 3.08].
%          - Distance (in pc). This is required for the blackbody fit.
%          - Vector of temperatures. Required for the BB fit.
%            See fit_blackbody.m.
%          - plotting flag. plt=1 plots SED for reviewing.
% Output : - Total Flux (F_lambda, by integration over relevant WL range)
%          - MinWL (Ang) - Eff. WL of the red-most filter
%            (red side integral limit).
%          - MaxWL (Ang) - Eff. WL of the blue-most filter
%            (blue side integral limit).
%          - Flux array - F_lambda of the provided magnitudes.
%
% Tested : Matlab 7.14
%     By : Ofer Yaron                         Dec 2013
%          modified by Nora to include redshift (May 2015)
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [TotFlux,MinWL,MaxWL,Flux,BestFitTemp,BestFitRadius]=fit_bolometric_flux({'SDSS','u','AB',5,0.1;'SDSS','g','AB',3,0.1},2,0.035,[],[],1);
% Reliable: 1
%--------------------------------------------------------------------------

Def.z = 0;  % redshift set to 0
Def.n = 0;  % interpolation
Def.Ebv = [0 3.08];
Def.d = 10; % Distance in pc    % (for BB fit, if requested)
Def.TempVec=logspace(3,6,30).'; % (for BB fit, if requested)
Def.plt=0;

if (nargin==1),
   z   = Def.z;
   n   = Def.n;
   Ebv = Def.Ebv;
   dist   = Def.d;   
   TempVec = Def.TempVec;
   plt = Def.plt;
elseif (nargin==2),
   n   = Def.n;
   Ebv = Def.Ebv;
   dist   = Def.d;   
   TempVec = Def.TempVec;
   plt = Def.plt;
elseif (nargin==3),
   Ebv = Def.Ebv;
   dist   = Def.d;    
   TempVec = Def.TempVec;
   plt = Def.plt;
elseif (nargin==4),
   dist   = Def.d;    
   TempVec = Def.TempVec;
   plt = Def.plt;
elseif (nargin==5),
   TempVec = Def.TempVec;
   plt = Def.plt;
elseif (nargin==6),
   plt = Def.plt;
else
   % do nothing
end

if (length(Ebv)==1),
    Ebv(2) = Def.Ebv(2);
end

Col.FILTER_FAMILY  = 1;
Col.FILTER_NAME    = 2;
Col.MAG_TYPE       = 3;
Col.MAG            = 4;
Col.MAG_ERR        = 5;

c         = get_constant('c');
pc        = get_constant('pc');     % parsec [cm]
sig       = get_constant('sigma');
Rsun      = get_constant('SolR');
%Dist      = 10;                     % 10 [pc]
RAD       = 180/pi;

Inn = find(~isnan([Data{:,Col.MAG}]));
Data = Data(Inn,:)

TotFlux       = 0;
BestFitTemp   = 0;
BestFitRadius = 0;

Nfilt    = size(Data,1);
F_lambda = zeros(Nfilt,1);
EffWL    = NaN(Nfilt,1);

ObsMag   = ones(1)*[Data{:,Col.MAG}];
ObsMagE  = ones(1)*[Data{:,Col.MAG_ERR}];

for Ifilt=1:1:Nfilt,
   Amag=0;
   FW = get_filter(Data{Ifilt,Col.FILTER_FAMILY},Data{Ifilt,Col.FILTER_NAME});
   EffW = FW.eff_wl{1};
   if Ebv(1)>0,
       Amag = optical_extinction(Ebv(1),'B','V',EffW./10000,'C',Ebv(2));
   end
   %EffWL(Ifilt)=EffW;
   Mag=ObsMag(Ifilt)-Amag;
   %F_lambda=convert_flux(Mag,'AB','cgs/A',EffW,'A');
   %F_lambda(Ifilt)=convert_flux(Mag,Data{Ifilt,Col.MAG_TYPE},'cgs/A',EffW,'A');
   F_lambda(Ifilt)=convert_flux(Mag,Data{Ifilt,Col.MAG_TYPE},'cgs/A',EffW,'A')*(1+z);
   EffWL(Ifilt)=EffW/(1+z);
end

MinWL=min(EffWL(:));
MaxWL=max(EffWL(:));

wls=ceil(MinWL):1:floor(MaxWL); % WLs vector over range

if n>0          % fit order n polynomial
    p=polyfit(EffWL,F_lambda,n);
    fluxes = polyval(p,wls);
elseif n==0     % interpolation
    [sorted_WL,I]=sort(EffWL);
    sorted_Flux=F_lambda(I);
    fluxes=interp1(sorted_WL,sorted_Flux,wls,'cubic');
else            % n<0 : fit blackbody
    [RMS,Factor,Res]=fit_blackbody_spec(TempVec,[EffWL,F_lambda],'F_lambda','A',dist,Ebv(1),1);
    %[BestFitTemp,BestFitDeltaMag,BestFitAngRad,MinChi2,Dof,Chi2,MinRMS,Error]=fit_blackbody(TempVec,Data,Ebv);
    index=find(RMS==min(RMS));
    BestFitTemp=TempVec(index);
    BestFitRadius=sqrt(Factor(index)*(dist*pc).^2);
    
    %BestFitTemp=BestFitTemp*(1+z);
    [~,~,IlA]=black_body(BestFitTemp,wls);
    %IlA=IlA/(1+z);
    %BestFitRadius=BestFitAngRad*d*pc/RAD/3600;
    fluxes=IlA.*Factor(index); % F_lam = I_lam * R^2/d^2    
end

if (plt)
    clf
    plot(EffWL,F_lambda,'*b'); hold on;
    plot(wls,fluxes,'r');
    drawnow
    pause(1)
end

TotFlux=trapz(fluxes);

end
