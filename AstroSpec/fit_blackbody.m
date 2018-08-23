function [BestFitTemp,BestFitDeltaMag,BestFitAngRad,MinChi2,Dof,Chi2,MinRMS,Error]=fit_blackbody(TempVec,Data,Ebv)
%--------------------------------------------------------------------------
% fit_blackbody function                                         AstroSpec
% Description: Fit set of magnitudes to a black-body spectrum, and derive
%              its best fit parameters.
%              OBSOLETE: Use fit_bb.m instead.
% Input  : - Vector of temperatures to test (e.g., [1000:500:30000]') [K].
%            If empty matrix the use default. Default is logspace(3,6,30).'
%          - Data to fit in the format:
%            {FILTER_FAMILY, FILTER_NAME, MAG_TYPE, MAG, MAG_ERR}.
%          - Known extinction E_{B-V} value with which to correct
%            the magnitude before the fit.
%            If two element vector is given, then the second value is
%            R_{V}. Default is [0 3.08].
% Output : - Best fit temperature.
%          - Best fit magnitude difference - by how much the best fit
%            object is fainter than a black-body with the Sun radius
%            at 10pc and the same temperature.
%          - Best fit angular radius of black-body [arcsec].
%          - \chi^2 of best fit temperature.
%          - Number of degrees of freedom.
%          - Vector of \chi^2 value for each temperature.
%          - best fit RMS [mag]
%          - Error on best fit temperature [down up]
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Feb 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [BestFitTemp,BestFitDeltaMag,BestFitAngRad,MinChi2,Dof,Chi2,MinRMS]=fit_blackbody([1000:500:30000]',{'SDSS','u','AB',5,0.1;'SDSS','g','AB',3,0.1});
% Reliable: 2
%--------------------------------------------------------------------------

Def.Ebv = [0 3.08];
if (nargin==2),
   Ebv    = Def.Ebv;
end
if (length(Ebv)==1),
    Ebv(2) = Def.Ebv(2);
end

if (isempty(TempVec)),
   TempVec = logspace(3,6,30).';
end

Col.FILTER_FAMILY  = 1;
Col.FILTER_NAME    = 2;
Col.MAG_TYPE       = 3;
Col.MAG            = 4;
Col.MAG_ERR        = 5;

RAD       = 180./pi;
Pc        = get_constant('pc');     % parsec [cm]
Radius    = get_constant('SolR');   % solar radius [cm]
Dist      = 10;                     % 10 [pc]

Nfilt    = size(Data,1);
Ntemp    = length(TempVec);
CalcMag  = zeros(Ntemp,Nfilt);
ObsMag   = ones(Ntemp,1)*[Data{:,Col.MAG}];
ObsMagE  = ones(Ntemp,1)*[Data{:,Col.MAG_ERR}];
for Ifilt=1:1:Nfilt,
   CalcMag(:,Ifilt) = blackbody_mag_c(TempVec,Data{Ifilt,Col.FILTER_FAMILY},...
                             Data{Ifilt,Col.FILTER_NAME},...
                             Data{Ifilt,Col.MAG_TYPE},...
                             Radius,Dist);
   if (Ebv(1)>0),         
      FW = get_filter(Data{Ifilt,Col.FILTER_FAMILY},Data{Ifilt,Col.FILTER_NAME});
      EffW = FW.eff_wl{1};
      %Amag = optical_extinction(Ebv(1),'B','V',Data{Ifilt,Col.FILTER_NAME},'C',Ebv(2));
      Amag = optical_extinction(Ebv(1),'B','V',EffW./10000,'C',Ebv(2));
      CalcMag(:,Ifilt) = CalcMag(:,Ifilt) + Amag;
   end
end

% DeltaMag(Temp,Filt)
DeltaMag = ObsMag - CalcMag;

Npar = 2;
% fit for each temperature... (need to find mag offset)
BestDeltaMag = zeros(Ntemp,1);
for Itemp=1:1:Ntemp,
   BestDeltaMag(Itemp) = wmean(DeltaMag(Itemp,:)',ObsMagE(Itemp,:)');
end
BestDeltaMag = BestDeltaMag * ones(1,Nfilt);

Chi2 = nansum(((ObsMag - CalcMag - BestDeltaMag)./ObsMagE).^2,2);
[MinChi2,MinInd] = min(Chi2);
Resid = ObsMag - CalcMag - BestDeltaMag;
RMS   = nanstd(Resid,[],2);
MinRMS= RMS(MinInd);

BestFitTemp = TempVec(MinInd);
BestFitDeltaMag = BestDeltaMag(MinInd,1);
Dof      = Nfilt - Npar;

BestDist = Dist.*10.^(0.2.*BestFitDeltaMag);
BestFitAngRad = Radius./(BestDist.*Pc) .*RAD.*3600;

% estimate errors on temp
plot(TempVec,Chi2-(MinChi2+chi2inv(0.68,2)),'k-')
Extram = find_local_zeros(TempVec,Chi2-(MinChi2+chi2inv(0.68,2)));
%I = find(Extram(:,3)>0)
MaxBound = max(Extram(:,1));
MinBound = min(Extram(:,1));
if (MaxBound>BestFitTemp)
    ErrorUp = MaxBound-BestFitTemp;
else
    ErrorUp = NaN;
end
if (MinBound<BestFitTemp)
    ErrorDown = BestFitTemp-MinBound;
else
    ErrorDown = NaN;
end

Error = [ErrorDown, ErrorUp];

%plot([2316;4718;6184;7499],CalcMag(166,:),'o')
%hold on
%plot([2316;4718;6184;7499],ObsMag(166,:)-BestFitDeltaMag,'ro')
%std(DeltaMag(166,:))
