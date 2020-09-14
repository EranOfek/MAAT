function Res=fit_polys_deltachi2(Data,varargin)
% Fit polynomials of various orders to a time series and calculate chi^2
% Package: timeseries
% Description: Fit polynomials of various orders to a time series and
%              calculate chi^2 and rms for each polynomial fit.
% Input  : - Data [time, magnitude, error]
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'MaxOrder'   - Maximum polynomial order. Default is 10.
%            'RejectProb' - Rejection probability for \Delta\chi^2
%                           Default is 0.95.
%            'ColT'       - Time colum index. Default is 1.
%            'ColM'       - magnitude colum index. Default is 2.
%            'ColE'       - Error colum index. Default is 3.
% Output : - Structure with fit information:
%            '
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Data=timeseries.fit_polys_deltachi2(Cat)
% Reliable: 
%--------------------------------------------------------------------------

DefV.MaxOrder             = 10;
DefV.RejectProb           = 0.95;
DefV.ColT                 = 1;
DefV.ColM                 = 2;
DefV.ColE                 = 3;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


T = Data(:,InPar.ColT);
M = Data(:,InPar.ColM);
E = Data(:,InPar.ColE);

MeanT  = mean(T);
RangeT = range(T);

T = (T-MeanT)./RangeT;

Npt = numel(T);
InPar.MaxOrder = min(Npt-1,InPar.MaxOrder);

VecOrder = (0:1:InPar.MaxOrder);
Norders  = numel(VecOrder);
Npar     = VecOrder + 1;


Mpred = zeros(Npt,Norders);
for I=1:1:Norders
    Deg   = VecOrder(I);
    
    Par   = polyfit(T,M,Deg);
    Mpred(:,I) = polyval(Par,T);
end

Resid = M - Mpred;
Res.Std   = std(Resid,[],1);
Chi2  = sum((Resid./E).^2,1);

Res.Chi2 = Chi2;

Chi2Factor = (Res.Chi2./(Npt-Npar));

Res.NewChi2   = Chi2./Chi2Factor;
Res.DeltaChi2 = chi2inv(InPar.RejectProb,Npar);

