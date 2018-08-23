function [XC,Res]=xcorr(TS1,TS2,Lag,varargin)
% Calculate the \chi2 and cross correlation between two time series.
% Package: timeseries
% Description: Calculate the \chi2 and cross correlation between two time
%              series.
% Input  : - First time series. Matrix of [Time, Mag, Error], where Error
%            is optional. If Error column is not provided then will be set
%            to 1. Note that Error is used only in the \chi2^2 calculation
%            and not for the cross-correlation.
%          - Second time series.
%          - Vector of lags to check.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SubMeanChi2'- For each lag, subtract mean value of
%                           difference between the two series.
%                           Default is true.
%            'SubMeanCorr'- For each lag, subtract the mean of each series
%                           before cross-correlation.
%                           Default is true.
%            'ErrConf'    - Percentiles for condidence intervals.
%                           Default is 1 - 2.*normcdf([1;2;3],0,1,'upper')
%            'ColT'       - Time column in input time series. Default is 1.
%            'ColM'       - Mag column in input time series. Default is 2.
%            'ColE'       - Err column in input time series. Default is 3.
%            'InterpMethod'- Series interpolation method.
%                           Default is 'linear'.
%            'InterpMethodErr'- Error interpolation method.
%                           Default is 'nearest'.
% Output : - A matrix with [Lag, Chi2, Dof, cross-correlation].
%          - A structure with the following information:
%            'Npar'      - Number of free parameters (2).
%            'MeanDiffM' - Vector of the Mean difference between the two
%                          series for each lag.
%            'MinChi2'   - Minimum \chi^2 value.
%            'MinChi2Ind'- Lag index of min \chi^2.
%            'MinChi2Lag'- Lag of min \chi^2.
%            'MinDof'    - Dof at min \chi^2.
%            'MaxCorr'   - Maximum cross-correlation.
%            'MaxCorrInd'- Lag index of max cross-correlation.
%            'MaxCorrLag'- Lag of max cross-correlation.
%            'DeltaChi2'
%            'AdjustedChi2'
%            'CI'
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [T,F,TotalF,X,Y]=AstroUtil.lensing.generate_timedelay_lc;
%          TS1=[T,F(:,1)];                                          
%          TS2=[T,F(:,2)];
%          [XC,R]=timeseries.xcorr(TS1,TS2,(-100:1:100)');     
% Reliable: 2
%--------------------------------------------------------------------------

DefV.SubMeanChi2          = true;
DefV.SubMeanCorr          = true;
DefV.ErrConf              = 1 - 2.*normcdf([1;2;3],0,1,'upper');  % 1,2,3 sigma probabilities
DefV.ColT                 = 1;
DefV.ColM                 = 2;
DefV.ColE                 = 3;
DefV.InterpMethod         = 'linear';
DefV.InterpMethodErr      = 'nearest';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% make Lag a row vector:
Lag = Lag(:).';

[Nrow1,Ncol1] = size(TS1);
[Nrow2,Ncol2] = size(TS2);

if (Ncol1==2)
    % set errors to 1
    TS1 = [TS1, ones(Nrow1,1)];
end
if (Ncol2==2)
    % set errors to 1
    TS2 = [TS2, ones(Nrow2,1)];
end

T1 = TS1(:,InPar.ColT);
T2 = TS2(:,InPar.ColT);
M1 = TS1(:,InPar.ColM);
M2 = TS2(:,InPar.ColM);
E1 = TS1(:,InPar.ColE);
E2 = TS2(:,InPar.ColE);

% use bsxfun in oldder versions:
% T1 is a column vector, Lag is a row vector
MatLagT1 = T1 + Lag;


LagM2 = interp1(T2,M2,MatLagT1,InPar.InterpMethod);
LagE2 = interp1(T2,E2,MatLagT1,InPar.InterpMethod);

%MeanM1    = nanmean(M1);
%MeanLagM2 = nanmean(LagM2);  % for each column

%StdM1     = nanstd(M1);
%StdLagM2  = nanstd(LagM2);   % for each column

% number of free parameters (Lag)
Res.Npar = 1;
% calculate the difference mean for each lag
DiffM    = LagM2 - M1;
Res.MeanDiffM = nanmean(DiffM); % for each lag
% subtract difference mean from difference series for each lag
if (InPar.SubMeanChi2)
    DiffM    = DiffM - Res.MeanDiffM;
    Res.Npar = Res.Npar + 1; % add mag diff as a free par
end

% calculate \chi2 for all lags
%Chi2ind = ((LagM2 - MeanLagM2) - (M1-MeanM1)).^2./(E1.^2 + LagE2.^2);
Chi2ind = DiffM.^2./(E1.^2 + LagE2.^2);
Dof     = sum(~isnan(Chi2ind));   % for each Lag
Chi2    = nansum(Chi2ind);        % for each lag

% calculate xcorr for all lags
if (InPar.SubMeanCorr)
   LagM2 = LagM2 - nanmean(LagM2);
   M1    = M1    - nanmean(M1);
end

Variance = nanstd(LagM2).*nanstd(M1);
Xcorr    = nansum(LagM2.*M1./Variance)./Dof;
%SF       = Variance.*(1-Xcorr);

XC = [Lag.', Chi2.', Dof.', Xcorr.'];
if (nargout>1)
    % calc chi2 minimum
    [Res.MinChi2,Res.MinChi2Ind] = min(Chi2);
    Res.MinChi2Lag = Lag(Res.MinChi2Ind);
    Res.MinDof = Dof(Res.MinChi2Ind);
    
    [Res.MaxCorr, Res.MaxCorrInd] = max(Xcorr);
    Res.MaxCorrLag = Lag(Res.MaxCorrInd);
    
    Res.DeltaChi2    = chi2inv(InPar.ErrConf,Res.Npar);
    NDC2             = numel(Res.DeltaChi2);
    Res.AdjustedChi2 = XC(:,2)+(XC(:,3)-Res.MinDof);
    for Indc2=1:1:NDC2
        %not sure about this
        List0 = Util.find.find_local_zeros(Lag.', (Res.AdjustedChi2 - Res.MinChi2)-Res.DeltaChi2(Indc2) );
        if (~isempty(List0))
            Res.CI{Indc2} = List0(:,1);
        end
    end
    % find local minimum near index minimum
    %Util.find.find_local_extramum(
end


