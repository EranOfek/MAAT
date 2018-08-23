function [DCF,SF]=dcf(Ser1,Ser2,BinSize,Abs)
%--------------------------------------------------------------------------
% ccf function                                                  timeseries
% Description: Discrete cross-correlation function and structure
%              function for two sets of unequaly spaced stationary
%              time series.
% Input  : - first series matrix:
%            [Time, Mag, Err], the third column (Err) is optional,
%            if not given, assumes Err=0.
%          - second series matrix:
%            [Time, Mag, Err], the third column (Err) is optional,
%            if not given, assumes Err=0.
%          - BinSize is the delta lag on which the CCF is calculated.
%            Alternatively this can be a colmn vector of lag edges.
%          - Use absolute value of time lag {'y' | 'n'}, default is 'y'.
%            In that case CCF will be symmetrical around zero lag.
% Output : - Structure containing the discrete correlation function
%            with the following fields:
%            .Lag      - Mean lag.
%            .LagEdge  - Lower and upper edges of Lag.
%            .N        - Number of points in each lag.
%            .DCF      - The discrete correlation function.
%            .medDCF   - The median discrete correlation function.
%            .StD      - StD of discrete correlation function.
%                        To estimate error divide .StD by .N
%            .BS_Err   - Bootstrap errors in discrete correlation function.
%          - Structure function structure with the following fields:
%            .Lag      - Mean Lag.
%            .SF       - Structure function [units if input Mag/flux]
%            .Err      - Error in structure function.
% Reference : Edelson & Krolik, 1988, MNRAS 333, 646-659.
%             Koen, C. & Lombard, F. 1993 MNRAS 263, 287-308.
% See Also  : ccf.m, ccf_o.m (old version - for equally spaced data)
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Sep 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------
MeanType = 'mean';
Nsim     = 100;

Def.Abs = 'y';
if (nargin==3),
   Abs = Def.Abs;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end
   

N1 = size(Ser1,1);
N2 = size(Ser2,1);

if (size(Ser1,2)==2),
  Ser1 = [Ser1, zeros(N1,1)];
end
if (size(Ser2,2)==2),
  Ser2 = [Ser2, zeros(N2,1)];
end

Col.T = 1;
Col.M = 2;
Col.E = 3;
T1 = Ser1(:,Col.T);
T2 = Ser2(:,Col.T);
M1 = Ser1(:,Col.M);
M2 = Ser2(:,Col.M);
E1 = Ser1(:,Col.E);
E2 = Ser2(:,Col.E);

switch lower(MeanType)
 case 'mean'
    Mean1 = nanmean(M1);
    Mean2 = nanmean(M2);
 case 'median'
    Mean1 = nanmedian(M1);
    Mean2 = nanmedian(M2);
 otherwise
    error('Unknown MeanType option');
end
Std1 = nanstd(M1);
Std2 = nanstd(M2);
CorrStd = sqrt( (Std1.^2 - nanmedian(E1).^2).*(Std2.^2 - nanmedian(E2).^2) );
%CorrStd = sqrt( (Std1.^2 ).*(Std2.^2 ) );

FunUDCF = @(M,Mean1,Mean2,CorrStd) (M(:,1)-Mean1).*(M(:,2)-Mean2)./CorrStd;
FunDCF  = @(M,Mean1,Mean2,CorrStd) nanmean((M(:,1)-Mean1).*(M(:,2)-Mean2)./CorrStd);
FunDCFm = @(M,Mean1,Mean2,CorrStd) nanmedian((M(:,1)-Mean1).*(M(:,2)-Mean2)./CorrStd);



if (length(BinSize)==1),
   % uvenly spaced tau
   MaxLag = max(abs(max(T1)-min(T2)),abs(max(T2)-min(T1)));

   LagEdgeVec = (BinSize.*0.5:BinSize:MaxLag).';
   LagEdgeVec = [-flipud(LagEdgeVec); LagEdgeVec];
else
   LagEdgeVec = BinSize;
end
Nlag = length(LagEdgeVec) - 1;

DCF.Lag    = (LagEdgeVec(1:end-1) + LagEdgeVec(2:end)).*0.5;
DCF.LagEdge= zeros(Nlag,2).*NaN;
DCF.N      = zeros(Nlag,1).*NaN;
DCF.DCF    = zeros(Nlag,1).*NaN;
DCF.medDCF = zeros(Nlag,1).*NaN;
DCF.StD    = zeros(Nlag,1).*NaN;
DCF.BS_Err = zeros(Nlag,1).*NaN;


% matrix of all combination of time differences between two series:
TimeDiff = repmat(T1,1,N2) - repmat(T2.',N1,1);

for Ilag=1:1:Nlag,
   switch lower(Abs)
    case 'y'
       [I,J] = find(abs(TimeDiff)>abs(LagEdgeVec(Ilag)) & abs(TimeDiff)<=abs(LagEdgeVec(Ilag+1)));
    case 'n'
       [I,J] = find(TimeDiff>LagEdgeVec(Ilag) & TimeDiff<=LagEdgeVec(Ilag+1));
    otherwise
       error('Unknown Abs option');
   end

   DCF.LagEdge(Ilag,1:2) = [LagEdgeVec(Ilag), LagEdgeVec(Ilag+1)];
   DCF.N(Ilag,1)         = length(I);

   if (DCF.N==0),
      DCF.DCF(Ilag)      = NaN;
      DCF.medDCF(Ilag)   = NaN;
      DCF.StD(Ilag)      = NaN;
      DCF.BS_Err(Ilag)   = NaN;
   else
      DCF.DCF(Ilag)         = FunDCF([M1(I),M2(J)],Mean1,Mean2,CorrStd);
      DCF.medDCF(Ilag)      = FunDCFm([M1(I),M2(J)],Mean1,Mean2,CorrStd);

      DCF.StD(Ilag)         = nanstd( FunUDCF([M1(I),M2(J)],Mean1,Mean2,CorrStd) - DCF.DCF(Ilag));
      DCF.BS_Err(Ilag)      = bootstrap_std([M1(I),M2(J)],FunDCF,Nsim,Mean1,Mean2,CorrStd);
   end
end

if (nargout>1),
   % Estimate the structure function:
   % sqrt(2 * Norm * (1- DCF))
   % where Norm = std(Series)
   Var   = Std1.*Std2;   % The variance
   SF.Lag = DCF.Lag;
   SF.SF  = sqrt(2.*Var.*(1-DCF.DCF));
   SF.Err = 2.*Var.*DCF.BS_Err./(2.*sqrt( 2.*Var.*(1-DCF.DCF)));
end
