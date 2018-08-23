function [CCFMat]=ccf(Ser1,Ser2,BinSize,Type,Correct,MaxDeltaT)
%------------------------------------------------------------------------------
% ccf function                                                      timeseries
% Description: Discrete cross-correlation function for two sets of
%              unequaly spaced stationary time series.
%              ccf.m is obsolte - use dcf.m instead.
% Input  : - first series matrix:
%            [Time, Mag, Err], the third column (Err) is optional,
%            if not given, assumes Err=0.
%          - second series matrix:
%            [Time, Mag, Err], the third column (Err) is optional,
%            if not given, assumes Err=0.
%          - BinSize is the delta lag on which the CCF is calculated.
%          - Type of correlation:
%            'normal' : normal correlation, error given by:
%                       (1-rho^2)/sqrt(N-1), (default).
%            'z'      : z-transform, rho is transformed to z
%                       (using z-transform), the errors are given
%                       by 1/sqrt(N-3) has more Gaussian shape.
%                       (Barlow 1989, p. 80)
%          - Correct structure function to measurments errors
%            {'y' | 'n'}, default is 'y'.
%          - Maximum lag to check, default is no limit (i.e., NaN).
% Output : - CCF matrix:
%            [Mid. Lag, CCF, Err, point in bin-lag, Mean Lag, Struct_fun, err_struct_fun].
%            Note that: The structure function has units of amplitude^2.
%                       The structure function is corrected for the
%                       measurments errors, so it has zero amplitude
%                       at zero lag.
% Reference : Edelson, R.A. & Krolik, J.H. 1988 MNRAS 333, 646-659.
%             Koen, C. & Lombard, F. 1993 MNRAS 263, 287-308.
% See Also  : ccf_o.m (old version - for equally spaced...)
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      June 1999
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------------------
MinEl = 5;

if (nargin==3),
   Type      = 'normal';
   Correct   = 'y';
   MaxDeltaT = 'NaN';
elseif (nargin==4),
   Correct   = 'y';
   MaxDeltaT = 'NaN';
elseif (nargin==5),
   MaxDeltaT = 'NaN';
elseif (nargin==6),
   % do nothing
else
   error('Illigal number of input arguments');
end

if (length(Ser1(:,1))<MinEl || length(Ser2(:,1))<MinEl),
   % return empty - do not calculate ccf
   CCFMat = [];
else
   N1 = size(Ser1,1);
   N2 = size(Ser2,1);
   
   T1 = Ser1(:,1);
   T2 = Ser2(:,1);
   M1 = Ser1(:,2) - mean(Ser1(:,2));
   M2 = Ser2(:,2) - mean(Ser2(:,2));
   if (size(Ser1,2)>=3 && size(Ser2,2)>=3),
      % Take errors into account
      E1 = Ser1(:,3);
      E2 = Ser2(:,3);
   
      TakeError = 1;
   else
      E1 = zeros(N1,1);
      E2 = zeros(N2,1);
      TakeError = 0;
   end
   
   
   
   StD1  = std(M1);
   StD2  = std(M2);
   Norm  = sqrt((StD1.^2 - mean(E1).^2).*(StD2.^2 - mean(E2).^2));
   
   N1 = length(T1);
   N2 = length(T2);
   
   MaxTimeSpan = max(T2)-min(T1);
   
   TimeBoundry = [[-rot90([BinSize:BinSize:MaxTimeSpan].',2)]; [0:BinSize:MaxTimeSpan].'];
   TimeLag     = TimeBoundry(1:end-1) + diff(TimeBoundry);
   N_Cell  = length(TimeLag);
   
   M          = zeros(N_Cell,1);
   LagOffset  = zeros(N_Cell,1);
   DCF        = zeros(N_Cell,1);
   UDCF       = zeros(N_Cell,1);
   UDCF2      = zeros(N_Cell,1);
   %CCFMat     = zeros(N_Cell,4);
   StrFun     = zeros(N_Cell,1);
   StrFunErr2 = zeros(N_Cell,1);
   
   
   for I=1:1:N1,
      for J=1:1:N2,
         DeltaT = T2(J) - T1(I);
         if (abs(DeltaT)<=MaxDeltaT),
            [MinVal, TimeLagInd]   = min(abs(DeltaT - TimeLag));
            LagOffset(TimeLagInd)  = LagOffset(TimeLagInd) + MinVal;
            M(TimeLagInd)          = M(TimeLagInd) + 1;
            UDCF(TimeLagInd)       = UDCF(TimeLagInd) + M1(I).*M2(J)./Norm;
            UDCF2(TimeLagInd)      = UDCF2(TimeLagInd) + (M1(I).*M2(J)./Norm).^2;
	        StrFun(TimeLagInd)     = StrFun(TimeLagInd) + (M1(I) - M2(J)).^2;
	        StrFunErr2(TimeLagInd) = StrFunErr2(TimeLagInd) + (E1(I).*2.*(M1(I)-M2(J))).^2 + (E2(J).*2.*(M1(I)-M2(J))).^2;
         end                       
      end
   end
   
   
   LagOffset   = LagOffset./M;
   DCF         = UDCF./M;
   StrFun      = StrFun./M; 
   StrFunErr2  = sqrt(StrFunErr2./(M.^2)); 

   % correct structure function, so it will start with zero amplitude (at zero lag).
   
   switch Correct
    case 'y'
       StrFunShift = (mean(E1).^2 + mean(E2).^2);
       StrFun      = StrFun - StrFunShift;
    case 'n'
       % do nothing
    otherwise
       error('Unknown Correct option');
   end


   % Edelson & Krolik formula:
   %ErrDCF = sqrt(UDCF2 - 2.*UDCF.*DCF + DCF.^2)./(M-1);
   
   % correct for abs(rho)>1
   I = find(abs(DCF)>1);
   DCFe = DCF;
   DCFe(I) = 1;
   
   switch Type
    case 'normal'
       ErrDCF = (1-DCFe.^2)./sqrt(M-1);
    case 'z'
       % z-transform
       DCF    = 0.5.*log((1+DCF)./(1-DCF));
       %ErrDCF = 1./sqrt(M-3);
       ErrDCF = 1./sqrt(sqrt(M)-3);
    otherwise
       error('Unknown CCF type');
   end
   
   CCFMat = [TimeLag, DCF, ErrDCF, M, TimeLag+LagOffset, StrFun, StrFunErr2];
   
   
   
end
