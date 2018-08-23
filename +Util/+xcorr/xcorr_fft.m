function [Lag,XC,Info]=xcorr_fft(S1,S2,Back,varargin);
%---------------------------------------------------------------------------
% xcorr_fft function                                                General
% Description: cross correlation of two 1-D serieses using FFT.
%              The shortest vector will be padded by zeros (at the end).
% Input  : - Vector of first series (equaly spaced).
%          - Second vector.
%          - Background subtraction method:
%            'none' - do not subtract background.
%            'mean' - subtract the mean.
%            'median' - subtract the median (default).
%            'medfilt' - subtract a median filter, which length is
%                        specified as the next input argument.
%          * Arbitrary number of input argument to pass to the
%            background subtraction algorithm.
% Output : - Vector of lags.
%          - Vector of cross-correlation between the two vectors,
%            corresponding to the lags.
%          - Structure containing information about te best shift
%            between the vectors. Fields are:
%            .Shift   - Shift of first vector relative to the second vector.
%            .Corr    - Correlation at shift.
%            .BestShift - Interpolated position of best shift.
%            .BestCorr  - Correlation at best shift.
%            .Best2dr   - 2nd derivative at best shift.
% See also: xcorr_fft_multi.m
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    August 2010
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: R1 = rand(1e5,1); R2 = rand(1e5,1);
%          R1(2001:2100)=100; R2(1001:1100)=100;
%          [XC,Lag,Info] = xcorr_fft(R1,R2);
% Reliable: 2
%---------------------------------------------------------------------------
Def.Back = 'median';
if (nargin==2),
   Back = Def.Back;
end


S1 = subtract_back1d(S1,Back,varargin{:});
S2 = subtract_back1d(S2,Back,varargin{:});


N1 = length(S1);
N2 = length(S2);
if (N1>N2),
   % pad S2
   S2 = [S2; zeros(N1-N2,1)];
end
if (N2>N1),
   % pad S1
   S1 = [S1; zeros(N2-N1,1)];
end
N1 = length(S1);
N2 = length(S2);

XC  = fftshift(real(ifft(fft(S1).*conj(fft(S2)))));
Lag = [-floor(N1./2):1:ceil(N1./2-1)]'; 
Factor = 1./(N1.*std(S1).*std(S2));
XC     = XC.*Factor;

%[Max,MaxInd] = max(XC);
Buffer = 6;
[Max,MaxInd] = max(XC(Buffer:end-Buffer));
MaxInd = MaxInd + Buffer - 1;
Info.Shift = Lag(MaxInd);
Info.Corr  = Max;


if (MaxInd<6 || MaxInd>(length(XC)-5) || isempty(find(~isnan(XC),1))),
   Info.BestShift = NaN;
   Info.BestCorr  = NaN;
   Info.Best2dr   = NaN;
else    
   Extram = find_local_extramum(Lag(MaxInd-5:MaxInd+5),XC(MaxInd-5:MaxInd+5));
   if (isempty(Extram)),
       Info.BestShift = NaN;
       Info.BestCorr  = NaN;
       Info.Best2dr   = NaN;
   else
      Imax = find(Extram(:,3)<0);
      if (isempty(Imax)),
          Info.BestShift = NaN;
          Info.BestCorr  = NaN;
          Info.Best2dr   = NaN;
      else
          [MaxE,MaxEInd] = max(Extram(Imax,2));
          Info.BestShift = Extram(Imax(MaxEInd),1);
          Info.BestCorr  = Extram(Imax(MaxEInd),2);
          Info.Best2dr   = Extram(Imax(MaxEInd),3);
      end
   end
end


