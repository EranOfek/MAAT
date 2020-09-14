function TimeSeries=rand_ps(Time,PowSpec,Errors)
% Generate a random time series with a given power spectrum.
% Package: Util.stat
% Description: Generate a random realization of a time series with a given
%              power spectrum (e.g., power-law) and optional gaussian
%              measuments errors.
% Input  : - Column vector of times in which to generate the time series.
%          - Power spectra.
%            * If two element vector is given, then the first is taken as
%            the power-law index of a power-law spectra (i.e., P=(1/f)^-Beta),
%            and the second is the power-spectrum normalization (amplitude^2)
%            at the minimum frequency (total time-span).
%            If the power-spectrum normalization is negative then
%            the standard deviation of the generated LC is taken
%            to be the absolute value of this normalization.
%            * In case that two column matrix is given, then the first
%            column is the frequency and the second column is the power
%            (i.e., amplitude^2), (e.g., [Freq, Power]).
%            Note that the power is taken as zero, outside the
%            frequency range given in this matrix.
%          - Optional column vector of errors (normal) to apply
%            to the time series, defualt is zeros.
% Output : - Time series, [Time, Value, Err].
%            Note that the error is optional.
% Reference: Timmer, J., Konig, M., 1995, A&A 300, 707.
% Tested : Matlab 5.3
%     By : Eran Ofek                      Nov 2002
%    URL : http://weizmann.ac.il/home/eran/matlab/
% Example: TimeSeries=Util.stat.rand_ps((0:1:100)',[2 1],randn(101,1).*0.1);
% Reliable: 2
%-----------------------------------------------------------------------------
FreqOverSampling = 2;                % frequency over sampling
OverRideNyq      = 'y';
Nt               = length(Time);

if (nargin==2)
   Errors = zeros(Nt,1);
elseif (nargin==3)
   % do nothing
else
   error('Illegal number of input arguments');
end

Time = sort(Time);

if (numel(PowSpec)==2)
   %----------------------------------------
   %--- assume a power-law power spectra ---
   %----------------------------------------
   Beta = PowSpec(1);    % power-law index (1/f)^Beta
   Norm = PowSpec(2);    % power-law normalization at lowest frequency

   if (Norm<0)
      NormStD = abs(Norm);
      Norm    = 1;
   else
      NormStD = 1;
      %Norm    = Norm; 
   end

   %--- estimate the Nyquist frequency ---
   MinDT    = min(diff(Time));
   switch OverRideNyq
    case 'n'
       NyqFreq  = 1./MinDT;
    case 'y'
       NyqFreq  = 1;
    otherwise
       error('Unknown OverRideNyq option');
   end
   %--- Time span (lowest) frequency ---
   MaxDT    = Time(end) - Time(1);
   LowFreq  = 1./MaxDT;
   FreqSamp = LowFreq./FreqOverSampling;   % frequency sampling
   FreqCorrNorm = Norm./((2.*pi.*LowFreq).^(-Beta)); % frequency corrected normalization

   % go over all frequencies
   FreqVec  = (LowFreq:FreqSamp:NyqFreq).';

   Nf       = length(FreqVec);
   Rand     = randn(Nf,2);     % generate gaussian random numbers in pairs
   % random number amplitude
   RandAmp  = sqrt(0.5.*FreqCorrNorm.*(2.*pi.*FreqVec).^(-Beta)); 
   % create the Fourier component
   FourComp = RandAmp.*Rand(:,1) + 1i.*RandAmp.*Rand(:,2);
   % negative frequencies f(-w) = f*(w) = conj(f(w))

   % inverse fourier transform (not evenly spaced)
   % slow / old
   %Value = zeros(Nt,1);
   %for I=1:1:Nt,
   %   Value(I) = sum(FourComp.*exp(+2.*pi.*1i.*Time(I).*FreqVec)) + ...,
   %              sum(conj(FourComp).*exp(-2.*pi.*1i.*Time(I).*FreqVec));
   %end

   % fast / new
   Tmp = bsxfun(@times,Time',FreqVec);
   Value = (sum(bsxfun(@times,FourComp,exp(+2.*pi.*1i.*Tmp))) + sum(bsxfun(@times,conj(FourComp),exp(-2.*pi.*1i.*Tmp)))).';
   
   Value = Value .* NormStD./std(Value);

else
   %------------------------------------
   %--- assume general power spectra ---
   %------------------------------------
   FreqVec = PowSpec(:,1);
   PowVec  = PowSpec(:,2);

   %--- estimate the Nyquist frequency ---
%   MinDT    = min(diff(Time));
%   NyqFreq  = 1./MinDT;
%   %--- Time span (lowest) frequency ---
%   MaxDT    = Time(end) - Time(1);
%   LowFreq  = 1./MaxDT;
%   FreqSamp = LowFreq./FreqOverSampling;   % frequency sampling
%   FreqCorrNorm = Norm./((2.*pi.*LowFreq).^(-Beta)); % frequency corrected normalization


   % go over all frequencies
   Nf       = length(FreqVec);
   Rand     = randn(Nf,2);     % generate gaussian random numbers in pairs
   % random number amplitude
   RandAmp  = sqrt(0.5.*PowVec); 
   % create the Fourier component
   FourComp = RandAmp.*Rand(:,1) + 1i.*RandAmp.*Rand(:,2);
   % negative frequencies f(-w) = f*(w) = conj(f(w))

   % inverse fourier transform (not evenly spaced)
   Value = zeros(Nt,1);
   for I=1:1:Nt
      Value(I) = sum(FourComp.*exp(+2.*pi.*1i.*Time(I).*FreqVec)) + ...,
                 sum(conj(FourComp).*exp(-2.*pi.*1i.*Time(I).*FreqVec));
   end



end

%--- add normaly distributed errors ---
ErrorsNoise = Errors.*randn(Nt,1);
Value       = Value + ErrorsNoise;

TimeSeries  = [Time, Value, Errors];
