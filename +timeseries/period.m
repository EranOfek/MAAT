function [PS,Peaks]=period(Data,Freq,varargin)
% Periodicity search in a non-equally spaced time series
% Package: timeseries
% Description: Periodicy search in a time series. This function can be
%              used to run different types of periodicity searches and
%              window functions, and automatically select the range
%              of frequencies to search.
% Input  : - Two column matrix containing the time series
%            [Time, measurment] or [Time, measurment, error].
%          - Frequency range and interval in which to calculate the
%            power spectrum.
%            This is a column vector of frequencies at which to
%            calculate the power spectrum.
%            Alternatively, this is a row vector of the following form
%            [MinFreq, MaxFreq, FreqInterval].
%            In this case MinFreq is the minimum frequency,
%            MaxFreq is the maximum frequency and FreqInterval
%            is the frequency interval. All given in units of 1/time.
%            Or [MaxFreq, FreqInterval]. In this case MinFreq=0.
%            Or [OverSampling] where oversampling is a scalar indicating
%            the over sampling to use for the FreqInterval. In this case
%            FreqInterval = MinFreq/OverSampling.
%            MinFreq      = 1/range(Time).
%            MaxFreq      = 1/min(diff(Time)).
%            If this parameter is an empty matrix (i.e., []) or
%            not specified, then the program use OverSampling=4
%            and calculate the rest of the parameters automatically.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            Where the possible keywords are available:
%            'Type'   - Type of power spectrum and algorithm:
%                       'Norm'     - Normal (Lomb) periodogram.
%                       'Norm2'    - Normal (Lomb) periodogram of the
%                                    second order (two equally weighted
%                                    harmonies).
%                       'complex'  - The Fourier transform of an unevenly
%                                    spaced data.
%                       'Scargle'  - Scargle periodogram (Default).
%                       'NormNL'   - Normal periodogram with
%                                    no loops. The no-loops version is
%                                    faster but requires a lot of memory.
%                       'ScargleNL'- Scargle periodogram with no loops.
%                       'fft'      - FFT for equaly spaced time series.
%                                    In this case the Freq vector is
%                                    overrid by the program.
%                       'Win'      - Calculate the window function.
%                                    Always normalized by amplitude.
%                       'WinNL'    - Calculate the window function using
%                                    no loops.
%                                    Always normalized by amplitude.
%           'Norm'    - Method of normalization:
%                       'var'      - Normalized by the variance (Default).
%                       'amp'      - Amplitude normalization.
%           'RmNaN'   - Remove NaNs from time series beofre analysis
%                       {'y'|'n'}, default is 'y'.
%           'Taper'   - Taper function name {'none'|'cosbell'|'trapz'}.
%                       Default is 'none'.
%           'TaperPar'- Taper parameters - see taper.m for details.
% Output : - Power spectrum [Frequency, Power].
%          - Structure containing information about the peaks in the power
%            spectrum. The following fields are available:
%            .FAB   - 1 - False alarma probability
%            .Freq  - Frequency
%            .Per   - Period (1/Frequency)
%            .PS    - Power.
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    May 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=sort(rand(100,1).*100); Y=sin(X.*2.*pi.*0.9)+randn(100,1);
%          [PS,Peaks]=timeseries.period([X,Y]);
%          [PS,Peaks]=timeseries.period([X,Y],[0 2 0.001]);
%          [PS,Peaks]=timeseries.period([X,Y],[0 2 0.001],'Type','Scargle');
%          [PS,Peaks]=timeseries.period([X,Y],[],'Taper','cosbell','TaperPar',0.5);
% Reliable: 2
%--------------------------------------------------------------------------
import timeseries.*

Col.T  = 1;
Col.X  = 2;
Col.Xe = 3;

Def.Freq = [];
if (nargin==1)
   Freq = Def.Freq;
end

if (size(Freq,1)>size(Freq,2))
   % Assuming Freq is a column vector of frequencies
else
   switch length(Freq)
    case 0
       OverSampling = 4;
       Time         = Data(:,Col.T);
       MinFreq      = 1./range(Time);
       MaxFreq      = 1./min(diff(Time));   
       FreqInterval = MinFreq./OverSampling;
    case 1
       OverSampling = Freq(1);
       Time         = Data(:,Col.T);
       MinFreq      = 1./range(Time);
       MaxFreq      = 1./min(diff(Time));   
       FreqInterval = MinFreq./OverSampling;
    case 2
       MinFreq      = 0;
       MaxFreq      = Freq(1);
       FreqInterval = Freq(2);
    case 3
       MinFreq      = Freq(1);
       MaxFreq      = Freq(2);
       FreqInterval = Freq(3);
    otherwise
       error('Freq have illegal size');
   end
   Freq = (MinFreq:FreqInterval:MaxFreq).';
end

DefV.Type   = 'Scargle';
DefV.Norm   = 'Var';
DefV.RmNaN  = 'y';
DefV.Taper  = 'none';
DefV.TaperPar = [];
Par = InArg.populate_keyval(DefV,varargin,mfilename);

switch lower(Par.Taper)
    case 'none'
        % do nothing
    otherwise
        TaperVal = taper(Data(:,Col.T),'TaperFun',Par.Taper,'TaperPar',Par.TaperPar);
        Data(:,Col.X) = Data(:,Col.X).*TaperVal;
end

switch lower(Par.RmNaN)
    case 'y'
        Data = Data(~isnan(Data(:,2)),:);
    otherwise
        % do nothing
end

%PS = [];
switch lower(Par.Type)
 case 'fft'
    PS = period_fft(Data,Par.Norm);
 case 'norm'
    PS = period_norm(Data,Freq,Par.Norm);
 case 'norm2'
    PS = period_norm_order2(Data,Freq,Par.Norm);
 case 'complex'
    PS = period_complex(Data,Freq,Par.Norm);  
 case 'scargle'
    PS = period_scargle(Data,Freq,Par.Norm);
 case 'normnl'
    PS = period_normnl(Data,Freq,Par.Norm);
 case 'scarglenl'
    error('ScargleNL is not yet available');
 case 'wincomplex'
    Nd = size(Data,1);
    PS = timeseries.period_complex([Data(:,1),ones(Nd,1)],Freq,'amp',false);  
 case 'win'
    Nd = size(Data,1);
    PS = period_norm([Data(:,1),ones(Nd,1)],Freq,'amp',false);
 case 'winnl'
    Nd = size(Data,1);
    PS = period_normnl([Data(:,1),ones(Nd,1)],Freq,'amp',false);
 otherwise
    error('Unknown Type option');
end


if (nargout>1)
   Ipeak = diff(sign(diff([0;PS(:,2);0])))==-2;
   Nd    = size(Data,1);
   Ni = -6.362 + 1.193.*Nd + 0.00098.*Nd.^2;  % number of independent data points : Horne, J.H. & Baliunas, S.L. ApJ 302, 757-763 (1986)

   Peaks.FAB  = (1-exp(-PS(Ipeak,2))).^Ni;
   %[Peaks.FAB, SI] = sort(Peaks.FAB);
   Peaks.Freq = PS(Ipeak,1);
   Peaks.Per  = 1./PS(Ipeak,1);
   Peaks.PS   = PS(Ipeak,2);
   Sorted = sortrows([Peaks.FAB, Peaks.Freq, Peaks.Per, Peaks.PS],4);
   Peaks.FAB  = Sorted(:,1);
   Peaks.Freq = Sorted(:,2);
   Peaks.Per  = Sorted(:,3);
   Peaks.PS   = Sorted(:,4);
  
end


