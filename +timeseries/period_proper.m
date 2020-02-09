function []=period_proper(Data,FrqVec,varargin)
% SHORT DESCRIPTION HERE
% Package: timeseries
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------



DefV. = 
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (nargin==0)
    % simulations mode
    
    Err    = 0.1;
    Period = 12.3;
    Time = (0:1:200).';
    Time = timeseries.random_time_sequence;
    Time = Time(Time<500);
    Flux = 1.*sin(2.*pi.*Time./Period) + Err.*randn(size(Time));
    Data = [Time, Flux];
    
    FreqVec = (0:0.001:1).';
    FreqVec = (-1:0.001:1).';
end


% calculate the window function
Nd = size(Data,1);
WhC = timeseries.period_complex([Data(:,1),ones(Nd,1)],FreqVecAll,'amp',false);  
Wh  = timeseries.period_normnl([Data(:,1),ones(Nd,1)],FreqVecAll,'amp',false);


rPS = timeseries.period_normnl([Time,Flux],FreqVec);

Nf = numel(FreqVec);
Nt = numel(Time);
eFtT = exp(2.*pi.*1i.*FreqVec.'.*Time);
eFTt = exp(2.*pi.*1i.*FreqVec.*Time.');


PS = zeros(Nf,1);
tic;
for If=1:1:Nf
    %tic;
    % FT of periodic function with frequency FreqVec(If) sampled at Time
    Freq = FreqVec(If);
    
    PerF = zeros(Nf,1);
    for Ifin=1:1:Nf
        % equivalent to convolution of the periodic function with the
        % window function
        % normalization????
        % note the use of + in the first term...
        PerF(Ifin) = sum(exp(+2.*pi.*1i.*Freq.*Time).*exp(-2.*pi.*1i.*FreqVec(Ifin).*Time))./Nt;
    end
    % with no loops
    %PerF = sum(exp(-2.*pi.*1i.*Freq.*Time).*exp(-2.*pi.*1i.*FreqVec.'.*Time),1).';
  %  PerF = sum(exp(-2.*pi.*1i.*Freq.*Time).*eFtT,1).';
    
  
    % instead:
    % cross correlate PerF (in frequency domain) with the observed power
    % spectrum...
    PS(If) = sum(abs(PerF).*rPS(:,2));
end
  
  
    % take the inverse FT of PSfin
     PerT = zeros(Nt,1);
     for It=1:1:Nt
         PerT(It) = sum(PerF.*exp(2.*pi.*1i.*FreqVec.*Time(It)))./Nt;
     end
    % with no loops
    %PerT1 = sum(PerF.*exp(2.*pi.*1i.*FreqVec.*Time.'));
%    PerT = sum(PerF.*eFTt).';
    
    % calculate the power spectrum with the windowed periodic function
    PS(If) = sum(Flux.*PerT)./Nt;
    
    %toc
    
end
    
toc






for If=1:1:Nf
    e.^(-1i.*2.*pi.*FrqVecAll(If).*Data(:,1));
end


T = Data(:,1);
M = Data(:,2);
N = numel(T);
V = var(M);

Nf = numel(FreqVec);
WF = zeros(Nf,1);
SinF_hat = zeros(Nf,1);
Mf_hat   = zeros(Nf,1);

Nf = numel(FreqVec);
for If=1:1:Nf
    % for each frequency
    SinF     = sin(2.*pi.*FreqVec(If).*T); % + sin(4.*pi.*FreqVec(If).*T);
    SinF_hat(If) = sum(SinF.*exp(-2.*pi.*1i.*FreqVec(If).*T));
    %CosF     = sin(2.*pi.*FreqVec(If).*T);
    %CosF_hat(If) = sum(CosF.*exp(-2.*pi.*1i.*FreqVec(If).*T));
    
    WF(If)   = sum(exp(-2.*pi.*1i.*FreqVec(If).*T));
    
   
end

PS     = zeros(Nf,1);
ExpF_W = zeros(Nf,1);
C      = zeros(Nf,1);
for If=1:1:Nf
    % periodic function at FreqVec(If)
    ExpF = exp(-2.*pi.*1i.*FreqVec(If).*T);
    % calc the FT of ExpF
    PS_ExpF = zeros(Nf,1);
    for Ifi=1:1:Nf
        PS_ExpF(Ifi) = sum(ExpF.*exp(-2.*pi.*1i.*FreqVec(Ifi).*T))  ./N;
    end
    % convolve PS_ExpF with WF
    Con = conv(PS_ExpF,WF,'same');
    % sample at T 
    C(If) = sum(Con.*exp(2.*pi.*1i.*FreqVec(If).*T));
    
    ExpF_W(If) = C(If);
    
    PS(If) = sum(M.*ExpF_W(If))./N;
end

    



