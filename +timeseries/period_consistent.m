function [PS]=period_consistent(TS,FreqVec,SubMean)
% SHORT DESCRIPTION HERE
% Package: timeseries
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: N=1000; T = sort(rand(N,1).*1000); Per=23.2; M=sin(2.*pi.*T./Per)+randn(N,1).*0.1;
%          p = T./30 - floor(T./30); F=p<0.3; T=T(F); M=M(F);
%          Freq = (0:0.0005:0.1)';
%          PSn=timeseries.period_norm([T,M],Freq);
% Reliable: 
%--------------------------------------------------------------------------

if (nargin<3)
    SubMean = true;
end


Col.T = 1;
Col.M = 2;

T = TS(:,Col.T);
M = TS(:,Col.M);
if (SubMean)
    M = M - mean(M);
end

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);


V = var(M);

Nf = numel(FreqVec);
WF = zeros(Nf,1);
SinF_hat = zeros(Nf,1);
Mf_hat   = zeros(Nf,1);
for If=1:1:Nf
    % for each frequency
    SinF     = sin(2.*pi.*FreqVec(If).*T); % + sin(4.*pi.*FreqVec(If).*T);
    SinF_hat(If) = sum(SinF.*exp(-2.*pi.*1i.*FreqVec(If).*T));
    %CosF     = sin(2.*pi.*FreqVec(If).*T);
    %CosF_hat(If) = sum(CosF.*exp(-2.*pi.*1i.*FreqVec(If).*T));
    
    WF(If)   = sum(exp(-2.*pi.*1i.*FreqVec(If).*T));
end

for If=1:1:Nf
    % for each frequency
    Mf_hat(If)   = sum(M.*exp(-2.*pi.*1i.*FreqVec(If).*T));
    
    %WinSinConv = conv(SinF_hat,WF,'same');
end
WinSinConv = conv(abs(SinF_hat),abs(WF),'same');
WinSinConv = conv((SinF_hat),(WF),'same');
MfSinConv  = conv(Mf_hat,SinF_hat,'same');

    
%????
PS = Mf_hat.*SinF_hat./(WinSinConv);
PS = abs(PS).^2./V;

PS = [FreqVec(:), PS];


%%
N=1000; T = sort(rand(N,1).*1000); Per=23.2; M=sin(2.*pi.*T./Per)+randn(N,1).*0.1;
p = T./30 - floor(T./30); F=p<0.3; T=T(F); M=M(F);
Freq = (0:0.0005:1)';

Tf = (1:1:N)';
Mf = sin(2.*pi.*Tf./Per)+randn(N,1).*0.1;

Nf   = numel(Freq);
Nt   = numel(T);

WF = zeros(Nf,1);
PS = zeros(Nf,1);
for If=1:1:Nf
    % for each frequency
    %SinF     = sin(2.*pi.*FreqVec(If).*T); % + sin(4.*pi.*FreqVec(If).*T);
    %SinF_hat(If) = sum(SinF.*exp(-2.*pi.*1i.*FreqVec(If).*T));
    %CosF     = sin(2.*pi.*FreqVec(If).*T);
    %CosF_hat(If) = sum(CosF.*exp(-2.*pi.*1i.*FreqVec(If).*T));
    
    WF(If)   = sum(exp(-2.*pi.*1i.*Freq(If).*T));
    PS(If)   = sum(M.*exp(-2.*pi.*1i.*Freq(If).*T));
    
    WFf(If)  = sum(exp(-2.*pi.*1i.*Freq(If).*Tf));
    PSf(If)  = sum(Mf.*exp(-2.*pi.*1i.*Freq(If).*Tf));
end



WF = WF./Nt;
PS = PS./Nt;

WF = WF./sum(abs(WF));

%PSw = timeseries.period([T,M],Freq, 'Type','Win');
%PSn = timeseries.period([T,M],Freq, 'Type','Norm');

[MatX,MatY] = meshgrid((1:1:Nf),(1:1:Nf));
MatInd = MatX - MatY;


MatW = WF(abs(MatInd)+1);

[U,S,V] = svd(MatW);

TT = V* (diag(1./diag(S)))* (U' * PS)



