function [PS]=period_np(Data,FreqVec,Norm)
% SHORT DESCRIPTION HERE
% Package: timeseries
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


Def.Norm = 'Var';
if (nargin==2)
   Norm = Def.Norm;
end

Col.T = 1;
Col.M = 2;
T       = Data(:,Col.T);
N       = numel(T);
Nf      = numel(FreqVec);

M       = Data(:,Col.M) - mean(Data(:,Col.M));
PS      = zeros(Nf,2);
PS(:,1) = FreqVec;
for FreqInd=1:1:Nf
   %Tau           = atan(sum(sin(4.*pi.*FreqVec(FreqInd).*T))./sum(cos(4.*pi.*FreqVec(FreqInd).*T)))./(4.*pi.*FreqVec(FreqInd));
   
   %Tmp = 2.*pi.*FreqVec(FreqInd).*(T-Tau);
   %PS(FreqInd,2) = 0.5.*( sum(M.*cos(Tmp) ).^2./sum( cos(Tmp  ).^2 )  + ...
   %                       sum(M.*sin(Tmp) ).^2./sum( sin(Tmp  ).^2 )    );
                      
  

   PS(FreqInd,2) = abs(sum(M.*exp(-2.*pi.*1i.*T.*FreqVec(FreqInd)) - 0.5.*exp(-4.*pi.*1i.*T.*FreqVec(FreqInd)) )).^2./N;
   
   PS(FreqInd,2) =  abs((sum(M.*exp(-2.*pi.*1i.*T.*FreqVec(FreqInd))  ) ./ sum(exp(-4.*pi.*1i.*T.*FreqVec(FreqInd))  ) )).^2./N;
   %PS(FreqInd,2) = PS(FreqInd,2)./abs(sum(exp(-2.*pi.*1i.*T.*FreqVec(FreqInd)) - 0.5.*exp(-4.*pi.*1i.*T.*FreqVec(FreqInd)) )).^2;
end

switch lower(Norm)
 case 'amp'
    % do nothing
 case 'var'
    PS(:,2) = PS(:,2)./var(M);
 otherwise
    error('Unknwon normalization option');
end
