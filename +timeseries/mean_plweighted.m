function [Mean,MeanSimple]=mean_plweighted(MatF,varargin)
% Power-law weighted mean of a time series.
% Package: timeseries
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: timeseries.mean_plweighted(MatF)
% Reliable: 
%--------------------------------------------------------------------------


%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (nargin==0)
   PL = 2;
   MatF=zeros(101,1000);
    for I=1:1:1000
        Tmp=Util.stat.rand_ps((0:1:100)',[PL 1],randn(101,1).*0.1);
        MatF(:,I)=Tmp(:,2);
    end
end

FT     = fft(MatF,[],1);
MeanFT = abs(mean(FT,2));

Mean = sum(ifft(FT.*MeanFT./sum(MeanFT)));
%Mean = sum(ifft(FT.*MeanFT))./sum((MeanFT));
MeanSimple = mean(MatF);






