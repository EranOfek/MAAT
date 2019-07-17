function P=period_mi(Data,Freq,varargin)
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
% Example: P=timeseries.period_mi([x,y],Freq);
% Reliable: 
%--------------------------------------------------------------------------


DefV.ColT                 = 1;
DefV.ColM                 = 2;
DefV.NbinX                = 10;
DefV.NbinY                = 10;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

T = Data(:,InPar.ColT);
M = Data(:,InPar.ColM);
XLim = [0 1];
YLim = [min(M)-10.*eps, max(M)+10.*eps];

Nfreq = numel(Freq);
P     = [Freq(:),nan(Nfreq,1)];
for Ifreq=1:1:Nfreq
    Phase = mod(T,1./Freq(Ifreq));
    
    P(Ifreq,2) = Util.stat.mutual_info(Phase,M,'XLim',XLim,'YLim',YLim,'NbinX',InPar.NbinX,'NbinY',InPar.NbinY);
end