function XC=dcf1(TS1,TS2,Lag,varargin)
%--------------------------------------------------------------------------
% dcf1 function                                                 timeseries
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: T1=sort(rand(1000,1).*1000); T2=sort(rand(1000,1).*1000);
%          X1=sin(2.*pi.*T1./100); X2=sin(2.*pi.*T2./100+0.5);
%          XC = dcf1([T1,X1],[T2,X2],2);
% Reliable: 
%--------------------------------------------------------------------------

ColT = 1;
ColF = 2;
ColE = 3;


DefV.MeanSub1           = 'mean';
DefV.MeanSub1Par        = 1;
DefV.MeanSub2           = 'mean';
DefV.MeanSub2Par        = 1;
DefV.StdNorm1           = @std;
DefV.StdNorm2           = @std;
DefV.MinMaxLag          = [-100 100];
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

% Subtract mean from 1st series
if (isa(InPar.MeanSub1,'function_handle')),
    Mean1 = InPar.MeanSub1(TS1(:,ColF));
else
    switch lower(InPar.MeanSub1)
        case 'mean'
            Mean1 = mean(TS1(:,ColF));
        case 'median'
            Mean1 = median(TS1(:,ColF));
        case 'poly'
            Par1  = polyfit(TS1(:,ColT),TS1(:,ColF),InPar.MeanSub1Par);
            Mean1 = polyval(Par1,TS1(:,ColT));
        otherwise
            error('Unknown MeanSub1 option');
    end
end

TS1(:,ColF) = TS1(:,ColF) - Mean1;
        

% Subtract mean from 2nd series
if (isa(InPar.MeanSub2,'function_handle')),
    Mean2 = InPar.MeanSub2(TS2(:,ColF));
else
    switch lower(InPar.MeanSub2)
        case 'mean'
            Mean2 = mean(TS2(:,ColF));
        case 'median'
            Mean2 = median(TS2(:,ColF));
        case 'poly'
            Par2  = polyfit(TS2(:,ColT),TS2(:,ColF),InPar.MeanSub2Par);
            Mean2 = polyval(Par2,TS2(:,ColT));
        otherwise
            error('Unknown MeanSu2b option');
    end
end

TS2(:,ColF) = TS2(:,ColF) - Mean2;


% normalize series
TS1(:,ColF) = TS1(:,ColF)./InPar.StdNorm1(TS1(:,ColF));


      
B1 = binning(TS1,Lag,InPar.MinMaxLag(1),InPar.MinMaxLag(2));
B2 = binning(TS2,Lag,InPar.MinMaxLag(1),InPar.MinMaxLag(2));
B1(isnan(B1(:,2)),2) = 0;
B2(isnan(B2(:,2)),2) = 0;

LagVector = (InPar.MinMaxLag(1):Lag:InPar.MinMaxLag(2))';
Nlv       = numel(LagVector);

XC.Lag = LagVector(1:end-1) + 0.5.*Lag;
XC.XC  = ifft(fft(B1(:,2)).*conj(fft(B2(:,2))))./numel(XC.Lag);

    





