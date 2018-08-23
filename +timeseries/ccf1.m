function [XC,Res]=ccf1(S1,S2,varargin)
%--------------------------------------------------------------------------
% ccf1 function                                                 timeseries
% Description: 
% Input  : - 
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Sep 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X1=(1:0.1:100).'; Y1=sin(X1.*2.*pi./50)+randn(size(X1)).*0.01;
%          X2=(1:0.1:100).'; Y2=sin((X2+10).*2.*pi./50)+randn(size(X2)).*0.01;
%          [XC,Res]=ccf1([X1,Y1,0.01.*ones(size(X1))],[X2,Y2,0.01.*ones(size(X2))],1);
% Reliable: 
%--------------------------------------------------------------------------

DefV.ColT    = 1;
DefV.ColM    = 2;
DefV.ColE    = 3;
DefV.MeanFun = 'mean';   % {mean|median|none}
DefV.MaxT    = 100;      % 
DefV.StepT   = 1;        % 
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

VecT = (-InPar.MaxT:InPar.StepT:InPar.MaxT).';

N1 = size(S1,1);
N2 = size(S2,1);

if (size(S1,2)<3 || size(S2,2)<3),
    InPar.ColE = [];
end

if (isempty(InPar.ColE)),
    S1 = [S1(:,[InPar.ColT, InPar.ColM]), zeros(N1,1)];
    S2 = [S2(:,[InPar.ColT, InPar.ColM]), zeros(N2,1)];
else
    S1 = S1(:,[InPar.ColT, InPar.ColM, InPar.ColE]);
    S2 = S2(:,[InPar.ColT, InPar.ColM, InPar.ColE]);
end

Range1 = range(S1(:,1));
Range2 = range(S2(:,1));
MinRange = min(Range1,Range2);


% subtract mean
switch lower(InPar.MeanFun)
    case 'none'
        % do nothing
    case 'mean'
        S1(:,2) = S1(:,2) - mean(S1(:,2));
        S2(:,2) = S2(:,2) - mean(S2(:,2));
    case 'median'
        S1(:,2) = S1(:,2) - median(S1(:,2));
        S2(:,2) = S2(:,2) - median(S2(:,2));
    otherwise
        error('Unknown MeanFun option');
end

StdS1 = std(S1(:,2));
StdS2 = std(S2(:,2));



% for each DT
Nt = length(VecT);
XC = zeros(Nt,4);
XC(:,1) = VecT + [(VecT(2:end) - VecT(1:end-1)).*0.5;0];

for I1=1:1:N1,
    for I2=1:1:N2,
        DiffT = S1(I1,1) - S2(I2,1);
        FlagI = find(DiffT>VecT,1,'last'); % & DiffT<=VecT;
        
        XC(FlagI,2) = XC(FlagI,2) + S1(I1,2).*S2(I2,2);
        XC(FlagI,3) = XC(FlagI,3) + (S1(I1,2).*S2(I2,3)).^2 + (S1(I1,3).*S2(I2,2)).^2;
        XC(FlagI,4) = XC(FlagI,4) + 1;  % counter
    end
end
XC(:,2) = XC(:,2)./(StdS1.*StdS2.*XC(:,4));   % normalize
XC(:,3) = sqrt(XC(:,3))./(StdS1.*StdS2.*XC(:,4));

Res.XC = XC;

% for I1=1:1:N1,
%     DiffT = S1(I1,1) - S2(:,1);
%     PopInd = round(DiffT./Lag) + Ind0;  % index to populate
%     Res.Counter(PopInd) = Res.Counter(PopInd) + 1;
%     if (InPar.StdLocal),
%         StdS1 = std(S1(,2));
%     Res.ProdM(PopInd)   = Res.ProdM(PopInd) + S1(I1,2).*S2(:,2)./(StdS1.*StdS2);
% end    

%Res.XCorr = Res.ProdM./Res.Counter;

    
    
    
    
    
    