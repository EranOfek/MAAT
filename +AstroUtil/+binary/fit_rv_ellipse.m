function [Res,FreqVec,BestFit]=fit_rv_ellipse(Time,RV,VecP)
% Fit radial velocity to ellipse as a function of period
% Package: AstroUtil.binary
% Description: Given radial velocity as a function of time and a trial
%              period, fit an ellipse and calculate the RMS as a function
%              of trial period.
% Input  : - Vector of time.
%          - Vector of radial velocity
%          - Vector of trial periods.
% Output : - Structure array of results. For each trial period return
%            the rms of the fit.
%          - Vector of frequencies (i.e., 1/period).
%          - Structure of best fit parameters.
% Reference: Bhattacharyya & Nityanada (2008; MNRAS 387, 273-278)
%     By : Eran O. Ofek            Dec 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res,FreqVec]=AstroUtil.binary.fit_rv_ellipse; % simulation mode
%          plot(FreqVec,[Res.RMS])
% Reliable: 2

import Util.fit.*


if (nargin==0)
    % simulation mode
    Time = (1:1:700)';
    P    = 365;
    T    = 0;
    q    = 1;
    e    = 0.6;
    i    = 1;
    omega = 3.9;
    Nt = numel(Time);

    % generate RV curve
    [RV,K2] = AstroUtil.binary.binary_rv(Time,P,T,q,e,i,omega);
    RV = RV./1e5;
    RV = RV + randn(Nt,1);

    VecP = (36:1:700)';
    %FreqVec = 1./VecP; %(0.001:0.00001:0.01)';
end

FreqVec = 1./VecP;
Np   = numel(FreqVec);

%%
for Ip=1:1:Np
    % for each trial period
    [RV_Even,RV_Odd]=AstroUtil.binary.rv2ellipse(Time,RV,VecP(Ip)); %1./FreqVec(Ip));

    Res(Ip)=Util.fit.fit_ellipse(RV_Odd,RV_Even,[true true false true false]);
    
    %plot(RV_Even,RV_Odd,'.')
    %hold on
    
end

Power = 1./[Res.RMS].';
%StD=stdfilt1(Power,3,[],0.68);
%plot(FreqVec,Power)
%hold on;
%plot(FreqVec,StD)


% In order to determine the sign of omega need to check if the
% ellipse is clockwise or counter clockwise
%...


 [RV_Even,RV_Odd,~,RV_Time]=AstroUtil.binary.rv2ellipse(Time,RV,VecP(Ip)); %1./FreqVec(Ip));

 MeanX = mean(RV_Odd);
 MeanY = mean(RV_Even);
 [~,PA] = Util.Geom.plane_dist(RV_Odd,RV_Even,MeanX,MeanY);
 [~,SI] = sort(RV_Time);
 Direction = median(sign(diff(PA(SI))));
 
    


% find minimum RMS
[~,MinIp] = min([Res.RMS]);

A = Res(MinIp).Par(1);
B = Res(MinIp).Par(2);
C = Res(MinIp).Par(3);

% Note there is an error in the Bhattacharyya & Nityanada (2008;
% MNRAS 387, 273-278) paper below Eq. A3:
% 1-b^2/d^2 should be 1-d^2/b^2 everywhere...
% Solving for: a,b,d:
%  syms A B C a b d
%  SS=solve(A==(1/a^2)/(1-d^2/b^2),B==(1/b^2)/(1-d^2/b^2),C==(-2*d/b^2)/(1-d^2/b^2),a,b,d)

a= ((C.^2 + 4.*B)./(A.*B)).^(1./2)./2;
b= (C.^2 + 4.*B).^(1./2)./(2.*B);
GMh = sqrt(a.^2 + b.^2);
d= -C./(2.*B);
e = d/b;


if (e<0)
    if (Direction>0)
        omega = 2.*pi - atan2(a,-b);
    else
        omega = atan2(a,Direction.*b);
    end
else
    if (Direction<0)
        omega = pi+atan2(a,Direction.*b);
    else
        omega = atan2(a,Direction.*b);
    end
end
e = abs(e);

BestFit.Period = VecP(MinIp);
BestFit.GMh   = GMh;
BestFit.e     = e;
BestFit.omega = omega;


