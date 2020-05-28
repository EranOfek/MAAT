function [Res]=fit_proper_motion(T,RA,Dec,ErrRA,ErrDec,Epoch)
% SHORT DESCRIPTION HERE
% Package: celestial
% Description: 
% Input  : - Time (days)
%          - RA [rad]
%          - Dec [rad];
%          - ErrRA ["]
%          - ErrDec ["]
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res]=celestial.coo.fit_proper_motion(T,RA,Dec,ErrRA,ErrDec,Epoch)
% Reliable: 
%--------------------------------------------------------------------------

RAD   = 180./pi;
JYEAR = 365.25;

if nargin<6
    Epoch = [];
    if nargin<5
        ErrDec = 1;
        if nargin<4
            ErrRA = 1;
        end
    end
end
            
if isempty(Epoch)
    Epoch = median(T);
end


Nt = numel(T);
if numel(ErrRA)==1
    ErrRA = ErrRA.*ones(Nt,1);
end
if numel(ErrDec)==1
    ErrDec = ErrDec.*ones(Nt,1);
end

% deal with zero crossing
if range(RA)>pi
    RA(RA<pi) = 2.*pi + RA(RA<pi);
end

MeanRA  = mean(RA);
MeanDec = mean(Dec);

% normalization
T   = (T - Epoch)./JYEAR;
RA  = RA - MeanRA;
Dec = Dec - MeanDec;
RA  = RA.*RAD.*3600;
Dec = Dec.*RAD.*3600;


H = [ones(Nt,1), T];
[AlphaPar, AlphaParErr] = lscov(H,RA,1./(ErrRA.^2));
[DeltaPar, DeltaParErr] = lscov(H,Dec,1./(ErrDec.^2));

ResidRA  = RA  - H*AlphaPar;
ResidDec = Dec - H*DeltaPar;


Res.Epoch    = Epoch;
Res.rmsRA    = std(ResidRA);   % arcsec
Res.rmsDec   = std(ResidDec);  % arcsec
Res.Chi2_RA  = sum((ResidRA./ErrRA).^2);
Res.Chi2_Dec = sum((ResidDec./ErrDec).^2);
Res.Nobs     = Nt;
Res.Alpha    = AlphaPar(1)./(RAD.*3600) + MeanRA;
Res.Delta    = DeltaPar(1)./(RAD.*3600) + MeanDec;
Res.AlphaErr = AlphaParErr(1).*1000;  % mas
Res.DeltaErr = DeltaParErr(1).*1000;  % mas
Res.PM_Alpha = AlphaPar(2).*1000  .*cos(Res.Delta);   % [mas/yr]
Res.PM_Delta = DeltaPar(2).*1000;   % [mas/yr]
Res.PM_AlphaErr = AlphaParErr(2).*1000  .*cos(Res.Delta);   % [mas/yr]
Res.PM_DeltaErr = DeltaParErr(2).*1000;   % [mas/yr]

% H0 - no PM fit
H = [ones(Nt,1)];
[AlphaPar, AlphaParErr] = lscov(H,RA,1./(ErrRA.^2));
[DeltaPar, DeltaParErr] = lscov(H,Dec,1./(ErrDec.^2));

ResidRA  = RA  - H*AlphaPar;
ResidDec = Dec - H*DeltaPar;
Res.H0_Chi2_RA  = sum((ResidRA./ErrRA).^2);
Res.H0_Chi2_Dec = sum((ResidDec./ErrDec).^2);

% close to 1 means that PM is real...
Res.ProbH1H0_RA = chi2cdf(Res.H0_Chi2_RA - Res.Chi2_RA,1);
Res.ProbH1H0_Dec = chi2cdf(Res.H0_Chi2_Dec - Res.Chi2_Dec,1);



