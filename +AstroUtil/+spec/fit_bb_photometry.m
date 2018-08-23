function [Chi2n,Res]=fit_bb_photometry(Phot,varargin)
% SHORT DESCRIPTION HERE
% Package: AstroUtil
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Chi2=AstroUtil.spec.fit_bb_photometry(tRes)
% Reliable: 
%--------------------------------------------------------------------------


DefV.Dist                 = 40e6;  % [pc]
DefV.VecT                 = logspace(log10(1000),log10(40000),200)';
DefV.VecRad               = logspace(14,16,150)';
DefV.Ebv                  = 0.11;
DefV.ErrSigma             = 1;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% Phot contains:
% Phot().Band
% Phot().System
% Phot().Mag
% Phot().MagErr



Nt        = numel(InPar.VecT);
Nr        = numel(InPar.VecRad);
Nphot     = numel(Phot);
Dof       = Nphot - 2;

ObsMag    = zeros(Nphot,1);
ObsMagErr = zeros(Nphot,1);
Mag       = zeros(Nphot,Nt);

for I=1:1:Nphot
    ObsMag(I)    = Phot(I).Mag;
    ObsMagErr(I) = Phot(I).MagErr;
    [Mag(I,:)] = AstroUtil.spec.blackbody_mag_c(InPar.VecT,Phot(I).Family,Phot(I).Band,Phot(I).System,InPar.VecRad(1),InPar.Dist,InPar.Ebv).';

end

Chi2 = zeros(Nt,Nr);
for Ir=1:1:Nr
    Offset = -2.5.*log10((InPar.VecRad(Ir)./InPar.VecRad(1)).^2);
    Chi2(:,Ir) = sum(((Mag+Offset - ObsMag)./ObsMagErr.^2).^2);
end


MinChi2 = min(Chi2(:));
% renormalize errors
ErrNorm = sqrt(sqrt(MinChi2./Dof));
ObsMagErr1 = ObsMagErr.*ErrNorm;

% re-calculate \chi^2 after error normalization
Chi2n = zeros(Nt,Nr);
for Ir=1:1:Nr
    Offset = -2.5.*log10((InPar.VecRad(Ir)./InPar.VecRad(1)).^2);
    Chi2n(:,Ir) = sum(((Mag+Offset - ObsMag)./ObsMagErr1.^2).^2);
end

MinChi2n = min(Chi2n(:));



%%

Likel = exp(-0.5.*Chi2n);
%Likel = Likel./sum(Likel(:));
%Util.math.int2d(InPar.VecT,InPar.VecRad,Likel');

% unit area in T/r surface
DT = diff(InPar.VecT);
DR = diff(InPar.VecRad);
DT = [DT(1);DT];
DR = [DR(1);DR];
Area = (DT.*DR.');

Likel = exp(-0.5.*Chi2n)./Area;

% Luminosity and its errors
Lum = 4.*pi.*constant.sigma.*InPar.VecT.^4.*(InPar.VecRad.^2).';

LL = [Lum(:),Likel(:)];
LL = sortrows(LL,1);
CumPL = cumsum(LL(:,2))./sum(LL(:,2));
CumPL = CumPL + 5.*eps.*(1:1:numel(CumPL))';
Llow  = interp1(CumPL,LL(:,1),normcdf(-InPar.ErrSigma,0,1));
Lhigh = interp1(CumPL,LL(:,1),normcdf(+InPar.ErrSigma,0,1));
Lmid  = interp1(CumPL,LL(:,1), 0.5);


% Temperature and its errors
LikeT = sum(Likel,2);
CumPT = cumsum(LikeT)./sum(LikeT);
CumPT = CumPT + 5.*eps.*(1:1:numel(CumPT))';
Tlow  = interp1(CumPT,InPar.VecT,normcdf(-InPar.ErrSigma,0,1));
Thigh = interp1(CumPT,InPar.VecT,normcdf(+InPar.ErrSigma,0,1));
Tmid  = interp1(CumPT,InPar.VecT, 0.5);

% Radius and its errors
LikeR = sum(Likel,1).';
CumPR = cumsum(LikeR)./sum(LikeR);
CumPR = CumPR + 5.*eps.*(1:1:numel(CumPR))';
Rlow  = interp1(CumPR,InPar.VecRad,normcdf(-InPar.ErrSigma,0,1));
Rhigh = interp1(CumPR,InPar.VecRad,normcdf(+InPar.ErrSigma,0,1));
Rmid  = interp1(CumPR,InPar.VecRad, 0.5);

Res.MinChi2  = MinChi2;
Res.Dof      = Dof;
Res.MinChi2n = MinChi2n; 
Res.ErrNorm  = ErrNorm;
Res.VecT     = InPar.VecT;
Res.VecRad   = InPar.VecRad;
Res.L        = Lmid;
Res.Lerr     = [Lmid-Llow, Lhigh-Lmid];
Res.T        = Tmid;
Res.Terr     = [Tmid-Tlow, Thigh-Tmid];
Res.R        = Rmid;
Res.Rerr     = [Rmid-Rlow, Rhigh-Rmid];





