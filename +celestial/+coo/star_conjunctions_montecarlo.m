function QuantD=star_conjunctions_montecarlo(varargin)
% SHORT DESCRIPTION HERE
% Package: celestial
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: QuantD=celestial.coo.star_conjunctions_montecarlo
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.Prob                  = normcdf([-1 1],0,1);
DefV.Nsim                  = 1000;
DefV.Epoch1                = 2015.5;
DefV.Epoch2                = 2015.5;
DefV.EpochOut              = (2000.0:0.01:2050.0)';
DefV.RA1                   = 0;
DefV.Dec1                  = 1./(RAD.*3600);
DefV.PMRA1                 = 0;
DefV.PMDec1                = -100;
DefV.Plx1                  = 0;
DefV.Covar1                = diag(ones(5,1));
DefV.Covar1(1,1) = 0.1./(RAD.*3600);
DefV.Covar1(2,2) = 0.1./(RAD.*3600);

DefV.RA2                   = 0;
DefV.Dec2                  = 0;
DefV.PMRA2                 = 0;
DefV.PMDec2                = 0;
DefV.Plx2                  = 0;
DefV.Covar2                = diag(ones(5,1));
DefV.Covar2(1,1) = 0.1./(RAD.*3600);
DefV.Covar2(2,2) = 0.1./(RAD.*3600);


InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Pos1 = [InPar.RA1, InPar.Dec1, InPar.PMRA1, InPar.PMDec1, InPar.Plx1].';
Pos2 = [InPar.RA2, InPar.Dec2, InPar.PMRA2, InPar.PMDec2, InPar.Plx2].';


% position realizations:
% row per realization
PosReal1 = mvnrnd(Pos1,InPar.Covar1,InPar.Nsim);
PosReal2 = mvnrnd(Pos2,InPar.Covar2,InPar.Nsim);

Nep = numel(InPar.EpochOut);
D = zeros(Nep,InPar.Nsim);
for Isim=1:1:InPar.Nsim
    [RA1, Dec1] = celestial.coo.proper_motion(InPar.EpochOut,InPar.Epoch1,InPar.Epoch1,...
                                PosReal1(Isim,1),PosReal1(Isim,2),PosReal1(Isim,3),PosReal1(Isim,4),PosReal1(Isim,5),0);
    [RA2, Dec2] = celestial.coo.proper_motion(InPar.EpochOut,InPar.Epoch2,InPar.Epoch2,...
                                PosReal2(Isim,1),PosReal2(Isim,2),PosReal2(Isim,3),PosReal2(Isim,4),PosReal2(Isim,5),0);
    D(:,Isim) = celestial.coo.sphere_dist_fast(RA1,Dec1,RA2,Dec2);
end

QuantD = quantile(D,InPar.Prob,2);
