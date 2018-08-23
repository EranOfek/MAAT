function [Out]=proper_motion_sdss_ps1(RA,Dec,varargin)
% SHORT DESCRIPTION HERE
% Package: VO.search
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Out]=VO.search.proper_motion_sdss_ps1(RA,Dec);
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

JYear = 365.25;

DefV.SearchRadius         = 7;  % [arcsec]
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

%RA=2.8307
%Dec= 0.060672

% PS1 query
[CatPS1,ColCellPS1] = catsHTM.cone_search('PS1',RA,Dec,InPar.SearchRadius);

% SDSS query
[RA1,RA2,Dec1,Dec2] = celestial.coo.coo2box(RA,Dec,InPar.SearchRadius./(3600.*RAD));
Q{1} = {'ra','dec','mjd','modelMag_g','modelMagErr_g','modelMag_r','modelMagErr_r','modelMag_i','modelMagErr_i'};
Q{2} = {'PhotoPrimary'};
Q{3} = sprintf('(ra between %f and %f) and (dec between %f and %f)',RA1.*RAD,RA2.*RAD,Dec1.*RAD,Dec2.*RAD);
[CatSDSS,Msg,FirstLine]=VO.SDSS.run_sdss_sql(Q);


SysErr = 0.015;

PS1.RA          = CatPS1(:,1);
PS1.Dec         = CatPS1(:,2);
PS1.Epoch       = CatPS1(:,5);
PS1.posMeanChi2 = CatPS1(:,6);
PS1.g           = CatPS1(:,7);
PS1.gErr        = sqrt(CatPS1(:,8).^2 + SysErr.^2);
PS1.r           = CatPS1(:,14);
PS1.rErr        = sqrt(CatPS1(:,15).^2 + SysErr.^2);
PS1.i           = CatPS1(:,21);
PS1.iErr        = sqrt(CatPS1(:,22).^2 + SysErr.^2);

SDSS.RA         = CatSDSS(:,1)./RAD;
SDSS.Dec        = CatSDSS(:,2)./RAD;
SDSS.Epoch      = CatSDSS(:,3);
SDSS.g          = CatSDSS(:,4);
SDSS.gErr       = CatSDSS(:,5);
SDSS.r          = CatSDSS(:,6);
SDSS.rErr       = CatSDSS(:,7);
SDSS.i          = CatSDSS(:,8);
SDSS.iErr       = CatSDSS(:,9);

Out.gDiff = (PS1.g - SDSS.g')./sqrt(PS1.gErr.^2 + (SDSS.gErr.').^2);
Out.rDiff = (PS1.r - SDSS.r')./sqrt(PS1.rErr.^2 + (SDSS.rErr.').^2);
Out.iDiff = (PS1.i - SDSS.i')./sqrt(PS1.iErr.^2 + (SDSS.iErr.').^2);

% Note PA measured Westward
[D,PA] = celestial.coo.sphere_dist_fast(SDSS.RA,SDSS.Dec,PS1.RA.',PS1.Dec.');

Out.Epoch_SDSS = SDSS.Epoch;
Out.Epoch_PS1  = PS1.Epoch;
Out.Dist  = D.*RAD.*3600;    % [arcsec]
Out.TimeDiff = ((PS1.Epoch - SDSS.Epoch)./JYear);   % [Julian years]
Out.PM    = Out.Dist./Out.TimeDiff;  % [arcsec/yr]
Out.PA    = PA.*RAD;
Out.posMeanChi2 = PS1.posMeanChi2;



