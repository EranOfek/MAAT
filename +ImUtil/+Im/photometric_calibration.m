function [Res]=photometric_calibration(Cat,varargin)
% SHORT DESCRIPTION HERE
% Package: ImUtil
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res]=ImUtil.Im.photometric_calibration(Cat);
% Reliable: 
%--------------------------------------------------------------------------

CatField     = AstCat.CatField;
ColCellField = AstCat.ColCellField;

DefV.CalibMethod          = 'simple';
DefV.CatName              = 'PS1'; %'SDSSDR10';
DefV.Band                 = 'g';
DefV.Color                = {'g','r'};
DefV.ColCell              = {};
DefV.ColRA                = {'RA','ALPHAWIN_J2000'};
DefV.ColDec               = {'Dec','DELTAWIN_J2000'};
DefV.ColMag               = {'Mag','MAG_PSF'};
DefV.ColMagErr            = {'MagErr','MAGERR_PSF'};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


switch lower(InPar.CatName)
    case 'gaiadr2'
        ColMagBand = 'Mag_G';
        ColMagCol1 = 'Mag_BP';
        ColMagCol2 = 'Mag_RP';
        
    case 'sdssdr10'
        MagBaseName = 'modelMag_%s';
        ColMagBand  = sprintf(MagBaseName,InPar.Band);
        ColMagCol1  = sprintf(MagBaseName,InPar.Color{1});
        ColMagCol2  = sprintf(MagBaseName,InPar.Color{2});
        
    case 'ps1'
        MagBaseName = '%sPSFMag';
        ColMagBand  = sprintf(MagBaseName,InPar.Band);
        ColMagCol1  = sprintf(MagBaseName,InPar.Color{1});
        ColMagCol2  = sprintf(MagBaseName,InPar.Color{2});
        
    otherwise
        error('Unknwon CatName option');
end



if (~AstCat.isastcat(Cat))
    Tmp = Cat;
    Cat = AstCat;
    Cat.(CatField) = Tmp;
    Cat.(ColCellField) = InPar.ColCell;
end
% else: Cat is already an AstCat object

% select relevant columns
% RA, Dec, and instrumental magnitude
[~,Col.RA,~]     = select_exist_colnames(Cat,InPar.ColRA(:));
[~,Col.Dec,~]    = select_exist_colnames(Cat,InPar.ColDec(:));
[~,Col.Mag,~]    = select_exist_colnames(Cat,InPar.ColMag(:));
[~,Col.MagErr,~] = select_exist_colnames(Cat,InPar.ColMagErr(:));

% match catalog against GAIA and another catalog

% for each star in GAIA
% 1. fit PS1 magnitudes to stellar spectral templates
% 2. Ft all PS1 mag with spectral template
% 3. select sources which have a good fit (small rms)
% 4. Do synPhot on selected targets: GAIA band + band(parameters)
% 5. calc
% mag_inst = ZP + Mag_GAIA-SynMag(BandGAIA/exact) +SynMagTemp(Band, BestFitTemp+airmass,Telluric,slope,...)
% and iterate parameters until convergence


CatM = catsHTM.sources_match(InPar.CatName,Cat);
Nsrc = size(Cat.(CatField),1);

Mag_inst    = Cat.(CatField)(:,Col.Mag);
MagErr_inst = Cat.(CatField)(:,Col.MagErr);


Mag_Band  = col_get(CatM,ColMagBand);
Mag_Color1 = col_get(CatM,ColMagCol1);
Mag_Color2 = col_get(CatM,ColMagCol2);

Flag   = ~isnan(Mag_Band) & ~isnan(Mag_Color1) & ~isnan(Mag_Color2) & ~isnan(Mag_inst) & ~isnan(MagErr_inst) & Mag_Band>14.5 & Mag_Band<17;
H = [ones(Nsrc,1), Mag_Color1-Mag_Color2];
Par = H(Flag,:)\[Mag_inst(Flag) - Mag_Band(Flag)];
Resid = [Mag_inst(Flag) - Mag_Band(Flag)] - H(Flag,:)*Par;
RStd  = Util.stat.rstd(Resid);

Res.CalibMethod          = InPar.CalibMethod;
Res.CatName              = InPar.CatName;
Res.Par                  = Par;
Res.Resid                = Resid;
Res.RStd                 = RStd;








