function Res=cat_fit_template2phot(Cat,varargin)
% SHORT DESCRIPTION HERE
% Package: AstroUtil
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=AstroUtil.spec.cat_fit_template2phot(S)
% Reliable: 
%--------------------------------------------------------------------------

CatField        = 'Cat';
ColField        = 'Col';

DefV.SearchRadius         = 2;
DefV.SearchRadiusUnits    = 'arcsec';
DefV.CooUnits             = 'deg';
DefV.ColMag               = {'MAG_PSF'};
DefV.ColMagErr            = {'MAGERR_PSF'};
DefV.ColRA                = {'RA','ALPHAWIN_J2000'};
DefV.ColDec               = {'Dec','DELTAWIN_J2000'};
DefV.CatName              = {'GAIADR2','PS1'};
DefV.CatColMag            = { {'Mag_G','Mag_BP','Mag_RP'}, {'gPSFMag','rPSFMag','iPSFMag','zPSFMag','yPSFMag'}};
DefV.CatColMagErr         = { {'MagErr_G','MagErr_BP','MagErr_RP'}, {'gPSFMagErr','rPSFMagErr','iPSFMagErr','zPSFMagErr','yPSFMagErr'}};
DefV.CatFilterFamily      = { {'GAIA','GAIA','GAIA'}, {'PAN-STARRS/PS1','PAN-STARRS/PS1','PAN-STARRS/PS1','PAN-STARRS/PS1','PAN-STARRS/PS1'}};
DefV.CatFilterName        = { {'G','BP','RP'}, {'g','r','i','z','y'}};
DefV.CatFilterType        = { {'Vega','Vega','Vega'}, {'AB','AB','AB','AB','AB'}};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Nsrc = size(Cat.(CatField),1);

[~,ColRA,~]     = select_exist_colnames(Cat,InPar.ColRA(:));
[~,ColDec,~]    = select_exist_colnames(Cat,InPar.ColDec(:));
[~,ColMag,~]    = select_exist_colnames(Cat,InPar.ColMag(:));
[~,ColMagErr,~] = select_exist_colnames(Cat,InPar.ColMagErr(:));

CooUnitsFactor = convert.angular(InPar.CooUnits,'rad');

%RA  = Cat.(CatField)(:,ColRA)  .* CooUnitsFactor;
%Dec = Cat.(CatField)(:,ColDec) .* CooUnitsFactor;
%[CX,CY,CZ] = celestial.coo.coo2cosined(RA,Dec);
%[MeanRA,MeanDec] = celestial.cosined2coo(median(CX),median(CY),median(CZ));
%D = celestial.sphere_dist_fast(MeanRA,MeanDec,RA,Dec);
%Radius = max(D).*(1+100.*eps);

Ncat = numel(InPar.CatName);
Nf   = nan(Ncat,1);
for Icat=1:1:Ncat
    ExtCat(Icat) = catsHTM.sources_match(InPar.CatName{Icat},Cat,'ColRA',InPar.ColRA,'ColDec',InPar.ColDec,'CooUnits',InPar.CooUnits,...
                                            'SearchRadius',InPar.SearchRadius,'SearchRadiusUnits',InPar.SearchRadiusUnits);
                                        
    Nf(Icat) = numel(InPar.CatColMag{Icat});  % number of filters in each external catalog
end



SynMag     = [];
ExtrapFlag = [];
for Isrc=1:1:Nsrc
    PhotData = cell(sum(Nf),5);
    K = 0;
    for Icat=1:1:Ncat
        for If=1:1:Nf(Icat)
            Mag     = ExtCat(Icat).(CatField)(Isrc,ExtCat(Icat).(ColField).(InPar.CatColMag{Icat}{If}));
            MagErr  = ExtCat(Icat).(CatField)(Isrc,ExtCat(Icat).(ColField).(InPar.CatColMagErr{Icat}{If}));
            Family  = InPar.CatFilterFamily{Icat}{If};
            Band    = InPar.CatFilterName{Icat}{If};
            MagType = InPar.CatFilterType{Icat}{If};
            K = K + 1;
            PhotData(K,:) = {Mag, MagErr, Family, Band, MagType};
        end
    end
    
    % in the first time will calculate syn mag
    [Res(Isrc),SynMag,ExtrapFlag]=AstroUtil.spec.fit_template2phot(PhotData,'SynMag',SynMag,'ExtrapFlag',ExtrapFlag);
    Isrc
    
    
end

            
     