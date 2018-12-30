function []=fit_phot_transmission(varargin)
% SHORT DESCRIPTION HERE
% Package: ImUtil.calib
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Dec 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.ColorCatAstF         = 'Pan-STARRS/PS1';
DefV.ColorCatBands        = {'g','r','i','z','y'};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

%%
% get spectral template
S=AstSpec.get_pickles;


% try also GAIA spec...
%S=AstSpec.get_all_gaia_synspec;

% do we need to apply extinction and Telluric absorption???

Opt = cats.spec.AtmoExtinction.SNfactory;
IR = cats.spec.AtmoExtinction.MaunaKea_IR_AM1;
IR(:,2)= -2.5.*log10(IR(:,2));


S = atmospheric_extinction(S,-1.2,'KPNO');


% synthetic photometry for spectral template in Color cat bands
Nfilt = numel(InPar.ColorCatBands);
for Ifilt=1:1:Nfilt
    Pick.(InPar.ColorCatBands{Ifilt}) = synphot(S,InPar.ColorCatAstF,InPar.ColorCatBands{Ifilt},'AB');
end

% get catalog
[CatColor] = catsHTM.cone_search('PS1',1,1,1000,'outtype','AstCat');
[CatMag] = catsHTM.cone_search('GAIADR2',1,1,1000,'outtype','AstCat');
MatchedCatColor = match(CatColor,CatMag);

% select stars with good photometry
FlagG = MatchedCatColor.Cat(:,8)<0.02 & ...
        MatchedCatColor.Cat(:,15)<0.02 & ...
        MatchedCatColor.Cat(:,22)<0.02 & ...
        MatchedCatColor.Cat(:,29)<0.02 & ...
        MatchedCatColor.Cat(:,36)<0.02 & ...
        CatMag.Cat(:,13)<0.5;
    
MatchedCatColor.Cat   = MatchedCatColor.Cat(FlagG,:);

% find best fit spectral template for each star
Nstar = size(MatchedCatColor.Cat,1);

%%
K = 0;
clear CatCalib
for Istar=1:1:Nstar
    % diff between catalog magnitude and spectral template mag
    Diff = [MatchedCatColor.Cat(Istar,7)-Pick.g, MatchedCatColor.Cat(Istar,14)-Pick.r, MatchedCatColor.Cat(Istar,21)-Pick.i, MatchedCatColor.Cat(Istar,28)-Pick.z]; % MatchedCatColor.Cat(Istar,35)-Pick.y];
    % rms per spectral template
    Std  = std(Diff,[],2);
    % best spectral template
    [Min(Istar), Imin(Istar)] = min(Std);
    
    if (Min(Istar)<0.04)
        % good star for calibration
        K = K + 1;
        CatCalib(K).CatColor = MatchedCatColor.Cat(Istar,:);
        CatCalib(K).BestTemplateInd = Imin(Istar);
        
        ScaledSpec = scale2mag(S(Imin(Istar)),CatMag.Cat(Istar,16),'GAIA','G','Vega');
        BP = synphot(ScaledSpec,'GAIA','BP','Vega');
        RP = synphot(ScaledSpec,'GAIA','RP','Vega');
        G  = synphot(ScaledSpec,'GAIA','G','Vega');
        
        %check that synthetic BP and RP are consistent with GAIA BP and RP
        CatCalib(K).Gdiff  = G  - CatMag.Cat(Istar,16);
        CatCalib(K).BPdiff = BP - CatMag.Cat(Istar,18);
        CatCalib(K).RPdiff = RP - CatMag.Cat(Istar,20);
        
        
    end
    
end



