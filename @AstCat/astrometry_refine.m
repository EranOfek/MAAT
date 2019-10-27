function []=astrometry_refine(AC,varargin)
% SHORT DESCRIPTION HERE
% Package: @AstCat
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SearchRadius' - Search radius. Default is 2.
%            'SearchRadiusUnits' - Search radius units.
%                       Default is 'arcsec'.
%            'ColCell' - Default is {}.
%            'ColRA' - Default is {'RA','ALPHAWIN_J2000'}.
%            'ColDec' - Default is {'Dec','DELTAWIN_J2000'}.
%            'CooUnits' - Input catalog coordinates units.
%                       Default is 'rad'.
%            'ColDecHTM' - Default is 2.
%            'ColRAHTM'  - Default is 1.
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.EpochJD              = 

DefV.ColX                 = 
DefV.ColY                 =
DefV.ColRA                = {'RA','ALPHAWIN_J2000'};
DefV.ColDec               = {'Dec','DELTAWIN_J2000'};
DefV.CooUnits             = 'rad';
DefV.AstromCat            = 'GAIADR2';
DefV.SearchRadius         = 2;
DefV.SearchRadiusUnits    = 'arcsec';

DefV.ApplyPM            = true; %true; %true;
DefV.ApplyParallax      = false;
DefV.RC_ColPM_RA        = 'PMRA';
DefV.RC_ColPM_Dec       = 'PMDec';
DefV.RC_ColPlx          = 'Plx';
DefV.RC_ColRV           = 'RV';
DefV.RC_EpochInRA       = 'Epoch';
DefV.RC_EpochInDec      = 'Epoch';
DefV.RC_EpochInUnits    = 'yr';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);



% for each AstCat element
Nac = numel(AC);
for Iac=1:1:Nac
    %
    
    % match catalog against GAIA
    CatM = catsHTM.sources_match(InPar.AstromCat,AC(Iac),...
                                'OutType','astcat',...
                                'SearchRadius',InPar.SearchRadius,...
                                'SearchRadiusUnits',InPar.SearchRadiusUnits,...
                                'ColRA',InPar.ColRA,...
                                'ColDec',InPar.ColDec,...
                                'CooUnits',InPar.CooUnits);
                            
    % apply proper motion
    if InPar.ApplyPM
        RefCat = apply_proper_motion(RefCat,'EpochInRA',InPar.RC_EpochInRA,...
                                            'EpochInDec',InPar.RC_EpochInDec,...
                                            'EpochInUnits',InPar.RC_EpochInUnits,...
                                            'EpochOut',InPar.EpochJD,...
                                            'EpochOutUnits','JD',...
                                            'ColPM_RA',InPar.RC_ColPM_RA,...
                                            'ColPM_Dec',InPar.RC_ColPM_Dec,...
                                            'ColPlx',InPar.RC_ColPlx,...
                                            'ColRV',InPar.RC_ColRV, ...
                                            'ApplyParallax', InPar.ApplyParallax);
    end
    
    % selected columns
    
    % projection
    [X,Y]=projection(RC,'tan',[RC_ColRA RC_ColDec],[RAD.*3600./Scale RA Dec],'rad');
    
    % fit
    
           % Special treatment for PV distortions
           % convert coordinates from pixels to deg
           CD = [1 0; 0 1].*Scale./3600;
           %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!
           %clear many appearence of the same object.
           [indexes,ia,ic]=unique(MatchedCat(:,Sim(Isim).Col.IndexSimYsorted));
           MatchedRef=MatchedRef(ia,:);
           MatchedCat=MatchedCat(ia,:);
           %!!!!!!!!!!!!!!!!!!!-----------------------!!!!!!!!!!!!!!!!!!!!!

           MatchedRefCD        = MatchedRef;
           MatchedRefCD(:,1:2) = [CD*MatchedRefCD(:,1:2)']';
           MatchedCatCD        = MatchedCat;
           MatchedCatCD(:,[ColXc, ColYc]) = [CD*MatchedCatCD(:,[ColXc, ColYc])']';
           
           ResAst(Isim) = ImUtil.pattern.fit_transform(MatchedRefCD,MatchedCatCD,TranC,'ImSize',ImSize(Isim,:),...
                                                                          'BlockSize',InPar.AnalysisBlockSize,...
                                                                          'PixScale',InPar.Scale,...
                                                                          'CooUnits','deg',...
                                                                          'NormXY',1,...
                                                                          'ColCatX',ColXc,...
                                                                          'ColcatY',ColYc,...
                                                                          'PolyMagDeg',3,...
                                                                          'StepMag',0.1,...
                                                                          'Niter',InPar.Niter,...
                                                                          'SigClip',InPar.SigClip,...
                                                                          'MaxResid',InPar.MaxResid,...
                                                                          'Plot',InPar.Plot);
           
    
end