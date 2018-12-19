function []=acis_search_flares(varargin)
% SHORT DESCRIPTION HERE
% Package: ImUtil
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

CatField = AstCat.CatField;
ColField = AstCat.ColField;

DefV.EnergyRange          = [200 10000];
DefV.SelectGoodTimes      = true;
DefV.NinBin               = 100;
DefV.FlareExp             = 1000;

DefV.BrightPixSize        = 4;
DefV.PixSize              = 1;

DefV.ChipSize             = [1024 1024];  % x,y
DefV.timefield            = 'time';
DefV.ccdfield             = 'ccd_id';
DefV.chipxfield           = 'chipx';
DefV.chipyfield           = 'chipy';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


[EvtT,MetaData,Image] = ImUtil.X.read_chandra_acis('acis.evt2',...
                            'EnergyRange',InPar.EnergyRange,...  
                            'SelectGoodTimes',InPar.SelectGoodTimes,...
                            'NinBin',InPar.NinBin');

EvtT = astcat_array2table(EvtT);                        
UniqueCCD = unique(EvtT.(CatField).(InPar.ccdfield));
Nccd      = numel(UniqueCCD);
for Iccd=1:1:Nccd
    FlagCCD = EvtT.(CatField).(InPar.ccdfield) == UniqueCCD(Iccd);
    
    % generate image in order to remove bright sources and estimate
    % background level
    VecX = (1:InPar.BrightPixSize:InPar.ChipSize(1))';
    VecY = (1:InPar.BrightPixSize:InPar.ChipSize(2))';

    ImageB = histcounts2(EvtT.(CatField).(InPar.chipxfield)(FlagCCD),...
                         EvtT.(CatField).(InPar.chipyfield)(FlagCCD),...
                         VecX,VecY);
    
    CntRate  = sum(ImageB(:))./MetaData.ExpTime./prod(InPar.ChipSize);
    Fback    = ImageB<max(2,median(ImageB(:))+std(ImageB(:)).*5);
    BackRate = sum(ImageB(Fback))./MetaData.ExpTime./prod(InPar.ChipSize);

    
    % generate an image for CCD
    VecX = (1:InPar.PixSize:InPar.ChipSize(1))';
    VecY = (1:InPar.PixSize:InPar.ChipSize(2))';

    ImageB = histcounts2(EvtT.(CatField).(InPar.chipxfield)(FlagCCD),...
                         EvtT.(CatField).(InPar.chipyfield)(FlagCCD),...
                         VecX,VecY);
    
    % filter image
    
    
    