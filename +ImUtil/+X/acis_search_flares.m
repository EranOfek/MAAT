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

DefV.EvtFileName          = 'acis.evt2';
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
DefV.xfield               = 'x';
DefV.yfield               = 'y';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


[EvtT,MetaData,Image] = ImUtil.X.read_chandra_acis(InPar.EvtFileName,...
                            'EnergyRange',InPar.EnergyRange,...  
                            'SelectGoodTimes',InPar.SelectGoodTimes,...
                            'NinBin',InPar.NinBin');

EvtT = astcat_array2table(EvtT);                        

%UniqueCCD = unique(EvtT.(CatField).(InPar.ccdfield));
%Nccd      = numel(UniqueCCD);


% This should work in sky coordinates (x,y)
MinX = min(EvtT.Cat.(InPar.xfield));
MinY = min(EvtT.Cat.(InPar.yfield));
MaxX = max(EvtT.Cat.(InPar.xfield));
MaxY = max(EvtT.Cat.(InPar.yfield));

% construct a full image
% on scale of ~100x100 pix:
% 1. get local exp map
% 2. estimate back
% 3. run MF
% 4. Repeat on 1000 s image
% save LC for all cand in ~90% circle.





for Iccd=1:1:Nccd
    FlagCCD = EvtT.(CatField).(InPar.ccdfield) == UniqueCCD(Iccd);
    
    % generate image in order to remove bright sources and estimate
    % background level
    % Note that sky coordinates (x,y) should be used
    
    VecX = (1:InPar.BrightPixSize:InPar.ChipSize(1))';
    VecY = (1:InPar.BrightPixSize:InPar.ChipSize(2))';

    ImageB = histcounts2(EvtT.(CatField).(InPar.chipxfield)(FlagCCD),...
                         EvtT.(CatField).(InPar.chipyfield)(FlagCCD),...
                         VecX,VecY);
    
                     
    % do this for FlareExp time blocks...
                     
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
    
    
    