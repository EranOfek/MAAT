function FlagCat=flag_overdense_src(Sim,varargin)
% Flag sources in catalog in regions with very high overdensity
% Package: @AstCat
% Description: Flag sources in catalog in regions with very high overdensity
% Input  : - A Sim class object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColX'  - A cell array of column names that may contains the X
%                      axis column name. Default is {'XWIN_IMAGE','X'}.
%            'ColY'  - A cell array of column names that may contains the Y
%                      axis column name. Default is {'YWIN_IMAGE','Y'}.
%            'Threshold' - Mark sources which numberof sources is that
%                      times above the rstd.
%                      Default is 5.
%            'BinWidth' - Column/rows bin width. Default is 32.
% Output : - A structure array with the following fields:
%            'Flag' - a vector of logical flag indicating sources that
%                     are found in overdense region.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: FlagCat=flag_overdense_src(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

CatField  = AstCat.CatField;

DefV.ColX                 = {'XWIN_IMAGE','X','xpos'};
DefV.ColY                 = {'YWIN_IMAGE','Y','ypos'};
DefV.Threshold            = 5;
DefV.MinStd               = 0.9;
DefV.BinWidth             = 32;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);



Nsim = numel(Sim);
for Isim=1:1:Nsim
    [~,ColInd,~]=select_exist_colnames(Sim(Isim),InPar.ColX(:));
    ColX = ColInd(find(~isnan(ColInd),1));
    [~,ColInd,~]=select_exist_colnames(Sim(Isim),InPar.ColY(:));
    ColY = ColInd(find(~isnan(ColInd),1));
    
    MinX = floor(min(Sim(Isim).(CatField)(:,ColX)));
    MaxX = ceil(max(Sim(Isim).(CatField)(:,ColX)));
    MinY = floor(min(Sim(Isim).(CatField)(:,ColY)));
    MaxY = ceil(max(Sim(Isim).(CatField)(:,ColY)));
    
    [Nxy,XE,YE] = histcounts2(Sim(Isim).(CatField)(:,ColX),Sim(Isim).(CatField)(:,ColY),...
                       'XBinLimits',[MinX MaxX],...
                       'YBinLimits',[MinY MaxY],...
                       'BinWidth',InPar.BinWidth);
    
    Std = max(Util.stat.rstd(Nxy(:)),InPar.MinStd);
    
    Inxy  = find(Nxy>(Std.*InPar.Threshold));
    [J,I] = ind2sub(size(Nxy),Inxy);
    Ni    = numel(Inxy); 
    Nsrc  = size(Sim(Isim).Cat,1);
    FlagCat(Isim).Flag = false(Nsrc,1);
    for Ii=1:1:Ni
        Ylim1 = YE(I(Ii));
        Ylim2 = YE(I(Ii)+1);
        Xlim1 = XE(J(Ii));
        Xlim2 = XE(J(Ii)+1);
        
        Flag = Sim(Isim).(CatField)(:,ColX)>Xlim1 & Sim(Isim).(CatField)(:,ColX)<Xlim2 & ...
               Sim(Isim).(CatField)(:,ColY)>Ylim1 & Sim(Isim).(CatField)(:,ColY)<Ylim2;
           
        FlagCat(Isim).Flag = FlagCat(Isim).Flag | Flag;
    end
                                                   
end
