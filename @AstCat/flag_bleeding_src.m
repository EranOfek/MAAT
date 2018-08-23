function FlagCat=flag_bleeding_src(Sim,varargin)
% Flag sources in catalog which are concenrrated along the same lines/rows
% Package: @AstCat
% Description: Flag sources in catalog with anomalous concenration along
%              the same lines/rows.
% Input  : - A Sim class object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColX'  - A cell array of column names that may contains the X
%                      axis column name. Default is {'XWIN_IMAGE','X'}.
%            'ColY'  - A cell array of column names that may contains the Y
%                      axis column name. Default is {'YWIN_IMAGE','Y'}.
%            'Threshold' - Mark sources which number in the same
%                      column/rows is Threshold times above the robust std
%                      of number of sources over all columns/rows.
%                      Default is 20.
%            'BinWidth' - Column/rows bin width. Default is 3.
% Output : - A structure array with the following fields:
%            'Flag' - a vector of logical flag indicating sources that
%                     belong to a large group of sources in the same
%                     lines/rows.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: FlagCat=flag_bleeding_src(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

CatField  = AstCat.CatField;

DefV.ColX                 = {'XWIN_IMAGE','X','xpos'};
DefV.ColY                 = {'YWIN_IMAGE','Y','ypos'};
DefV.Threshold            = 20;
DefV.BinWidth             = 3;
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
    
    [Nx,Ex] = histcounts(Sim(Isim).(CatField)(:,ColX),'BinLimits',[MinX MaxX],'BinWidth',InPar.BinWidth);
    StdNx   = max(Util.stat.rstd(Nx(:)),1);
    [Ny,Ey] = histcounts(Sim(Isim).(CatField)(:,ColY),'BinLimits',[MinX MaxX],'BinWidth',InPar.BinWidth);
    StdNy   = max(Util.stat.rstd(Ny(:)),1);
    
    IBx = find(Nx>(InPar.Threshold.*StdNx));
    IBy = find(Ny>(InPar.Threshold.*StdNy));
    
    Nsrc = size(Sim(Isim).(CatField),1);
    FlagCat(Isim).Flag = false(Nsrc,1);
    for Ix=1:1:numel(IBx)
        X1 = Ex(IBx(Ix));
        X2 = X1+InPar.BinWidth;
        
        FlagCat(Isim).Flag = FlagCat(Isim).Flag | (Sim(Isim).(CatField)(:,ColX)>X1 & Sim(Isim).(CatField)(:,ColX)<X2 );
    end
    for Iy=1:1:numel(IBy)
        Y1 = Ey(IBy(Iy));
        Y2 = Y1+InPar.BinWidth;
        
        FlagCat(Isim).Flag = FlagCat(Isim).Flag | (Sim(Isim).(CatField)(:,ColY) > Y1 & Sim(Isim).(CatField)(:,ColY)<Y2 );
    end
    
end
