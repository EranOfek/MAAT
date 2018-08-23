function Sim=add_fitscat_to_sim(Sim,CatFiles,varargin)
% Add a FITS catalogs into SIM images
% Package: ImUtil.pipe
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

CatField     = AstCat.CatField;
ColCellField = AstCat.ColCellField;

DefV.Section              = [];
DefV.ColX                 = {'XWIN_IMAGE','X_IMAGE','X','xpos'};
DefV.ColY                 = {'YWIN_IMAGE','Y_IMAGE','Y','ypos'};
DefV.NewColNameX          = '';
DefV.NewColNameY          = '';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Nsim = numel(Sim);
for Isim=1:1:Nsim
    [Out,~] = FITS.read_table(CatFiles{Isim});
    
    [~,ColX,~]=select_exist_colnames(Out,InPar.ColX(:));
    [~,ColY,~]=select_exist_colnames(Out,InPar.ColY(:));
    
    if (~isempty(InPar.Section))
        % select sources in image section
        FlagIn = Out.(CatField)(:,ColX)>=InPar.Section(Isim,1) & ...
                 Out.(CatField)(:,ColX)<=InPar.Section(Isim,2) & ...
                 Out.(CatField)(:,ColY)>=InPar.Section(Isim,3) & ...
                 Out.(CatField)(:,ColY)<=InPar.Section(Isim,4);
        % remove from catalog sources outside image section
        Out.(CatField) = Out.(CatField)(FlagIn,:);
        
        % correct X/Y coordinate in catalog to match the image section
        Out.(CatField)(:,ColX) = Out.(CatField)(:,ColX) - (InPar.Section(Isim,1)-1);
        Out.(CatField)(:,ColY) = Out.(CatField)(:,ColY) - (InPar.Section(Isim,3)-1);
    end
    % replace X/Y column names in catalog
    if (~isempty(InPar.NewColNameX))
        Out.(ColCellField){ColX} = InPar.NewColNameX;
    end
    if (~isempty(InPar.NewColNameY))
        Out.(ColCellField){ColY} = InPar.NewColNameY;
    end
    Sim(Isim).(CatField)     = Out.(CatField);
    Sim(Isim).(ColCellField) = Out.(ColCellField);
    
    Sim(Isim)                = colcell2col(Sim(Isim));
    
end
    
    
    