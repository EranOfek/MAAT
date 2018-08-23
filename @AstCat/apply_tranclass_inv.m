function Cat=apply_tranclass_inv(Cat,TranC,varargin)
% Apply inverse transformation to X,Y coordinates in an AstCat object
% Package: @AstCat
% Description: Apply inverse transformation to X,Y coordinates in an
%              AstCat object. See also TranClass/apply_tranclass_inv.
% Input  : - An AstCat object.
%          - A TranClass object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Thresh'- Accuracy for convergence. Default is 1e-3.
%            'ColX' - A cell array of column names containing X coordinates
%                     to be transformed. Transformation is applied only for
%                     existing columns.
%                     Default is
%                     {'XWIN_IMAGE','X','XPEAK_IMAGE','X_IMAGE'}.
%            'ColY' - A cell array of column names containing Y coordinates
%                     to be transformed. Transformation is applied only for
%                     existing columns.
%                     Default is
%                     {'YWIN_IMAGE','Y','YPEAK_IMAGE','Y_IMAGE'}.
% Output : - An AstCat object with the transformed coordinates.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=apply_tranclass_inv(Cat,TranC);
% Reliable: 2
%--------------------------------------------------------------------------

CatField = AstCat.CatField;
            
DefV.Thresh               = 1e-3;
DefV.ColX                 = {'XWIN_IMAGE','X','XPEAK_IMAGE','X_IMAGE'};
DefV.ColY                 = {'YWIN_IMAGE','Y','YPEAK_IMAGE','Y_IMAGE'};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.ColX))
    InPar.ColX = {InPar.ColX};
end
if (~iscell(InPar.ColY))
    InPar.ColY = {InPar.ColY};
end

if (numel(Cat)>1)
    error('Cat should contain a single element');
end

% for each X/Y coordinate column in the AstCat object
NcolX = numel(InPar.ColX);
NcolY = numel(InPar.ColY);
Ncol  = max(NcolX,NcolY);
for Icol=1:1:Ncol
    IcolX = min(NcolX,Icol);
    IcolY = min(NcolY,Icol);

    ColX  = colname2ind(Cat,InPar.ColX{IcolX});
    ColY  = colname2ind(Cat,InPar.ColY{IcolY});
    if (~isnan(ColX) && ~isnan(ColY))
        X = col_get(Cat,ColX);
        Y = col_get(Cat,ColY);
        
        [OutX,OutY] = apply_tranclass_inv(TranC,X,Y,InPar.Thresh);
        Cat.(CatField)(:,ColX) = OutX;
        Cat.(CatField)(:,ColY) = OutY;
    end
end



