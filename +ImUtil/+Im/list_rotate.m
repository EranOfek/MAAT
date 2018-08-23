function Cat=list_rotate(Cat,RotAng,RefCoo,varargin)
% Rotate coordinates in a list
% Package: @ImUtil.Im
% Description: Rotate X/Y coordinates in each catalog in a list.
%              See also AstCat/cat_rotate
% Input  : - A matrix in which two of the columns are X, Y coordinates.
%          - Scalar of rotation angle [deg].
%          - Two element vector of [X,Y] centers around to rotate the
%            coordinates.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColX' - Cell array of X axis column name. Use the first
%                     existing column. Default is {'XWIN_IMAGE','X'}.
%            'ColY' - Cell array of Y axis column name. Use the first
%                     existing column. Default is {'YWIN_IMAGE','Y'}.
% Output : - AStCat object with rotated coordinates.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [CatRot]=ImUtil.Im.list_rotate(AstC,RotAng,RefCoo);
% Reliable: 
%--------------------------------------------------------------------------


if (nargin<3)
    RefCoo = [0 0];
    if (nargin<2)
        RotAng = 0;
    end
end

DefV.ColX                 = 1;
DefV.ColY                 = 2;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


X = Cat(:,InPar.ColX);
Y = Cat(:,InPar.ColY);

RefX  = RefCoo(1);
RefY  = RefCoo(2);

Cat(:,InPar.ColX) = (X - RefX).*cosd(RotAng) - (Y - RefY).*sind(RotAng);
Cat(:,InPar.ColY) = (X - RefX).*sind(RotAng) + (Y - RefY).*cosd(RotAng);
