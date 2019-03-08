function [Res,OutCat,Cat]=test_match(Cat,Ref,TranC,varargin)
% Apply transformation to catalog, match two lists and calc residuls.
% Package: ImUtil.pattern
% Description: Given two lists of [X,Y] coordinates (Cat and Ref). 
%              Apply the transformation to the Cat list, match it with
%              Ref and calculate the residuls.
% Input  : - Cat matrix with [X,Y] coordinates.
%          - Ref matrix with [X,Y] coordinates.
%          - A TranClass object or a vectotr of
%            [ShiftX, ShiftY]  (shift)
%            or [ShiftX, ShiftY, Rot(rad)]  (rotation)
%            or [ShiftX, ShiftY, RotX(rad), RotY(rad)]  (affine)
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'RotationUnits' - Rotation units (if TranC is a vector).
%                      Default is 'rad'.
%            'Flip' - [FlipX, FlipY] sign. Default is [1 1].
%                     Flip done before transiformation.
% Output : - Structure of residuls  information:
%            'IndRef' - indices of Ref catalog matched with OutCat.
%            'IndCat' - indices of OutCat catalog matched with Ref.
%            'ResidX' - X axis residuls of matched sources.
%            'ResidY' - Y axis residuls of matched sources.
%            'StdX'   - Std of X axis residuls.
%            'StdY'   - Std of Y axis residuls.
%          - Transformed Cat matrix
%          - Original input Cat but sorted like OutCat.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=sortrows(rand(1000,2),2); Ref=rand(1000,2);
%          Cat = Cat.*1000; Ref = Ref.*1000;
%          Cat(1:50,:) = Ref(1:50,:) + [10 -17];
%          Cat = sortrows(Cat,2);
%          [Res,OutCat]=ImUtil.pattern.test_match(Cat,Ref,-[10 -17]);
% Reliable: 
%--------------------------------------------------------------------------

DefV.Radius           = 2;
DefV.ColXc            = 1;            
DefV.ColYc            = 2;
DefV.ColXr            = 1;
DefV.ColYr            = 2;
DefV.RotationUnits    = 'rad';
DefV.Flip             = [1 1];
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% apply transformation
% first apply rotation and only than apply the shift

% Rotation = convert.angular(InPar.RotationUnits,'rad',Rotation);  % [rad]
% 
% XcatRot = Xcat.*cosd(Rotation(Irot)) - Ycat.*sind(Rotation(Irot));
% YcatRot = Xcat.*sind(Rotation(Irot)) + Ycat.*cosd(Rotation(Irot));


[OutCatX,OutCatY] = TranClass.apply_tran(Cat(:,InPar.ColXc).*InPar.Flip(1),...
                                         Cat(:,InPar.ColYc).*InPar.Flip(2),...
                                         [TranC(1:2),TranC(3)], InPar.RotationUnits);
OutCat = [OutCatX, OutCatY];

if (~issorted(OutCat(:,2)))
    [OutCat,SI] = sortrows(OutCat,2);
    Cat         = Cat(SI,:);
end

% match catalogs
%[~,Matched] = VO.search.match_cats_pl(OutCat(:,[InPar.ColXc, InPar.ColYc]),Ref(:,[InPar.ColXr, InPar.ColYr]),'Radius',InPar.Radius);
[~,Matched] = VO.search.match_cats_pl(OutCat(:,[1, 2]),Ref(:,[InPar.ColXr, InPar.ColYr]),'Radius',InPar.Radius);


% [Res11,~,~]=search_cat(OutCat(:,[InPar.ColXc, InPar.ColYc]),Ref(:,[InPar.ColXr, InPar.ColYr]),[],...
%                                     'CooType','plane',...
%                                     'SearchRad',InPar.Radius);
%                                     

% remove sources with more than one match
FlagN1 = [Matched.Num]==1;
IndCat = [Matched(FlagN1).IndCat];
IndRef = [Matched(FlagN1).IndRef];

% output
Res.IndCat = IndCat;
Res.IndRef = IndRef;
Res.ResidX = OutCat(IndCat,1) - Ref(IndRef,InPar.ColXr);
Res.ResidY = OutCat(IndCat,2) - Ref(IndRef,InPar.ColYr);

Res.StdX   = nanstd(Res.ResidX);
Res.StdY   = nanstd(Res.ResidY);


