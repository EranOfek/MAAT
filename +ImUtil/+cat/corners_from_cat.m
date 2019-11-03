function [Res]=corners_from_cat(Cat,varargin)
% Estimate catalog corners RA/Dec position from catalog.
% Package: ImUtil
% Description: 
% Input  : - A matrix with at least 4 colums, of RA, Dec, X, Y coordinates.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'CCDSEC'- [xmin, xmax, ymin, ymax] of corners.
%                      Default is [3172 3172].
%            'ColRA' - Index of column containing the Longitude/RA.
%                      Default is 1.
%            'ColDec'- Index of column containing the Latitude/Dec.
%                      Default is 2.
%            'ColX'  - Index of column containing the X.
%                      Default is 3.
%            'ColY'  - Index of column containing the Y.
%                      Default is 4.
%            'Units' - RA/Dec coordinates units {'rad'|'deg'}.
%                      Default is 'rad'.
%            'FitFun'- Design matrix funcruin generator.
%                      Default is @(x,y) [x, y, x.*y]
% Output : - A structure with the corners/center information.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Res]=ImUtil.cat.corners_from_cat(Cat);
% Reliable: 
%--------------------------------------------------------------------------

DefV.CCDSEC               = [1 3080 1 3072];
DefV.Units                = 'rad';
DefV.ColRA                = 1;
DefV.ColDec               = 2;
DefV.ColX                 = 3;
DefV.ColY                 = 4;
DefV.FitFun               = @(x,y) [ ones(size(x)), x, y, x.*y, x.^2, y.^2];
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Xcorner = InPar.CCDSEC([1 1 2 2]).'; 
Ycorner = InPar.CCDSEC([3 4 3 4]).';


% convert coordinates to radiuans
Factor = convert.angular(InPar.Units,'rad');
Cat(:,[InPar.ColRA, InPar.ColDec]) = Cat(:,[InPar.ColRA, InPar.ColDec]).*Factor;

% convert RA/Dec to cosine directions
[CD1,CD2,CD3] = celestial.coo.coo2cosined(Cat(:,InPar.ColRA),Cat(:,InPar.ColDec));

% fit CDn = a*X +b*Y + c*X*Y + d*X^2 + e*Y^2
H = InPar.FitFun(Cat(:,InPar.ColX), Cat(:,InPar.ColY));
Hc= InPar.FitFun(Xcorner, Ycorner);
P1 = H\CD1;
P2 = H\CD2;
P3 = H\CD3;

Resid1 = CD1 - H*P1;
Resid2 = CD2 - H*P2;
Resid3 = CD3 - H*P3;
Std1   = std(Resid1);
Std2   = std(Resid2);
Std3   = std(Resid3);

CornerCD1 = Hc*P1;
CornerCD2 = Hc*P2;
CornerCD3 = Hc*P3;

[CornerRA, CornerDec] = celestial.coo.cosined2coo(CornerCD1, CornerCD2, CornerCD3);

CenterCD1 = mean(CornerCD1);
CenterCD2 = mean(CornerCD2);
CenterCD3 = mean(CornerCD3);
[CenterRA, CenterDec] = celestial.coo.cosined2coo(CenterCD1, CenterCD2, CenterCD3);

Radius = celestial.coo.sphere_dist_fast(CenterRA, CenterDec, CornerRA, CornerDec);

Res.CornerRA  = CornerRA;
Res.CornerDec = CornerDec;
Res.CornerCD  = [CornerCD1, CornerCD2, CornerCD3];
Res.CenterRA  = CenterRA;
Res.CenterDec = CenterDec;
Res.CenterCD  = [CenterCD1, CenterCD2, CenterCD3];
Res.Radius    = max(Radius);
Res.FitStd    = [Std1, Std2, Std3];
