function Cat=update_coordinates(Cat,varargin)
% Update RA/Dec coordinates in an SIM/AstCat object based on a new WCS
% Package: @SIM
% Description: Given a SIM object in which the catalog and header are
%              populated, use the header information to update the RA/Dec
%              coordinates based on the X/Y coordinates and header WCS.
% Input  : - A SIM object containing a catalog and header.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColNameRA'  - RA column to update.
%                           Default is 'ALPHAWIN_J2000'.
%            'ColNameDec' - Dec column to update.
%                           Default is 'DELTAWIN_J2000'.
%            'ColNameX'   - Cell of possible column names in which to look
%                           for the X coordinates. Use the first existing
%                           column. Default is
%                           {'XWIN_IMAGE','X','X_IMAGE'}.
%            'ColNameY'   - Cell of possible column names in which to look
%                           for the Y coordinates. Use the first existing
%                           column. Default is
%                           {'YWIN_IMAGE','Y','Y_IMAGE'}.
%            'OutUnits'   - Output units. Default is 'rad'.
% Output : - The input SIM object with an update RA/Dec in the catalog.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=update_coordinates(Cat);
% Reliable: 2
%--------------------------------------------------------------------------

CatField = AstCat.CatField;


DefV.ColNameRA            = 'ALPHAWIN_J2000';
DefV.ColNameDec           = 'DELTAWIN_J2000';
DefV.ColNameX             = {'XWIN_IMAGE','X','X_IMAGE'};
DefV.ColNameY             = {'YWIN_IMAGE','Y','Y_IMAGE'};
DefV.OutUnits             = 'rad';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Ncat = numel(Cat);
for Icat=1:1:Ncat
    W = ClassWCS.populate(Cat(Icat));
    
    [RA,Dec] = xy2coo(W,Cat(Icat),'ColNameX',InPar.ColNameX,'ColNameY',InPar.ColNameY,'OutUnits',InPar.OutUnits);
    ColRA  = colname2ind(Cat(Icat),InPar.ColNameRA);
    ColDec = colname2ind(Cat(Icat),InPar.ColNameDec);
    Cat(Icat).(CatField)(:,ColRA) = RA;
    Cat(Icat).(CatField)(:,ColDec) = Dec;
    
end