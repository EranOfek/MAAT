function Prop=irsa_table2prop(Cat,varargin)
% Table of queried ZTF image to properties structure
% Package: VO.ZTF
% Description: Table of queried ZTF image to properties structure.
% Input  : - Catalog (output of VO.ZTF.query_ztf_images)
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ImType'  - Image type. Default is 'sci'.
%            'Product' - Product type. Default is 'sci'.
% Output : - Structure array of images.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Prop=VO.ZTF.irsa_table2prop(T)
% Reliable: 2
%--------------------------------------------------------------------------

CatField = AstCat.CatField;

DefV.ImType               = 'sci';
DefV.Product              = 'image';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Date  = celestial.time.jd2date(Cat.(CatField).obsjd);
Field = Cat.(CatField).field;
Ncat = numel(Field);
Prop = Util.struct.struct_def({'ImType','Year','Month','Day','FracDay','FieldID','FilterCode','CCDID','ImTypeCode','QuadID','Product'},Ncat,1);

for Icat=1:1:Ncat
    Prop(Icat).ImType     = InPar.ImType;
    Prop(Icat).Year       = Date(Icat,3);
    Prop(Icat).Month      = Date(Icat,2);
    Prop(Icat).Day        = Date(Icat,1);
    Prop(Icat).FracDay    = Cat.(CatField).filefracday(Icat) - str2double(sprintf('%04d%02d%02d000000',Prop(Icat).Year,Prop(Icat).Month, Prop(Icat).Day));
    Prop(Icat).FieldID    = Cat.(CatField).field(Icat);
    Prop(Icat).FilterCode = Cat.(CatField).filtercode{Icat};
    Prop(Icat).CCDID      = Cat.(CatField).ccdid(Icat);
    Prop(Icat).ImTypeCode = Cat.(CatField).imgtypecode{Icat};
    Prop(Icat).QuadID     = Cat.(CatField).qid(Icat);
    Prop(Icat).Product    = InPar.Product;
end