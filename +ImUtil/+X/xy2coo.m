function [RA,Dec]=xy2coo(File,X,Y,Tel)
% X-Ray event file X/Y to J2000 equatorial coordinates.
% Package: AstroX
% Description: Given a Swift XRT event file (image) convert X/Y (image)
%              position to J2000.0 equatorial coordinates.
%              see also: swxrt_coo2xy.m
% Input  : - FITS binary table of events file name.
%          - X coordinates. Default is to use all the X coordinates in
%            the FITS events file name.
%          - Y coordinates. Default is to use all the Y coordinates in
%            the FITS events file name.
%          - {'swift' | 'chandra' | 'nustar'}. Default is 'swift'.
% Output : - J2000.0 R.A. [radians].
%          - J2000.0 Dec. [radians].
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Nov 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%------------------------------------------------------------------------------

RAD = 180./pi;

Def.X = [];
Def.Y = [];
Def.Tel = 'swift';
if (nargin==1)
    X = Def.X;
    Y = Def.Y;
    Tel = Def.Tel;
elseif (nargin==2)
    Y = Def.Y;
    Tel = Def.Tel;
elseif (nargin==3)
    Tel = Def.Tel;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments');
end

switch lower(Tel)
 case 'swift'
    Keys       = {'TCRPX2','TCRVL2','TCDLT2','TCRPX3','TCRVL3','TCDLT3'};
 case 'chandra'
    Keys       = {'TCRPX11','TCRVL11','TCDLT11','TCRPX12','TCRVL12','TCDLT12'};
 case 'nustar'
    Keys       = {'TCRPX37','TCRVL37','TCDLT37','TCRPX38','TCRVL38','TCDLT38'};
 otherwise
    error('Unsupported telescope option');
end
%[KeywordVal,KS] = get_fits_keyword(File,Keys,1,'BinaryTable');
[KeywordVal,~,~] = FITS.get_keys(File,Keys,2); %,'BinaryTable');

KS.TCRPX2  = str2double(KeywordVal{1});
KS.TCRVL2  = str2double(KeywordVal{2});
KS.TCDLT2  = str2double(KeywordVal{3});
KS.TCRPX3  = str2double(KeywordVal{4});
KS.TCRVL3  = str2double(KeywordVal{5});
KS.TCDLT3  = str2double(KeywordVal{6});


if (isempty(X) || isempty(Y))
    %[~,~,~,Col]=get_fitstable_col(File,'BinTable');
    TT = FITS.read_table(File);
    Col = TT.Col;
    

    Table = fitsread(File,'BinTable');
    switch lower(Tel)
     case {'swift','nustar'}
        X = Table{Col.X};
        Y = Table{Col.Y};
     case 'chandra'
        X = Table{Col.x};
        Y = Table{Col.y};
     otherwise
        error('Unsupported telescope option');
    end
end


DX = (X - KS.TCRPX2).*KS.TCDLT2./RAD;
DY = (Y - KS.TCRPX3).*KS.TCDLT3./RAD;
[RA,Dec] = celestial.proj.pr_ignomonic(DX,DY,[KS.TCRVL2, KS.TCRVL3]./RAD);

