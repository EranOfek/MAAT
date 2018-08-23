function [AstC,ColCell,Col,ColUnits]=get_ucac4(RA,Dec,varargin)
%--------------------------------------------------------------------------
% get_ucac4 function                                             Catalogue
% Description: Search the local UCAC4 catalog (in HTM format).
% Input  : - J2000.0 RA [rad].
%          - J2000.0 Dec [rad].
%          - Search radius [rad]. Default is 1 deg.
%          - Search shape 'circ'|'box'. Default is 'circ'.
%          - Catalog out type: 'mat|'astcat'. Default is 'mat'.
% Output : - Catalog of UCAC4 sources within search region.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AstC]=get_ucac4(RA,Dec,1./RAD)
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

CatField = 'Cat';

NumVarargs = numel(varargin);
if NumVarargs > 3
     errId = 'get_ucac4.m:TooManyInputArguments';
     errMsg = 'RA, Dec, [SearchRad, Shape, OutType]';
     error(errId, errMsg);
end
Gaps = cellfun(@isempty, varargin);
DefArgs = {1./RAD, 'circ', 'mat'};    % default input arguments
Suboptargs = DefArgs(1 : NumVarargs);
varargin(Gaps) = Suboptargs(Gaps);
DefArgs(1 : NumVarargs) = varargin;
[SearchRad, Shape, OutType] = DefArgs{:};


FileBaseName = 'UCAC4htm_';
VarBaseName  = 'UCAC4_htm';
Suffix       = 'hdf5';
HTMcatName   = 'UCAC4htm.mat';


ColCell = {'RA','Dec','MagModel','MagAper','ErrMag','ObjType','DSflag','ErrRA','ErrDec',...
           'Nim','Nobs','Nep','EpochRA','EpochDec','PM_RA','PM_Dec','ErrPM_RA','ErrPM_Dec',...
           'ID2MASS','MagJ','MagH','MagK','FlagJ','FlagH','FlagK','ErrMagJ','ErrMagH','ErrMagK',...
           'MagB','MagV','Magg','Magr','Magi','ErrMagB','ErrMagV','ErrMagg','ErrMagr','ErrMagi',...
           'FlagYale','FlagExtCat','FlagLEDA','Flag2MASSext','ID','ZoneUCAC2','RecordUCAC2'};
ColUnits = {'rad','rad','mag','mag','mag','','','mas','mas','','','','year','year','mas/yr','mas/yr','mas/yr','mas/yr',...
            '','mag','mag','mag','','','','mag','mag','mag','mag','mag','mag','mag','mag','mag','mag','mag','mag','mag',...
            '','','','','','',''};
Ncol     = numel(ColCell);
Col      = cell2struct(num2cell(1:1:Ncol),ColCell,2);
        
        



UCAC4htm = Util.IO.load_check(HTMcatName);

Ncen = size(UCAC4htm.(CatField),1);  % number of HTM in lower level
Nlev = round(log(Ncen.*0.5)./log(4));  % number of levels
TriSide = 180./(2.^Nlev);              % estimated size of the triangle side
TriCM   = TriSide.*(sqrt(2)./2).*1.01./RAD;         % estimated size of distance from tri edge to center of mass [rad]



ColH    = UCAC4htm.Col;
Dist   = celestial.coo.sphere_dist(RA,Dec,UCAC4htm.Cat(:,ColH.RA),UCAC4htm.Cat(:,ColH.Dec));
Flag   = Dist<(SearchRad+TriCM);
SubCat = UCAC4htm.Cat(Flag,:);
Nhtm   = size(SubCat,1);

%HTMind = search_htm_coocat(RA,Dec,SearchRad,UCAC4htm,false);
%Nhtm = numel(HTMind);
Cat  = zeros(0,Ncol);
for Ihtm=1:1:Nhtm
    Line = UCAC4htm.(CatField)(UCAC4htm.Cat(:,ColH.Ptr)==SubCat(Ihtm,ColH.Ptr),:);
    
    FileName = sprintf('%s%04d.%s',FileBaseName,Line(3),Suffix);
    VarName  = sprintf('%s%07d',VarBaseName,Line(4));
    Cat      = [Cat; Util.IO.loadh(FileName,VarName)];
end


[Flag,~,~] = celestial.coo.sphere_dist_thresh(Cat(:,1),Cat(:,2),RA,Dec,SearchRad,Shape);
Cat  = Cat(Flag,:);

% switch lower(Shape)
%     case 'circ'
%         Dist = sphere_dist(RA,Dec,Cat(:,1),Cat(:,2));
%         Cat  = Cat(Dist<SearchRad,:);
%         
%     
%     otherwise
%         error('Unknown Shape option');
% end

switch lower(OutType)
    case 'astcat'
        AstC            = AstCat;
        AstC.(CatField) = Cat;
        AstC.ColCell    = ColCell;
        AstC            = colcell2col(AstC);
        AstC.ColUnits   = ColUnits;
        AstC.Name       = 'UCAC4 local catalog';
    case 'mat'
        AstC = Cat;
    otherwise
        error('Unknown OutType option');
end


    
