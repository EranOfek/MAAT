function [Cat,ColCell]=htmcat_search(CatName,RA,Dec,Radius,varargin)
% Cone earch on local HDF5/HTM catalog (OBSOLETE).
% Package: VO.search
% Description: Perform a cone search around RA/Dec on a local catalog in
%              HDF5 format sorted into HTM.
%              OBSOLETE: use catsHTM.cone_search instead.
% Input  : - Catalog name (e.g., 'GAfADR1').
%            see VO.search.htmcat_names for options.
%          - J2000.0 R.A. [radians, [H M S], or sexagesimal string].
%          - J2000.0 Dec. [radians, [sign D M S], or sexagesimal string].
%          - Search radius [arcsec].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Dec 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=VO.search.htmcat_search('UCAC4',1,1,10);
%          Cat=VO.search.htmcat_search('GAIADR1',1,1,10);
%          Cat=VO.search.htmcat_search('GALEX',1,1,10);
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.RadiusUnits          = 'arcsec';  % do not change this default!
DefV.IndexFileTemplate    = '%s_htm.hdf5';
DefV.CatFileTemplate      = '%s_htm_%06d.hdf5';
DefV.htmTemplate          = 'htm_%06d';
DefV.NcatInFile           = 100;
DefV.IndexVarName         = [];
DefV.ColRA                = 1;
DefV.ColDec               = 2;
DefV.OnlyCone             = true;
DefV.ColCellFile          = '%s_htmColCell.mat';
DefV.OutType              = 'mat';

if (isempty(varargin))
    InPar  = DefV;
    Radius = Radius./(RAD.*3600);  % arcsec to [radians]
else
    InPar  = InArg.populate_keyval(DefV,varargin,mfilename);
    Radius = convert.angular(InPar.RadiusUnits,'rad',Radius);  % [radians]
end

InPar.ColCellFile = sprintf(InPar.ColCellFile,CatName);

%Ncol = 45;
%Ncol = 8;
load(InPar.ColCellFile);
Ncol  = numel(ColCell);


MinDec = Dec - Radius;
MaxDec = Dec + Radius;

IndexFileName = sprintf(InPar.IndexFileTemplate,CatName);
ID     = catsHTM.search_htm_ind(IndexFileName,InPar.IndexVarName,RA,Dec,Radius);
FileID = floor(ID./InPar.NcatInFile).*InPar.NcatInFile;
Nid = numel(ID);
Cat = zeros(0,Ncol);
%C = Util.struct.struct_def({'Cat'},Nid,1);
for Iid=1:1:Nid
    
    %FileID    = floor(ID(Iid)./InPar.NcatInFile).*InPar.NcatInFile;
    FileName  = sprintf(InPar.CatFileTemplate,CatName,FileID(Iid));
    DataName  = sprintf(InPar.htmTemplate,ID(Iid));

    %Cat = [Cat; HDF5.load_cat(FileName,DataName,[MinDec, MaxDec],Ncol)];
    if (Iid==1)
        Cat = HDF5.load(FileName,DataName);
        %Cat = HDF5.load_cat(FileName,DataName,[MinDec, MaxDec],Ncol);
        %Ncol = size(Cat,2);
    else
        Cat = [Cat; HDF5.load(FileName,DataName)];
        %Cat = [Cat; HDF5.load_cat(FileName,DataName,[MinDec, MaxDec],Ncol)];
    end
    
    %C(Iid).Cat = HDF5.load_cat(FileName,DataName,[MinDec, MaxDec],Ncol).';
    
end
    
%Cat = [C.Cat]';
    
if (InPar.OnlyCone)
    D = celestial.coo.sphere_dist_fast(RA,Dec,Cat(:,InPar.ColRA),Cat(:,InPar.ColDec));
    Cat = Cat(D<Radius,:);
end



switch lower(InPar.OutType)
    case 'mat'
        % do nothing
    case 'astcat'
        AstC = AstCat;
        AstC.Cat = Cat;
        AstC.ColCell = ColCell;
        AstC = colcell2col(AstC);
        Cat  = AstC;
    otherwise
        error('Unknown OutType option');
end

