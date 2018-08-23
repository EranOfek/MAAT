function prep_atlas_htm(varargin)
% SHORT DESCRIPTION HERE
% Package: VO.prep
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.prep.prep_atlas_htm
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.CatName             = 'VSTatlas';
DefV.FileBaseName        = 'atlas';
DefV.FileExtName         = '.fit';
DefV.FileSplit           = '_';
DefV.FileType            = 'fits';
DefV.DecSize             = 1;  % deg
DefV.HTM_Level           = 9;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Files=dir(sprintf('%s*%s',InPar.FileBaseName,InPar.FileExtName));
Nf   = numel(Files);

RegExp = regexp({Files.name},InPar.FileSplit,'split');
Dec1   = nan(Nf,1);
Dec2   = nan(Nf,1);
for If=1:1:Nf
    StrDec1 = RegExp{If}{2};
    StrDec2 = RegExp{If}{3};
    StrDec2 = StrDec2(1:end-numel(InPar.FileExtName));
    
    if (strcmp(StrDec1(1),'m'))
        Dec1(If) = -str2double(StrDec1(2:end));
    else
        Dec1(If) = str2double(StrDec1);
    end
    if (strcmp(StrDec2(1),'m'))
        Dec2(If) = -str2double(StrDec2(2:end));
    else
        Dec2(If) = str2double(StrDec2);
    end
    
end

% sort by declination
[~,SI] = sort(Dec1);
Dec1   = Dec1(SI);
Dec2   = Dec2(SI);
Files  = Files(SI);

for If=1:1:Nf
    tic;
    clear Cat;
    
    Id = find(Dec1>(Dec1(If)-InPar.DecSize) & Dec2<(Dec2(If)+InPar.DecSize));
    Id = [Id; If-1; If+1];
    Id = unique(Id);
    Id = Id(Id>0 & Id<=Nf);
    
    Nid = numel(Id);
    for Iid=1:1:Nid
        Cat(Iid) = load_file(Files(Id(Iid)).name,InPar.FileType);
    end
    Cat = merge(Cat);
  
    % reorganize the catalog (specific)
    ColCell  = Cat.ColCell(3:end);
    ColUnits = Cat.ColUnits(3:end);
    Cat.Cat = Cat.Cat(:,3:end);
    Cat.Cat(:,1:2) = Cat.Cat(:,1:2)./RAD;
    Tmp = Cat.Cat(:,5:end);
    FF=Tmp>-1e-30 & Tmp<0;
    Tmp(FF) = NaN;
    Cat.Cat(:,5:end) = Tmp;
    
    ColCell{1} = 'RA';
    ColCell{2} = 'Dec';
    ColUnits{1} = 'rad';
    ColUnits{2} = 'rad';
    Cat         = sortrows(Cat,2);
    
    DecRange    = [Dec1(If), Dec2(If)]./RAD;
    VO.prep.build_htm_catalog(Cat.Cat,'HTM_Level',InPar.HTM_Level,'CatName',InPar.CatName,'SaveInd',false,'DecRange',DecRange);
    [If, toc]
end

% prep ind file
% save HTM index file
IndFileName = sprintf('%s_htm.hdf5',InPar.CatName);
delete(IndFileName);
Nsrc=HDF5.get_nsrc(InPar.CatName);
catsHTM.save_htm_ind(InPar.HTM_Level,IndFileName,[],ColCell,Nsrc)

catsHTM.save_cat_colcell(InPar.CatName,ColCell,ColUnits);


end


function Cat=load_file(File,Type)

switch lower(Type)
    case 'fits'
        Cat = FITS.read_table(File,'ModColName',true);
        Ncol = numel(Cat.ColCell);
        for Icol=1:1:Ncol
            Cat.ColCell{Icol}=Cat.ColCell{Icol}(2:end);
        end
    otherwise
        error('Unknown Type option');
end

end