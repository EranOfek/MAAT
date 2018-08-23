function [Cat,ColCell,Col]=get_cat(InCat,RA,Dec,varargin)
%--------------------------------------------------------------------------
% get_cat function                                               Catalogue
% Description: Search an astronomical catalog structure by coordinates.
% Input  : - Catalog. This can be a structure array or a file name string  
%            containing an astronomical catalog (e.g., 'FIRST.mat') or an
%            HTM index of a catalog (e.g., 'NVSS_htm.mat') or an array.
%            Alternatively, this can be a function handle that query
%            a catalog. The function input/output should looks like:
%            [Cat,ColCell,Col] = Fun(RA,Dec,SearchRad,Shape);
%            Examples for such functions include: @wget_sdss, @wget_2mass,
%            @wget_usnob1, @wget_ucac4.
%          - J2000.0 RA to search (rad, sexagesimal string, [H M S]).
%          - J2000.0 Dec to search (rad, sexagesimal string, [Sign D M S]).
%          * Search radius in radians. Default is 1 deg.
%            or a search radius followed by ...key,val,... arguments,
%            or a ...key,val,... arguments. The following keys are
%            available:
%            'SearchRad'- Search radius [radians]. Default is 1 deg.
%            'Shape'    - Search shape {'circ'}. Default is 'circ'.
%            'ColCell'  - Column names or indices of columns to return.
%            'CooOut'   - Units of output coordinates {'rad','deg'}.
%                         Default is 'rad'.
%            'StoreWS'  - Store catalog in workspace for fast retrieval
%                         in the next run. Default is true.
%            'ColRA'    - RA column name or index in catalog.
%                         Default is 'RA'.
%            'ColDec'   - Dec column name or index in catalog.
%                         'Default is 'Dec'.
%            'OutType'  - Output type. Options are:
%                         'mat' - return an array (default).
%                         'struct' - return a struct with .Cat, .Col
%                                 and .ColCell fields.
%            'Use_search_cat' - Use search_cat.m (true) or
%                         sphere_dist_fast.m (false) in the HTM index
%                         search. Default is false.
% Output : - Catalog.
%          - Cell array of returned column names.
%          - Structure of column indices.
% See also: get_first.m, get_nvss.m
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Cat,ColCell,Col]=get_cat('FIRST','12:00:00','+10:00:00');
%          [Cat]=get_cat('FIRST','12:00:00','+10:00:00',0.3./RAD,'OutType','struct');
%          [Cat]=get_cat('NVSS_htm','12:00:00','+10:00:00',0.3./RAD,'OutType','struct','ColCell',{'RA','Dec','Flux'});
%          Cat=get_cat('ROSAT_faint.mat',1,1);
%          [Cat]=get_cat(@wget_sdss,'12:00:00','+20:00:00',0.1./RAD,'OutType','struct');
% Reliable: 2
%--------------------------------------------------------------------------
import celestial.coo.*
import celestial.htm.*

RAD          = 180./pi;
InvRAD       = pi./180;
Field_isHTM  = 'isHTM';

DefV.SearchRad = InvRAD;
DefV.Shape     = 'circ';
DefV.ColCell   = [];  % return all columns
DefV.CooOut    = 'rad';
DefV.StoreWS   = true;
DefV.ColRA     = 'RA';    % or column index
DefV.ColDec    = 'Dec';
DefV.OutType   = 'mat';   % {'mat'|'struct'|'astcat'}
DefV.Use_search_cat = false;

% deal with the case that the first varargin is SearchRad
SearchRad = [];
if (numel(varargin).*0.5~=floor(numel(varargin).*0.5)),
    SearchRad = varargin{1};
    if (numel(varargin)>1),
        varargin  = varargin(2:end);
    else
        varargin  = {};
    end
end
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% override the SearchRad in varargin
if (isempty(SearchRad)),
    SearchRad = InPar.SearchRad;
end


if (isa(InCat,'function_handle')),
    % If InCat is a function handle then assume that
    % this function will retrieve catalog
    %error('InCat function is not supported yet');
    [Cat,ColCell,Col] = InCat(RA,Dec,SearchRad,InPar.Shape);
    % select columns
    [ColumnsInd,ColCell] = col_name2indvec(ColCell,InPar.ColCell,false);
    Cat = Cat(:,ColumnsInd);
    
else

    % input catalog
    InCatFN = [];  % InCat file name
    if (ischar(InCat)),
        % load catalog and store it in work space
        if (isempty(strfind(InCat,'.mat'))),
            % assume catalog is a mat file
            InCat = sprintf('%s.mat',InCat);
        end
        InCatFN = InCat;  % InCat file name
        InCat   = Util.IO.load_check(InCat,InPar.StoreWS);
    end

    if (isstruct(InCat)),
        % the input is a catalog structure
        [~,InCat.Col.RA,InCat.Col.Dec] = col_name2ind(InCat.ColCell, InPar.ColRA, InPar.ColDec);

    elseif (isnumeric(InCat)),
        InCat.Cat      = InCat;
        if (ischar(InPar.ColRA) || ischar(InPar.ColDec)),
            error('If catalog is numeric then column indices must be numeric');
        end
        InCat.Col.RA   = InPar.ColRA;   % in this case ColRA must be an index
        InCat.Col.Dec  = InPar.ColDec;  % in this case ColDec must be an index
    else
        error('Unsupported catalog format');
    end

    % check if input catalog is sorted by declination
    if (~issorted(InCat.Cat(:,InCat.Col.Dec))),
        error('Input catalog must be sorted by declination');
    end


    % check if catalog is an HTM catalog
    if (isfield_check(InCat,Field_isHTM,@all)),
        % deal with HTM catalog    

        HTM_BASE_NAME = InCat.FileBaseName; %sprintf('%s%%05d.mat',InCatFN(1:end-4));
        DataDir       = regexp(InCatFN,'.mat','split');
        DataDir       = DataDir{1};

        % convert coordinates to radians
        RA  = convertdms(RA,'gH','r');
        Dec = convertdms(Dec,'gD','R');

        DataPath = sprintf('%s%s%s%s',Util.files.which_dir(InCatFN),filesep,DataDir,filesep);

        % pointers to HTMs within search region
        HTMind = search_htm_coocat(RA,Dec,SearchRad,InCat,InPar.Use_search_cat);
        % construct HTM file names
        Files  = Util.cell.sprintf2cell(HTM_BASE_NAME,HTMind);
        % load files
        Cat.Cat = zeros(0,length(InCat.CatHTM.ColCell));
        for Ifiles=1:1:numel(Files),
            Cat1 = Util.IO.load2(sprintf('%s%s',DataPath,Files{Ifiles}));
            Cat.Cat = [Cat.Cat; Cat1];
        end

        % select columns
        Cat.ColCell          = InCat.CatHTM.ColCell;
        Cat.Col              = InCat.CatHTM.Col;
        Cat.ColUnits         = InCat.CatHTM.ColUnits;
        Cat.SortedBy         = InCat.CatHTM.SortedBy;
        Cat.SortedByCol      = InCat.CatHTM.SortedByCol;
        
        [ColumnsInd,ColCell] = col_name2indvec(Cat.ColCell,InPar.ColCell,false);

        % return sources within search radius
        Col = InCat.Col;
        D   = sphere_dist_fast(Cat.Cat(:,Cat.Col.RA),Cat.Cat(:,Cat.Col.Dec),RA,Dec);
        Cat = Cat.Cat(D<SearchRad,ColumnsInd);
        Cat = sortrows(Cat,Col.Dec);

    else

        % select columns
        [ColumnsInd,ColCell] = col_name2indvec(InCat.ColCell,InPar.ColCell,false);

        % deal with normal catalog
        Res = search_cat(InCat.Cat(:,[InCat.Col.RA, InCat.Col.Dec]),RA,Dec,...
                         'SearchRad',SearchRad,'SearchMethod','binms','CooType','sphere');
        Cat = InCat.Cat(Res.IndCat,ColumnsInd);
        Col = cell2struct(num2cell(1:1:length(ColCell)),ColCell,2);
    end


end

switch lower(InPar.CooOut)
    case 'rad'
        % do nothing
    case 'deg'
        Cat(:,[Col.RA,Col.Dec]) = Cat(:,[Col.RA,Col.Dec]).*RAD;
    otherwise
        error('Unknown CooOut option');
end



switch lower(InPar.OutType)
    case 'mat'
        % do nothing
    case 'struct'
        CatN = Cat;
        clear Cat;
        Cat.Cat     = CatN;
        Cat.Col     = Col;
        Cat.ColCell = ColCell;
    case 'astcat'
        CatN = Cat;
        clear Cat;
        Cat = AstCat;
        Cat.Cat = CatN;
        Cat.Col     = Col;
        Cat.ColCell = ColCell;
        
    otherwise
        error('Unknown OutType option');
end



