function [AstOut,AstUM]=match(AstC,AstR,varargin)
% Match catalogs in an AstCat object against a reference catalog
% Package: @AstCat
% Description: Match a collection of catalogs in an AstCat object against
%              a reference catalog. The matching is done either by X,Y
%              coordinates or spherical coordinates (e.g., RA,Dec).
%              The output is an AstCat object (corresponding to the
%              input catalogs) in which the number of lines in each output
%              catalog is like the number of lines in the reference
%              catalogs, and these lines corresponds (in the same order as)
%              to the reference catalog.
%              Lines in the catalog that don't have a match in the
%              reference catalog are populated with NaN.
%              The columns in the output catalog contains any requested
%              subset of columns of the input catalog and the reference
%              catalog. In addition, several other columns like the number
%              of matches, match distance, and distance from image edge are
%              provided.
%              In addition, the program returns an AstCat object that
%              contains the sources in the input catalog without a match in
%              the reference catalog.
% Input  : - An AstCat object containing multiple catalogs to match against
%            the reference catalog. The catalog must contain an a 2D
%            coordinates in pix, radians or degres units.
%            By default, this catalog will be sorted by Y or declination.
%          - An AstCat object containing a single reference catalog against
%            all the the input catalogs will be matched.
%            The catalog must contain an a 2D
%            coordinates in pix, radiaons or degres units.
%            If empty or not provided, one of the catalogs in the first
%            input argument will be selected as a reference catalog.
%            This selection is done using the 'RefSelectMetod' option (see
%            below). Defult is empty.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'RefSelectMetod' - In case the reference catalog is not
%                        provided it will be selected from the first input
%                        catalog. This option specify how to select the
%                        reference catalog from the first input catalog.
%                        Options are:
%                        'first' - First catalog.
%                        'last'  - Last catalog.
%                        'max'   - Catalog with maximum number of lines.
%                        Default is 'max'.
%           'SkipWCS'     - If true than will change CatColNames,
%                         RefColNames, CatUnits, RefUnits, CooType,
%                         SearchRadUnits to start attempt to match using
%                         X/Y coordinates rather than RA/Dec.
%                         Default is false.
%           'CatColNames' - Either a vector of the column indices in the
%                        catalog containing the 2D coordinates, or a row
%                        cell array of column names (e.g., {'X','Y'}).
%                        If multiple rows are given then the program will
%                        check if the column names exist in the catalog and
%                        will use the first column names found in the
%                        catalog. Default is"
%                        {'ALPHAWIN_J2000','DELTAWIN_J2000';...
%                         'RA','Dec';...
%                         'XWIN_IMAGE','YWIN_IMAGE';...
%                         'X_IMAGE','Y_IMAGE';...
%                         'X','Y'}.
%                        This means that by default the program will try to
%                        match the coordinates using the 'ALPHAWIN_J2000'
%                        and 'DELTAWIN_J2000' columns.
%           'RefColNames' - Like 'CatColNames' but for the reference catalog.
%                        If NaN then will use the content of 'CatColNames'.
%                        Default is NaN.
%           'CatUnits'    - A string indicating the coordinate units in the
%                        first input catalog. If a cell array of strings
%                        then each element corresponds to each line in the
%                        'CatColNames' input. For example, if the 2nd row
%                        of 'CatColNames' is selected than the second
%                        element of 'CatUnits' will be used.
%                        Default is {'deg'; 'deg'; 'pix';'pix'; 'pix'}.
%           'RefUnits'    - Like 'CatUnits' but for the reference catalog.
%                        If NaN then will use the content of 'CatUnits'.
%                        Default is NaN.
%           'CooType'     - A string indicating the coordinate type
%                        ('sphere' | 'plane'). If a cell array of strings
%                        then each element corresponds to each line in the
%                        'CatColNames' input. For example, if the 2nd row
%                        of 'CatColNames' is selected than the second
%                        element of 'CooType' will be used.
%                        Default is {'sphere';'sphere';'plane';'plane';'plane'}.
%           'SearchRad'   - Search radius. Default is 2.
%           'SearchRadUnits'- Units of search radius.
%                        A string indicating the search radius units
%                        ('arcsec'|'rad'|'deg'|'pix').
%                        If a cell array of strings
%                        then each element corresponds to each line in the
%                        'CatColNames' input. For example, if the 2nd row
%                        of 'CatColNames' is selected than the second
%                        element of 'SearchRadUnits' will be used.
%                        Default is {'arcsec';'arcsec';'pix';'pix';'pix'}.
%           'SortCat'     - Sort input catalogs by Y/Dec. Default is false.
%                        However this may be override if 'CheckSortCat' is
%                        true.
%           'SortRef'     - Sort reference catalog by Y/Dec.
%                        Default is false.
%           'CheckSortCat' - Sort the input catalog if it is not sorted.
%                        This overrides the 'SortCat' option.
%                        Default is true.
%           'SearchMethod' - Search method. See search_cat.m for options.
%                        Default is 'binms'.
%                        This option requires that the input catalog will
%                        be sorted by Y/Dec.
%           'CatOutCol'    - Indices or a cell array of column names in the
%                        input catalog to include in the output catalog.
%                        If empty then include all columns.
%                        If NaN then no columns.
%                        Default is [].
%           'RefOutCol'    - Like 'CatOutCol' but for the reference
%                        catalog. Default is NaN.
%           'RefColPrefix' - Prefix for reference catalog column names
%                        appended to the output catalog. Default is 'REF_'.
%           'NmatchName'   - Name of column indicating the number of
%                        matches. Default is 'N_Match'.
%           'DmatchName'   - Name of column indicating the distance between
%                        matches. Default is 'Dist_Match'.
%           'DistUnits'    - Match distance units ('arcsec'|'deg'|'rad').
%                        Default is 'arcsec'.
%           'AuxCol'       - Vector of flags indicating if to include in
%                        the output catalog the [Nfound, Dist] columns.
%                        Default is [true true].
%           'EdgeCol'      - Vector of flags indicating if to include in
%                        the output catalog the
%                        [Flag, DistEdge, DistCenter, DistPA] columns.
%                        Default is [true true true true].
%                        Note that these columns remains NaN in the AstOut
%                        and are populated only in AstUM.
%           'Verbose'      - Verbose. Default is false.
% Output : - AstCat object with the combined matched catalogs.
%            The number of lines is the number of sources in the reference
%            catalog. If there is no match than Nfound will be zero, and
%            the columns of the other catalog will be set to NaN.
%          - An AstCat object with sources in the other catalogs that were
%            unmatched to the reference catalog.
%            The number of elements in this object is equal to the number
%            of other catalogs.
%            Also appended to the columns in these catalogs is a flag
%            indicating if the sources is in the footprint of the
%            reference catalog and its distance from the image edge [pix].
% See also: astcat2matched_array.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [AstOut,AstUM]=match(S1s,S1m);
%          [AstOut,AstUM]=match(S1s,S1m,'CooType','sphere');
%          [AstOut,AstUM]=match(S1s,S1m,'CatOutCol',{'XWIN_IMAGE','YWIN_IMAGE'},'RefOutCol',{'XWIN_IMAGE','YWIN_IMAGE'});
%          [AstOut,AstUM]=match(S1s,S1m,'CatOutCol',{'XWIN_IMAGE','YWIN_IMAGE'},'RefOutCol',NaN,'EdgeCol',false(1,4));
% Reliable: 2
%--------------------------------------------------------------------------
RAD        = 180./pi;
ARCSEC_DEG = 3600;

CatField       = 'Cat';
ColCellField   = 'ColCell';

if (nargin==1)
    AstR = [];
end

% default parameters
DefV.RefSelectMethod   = 'max';   % 'first'|'last'|'max'
DefV.SkipWCS           = false;
DefV.CatColNames       = {'ALPHAWIN_J2000','DELTAWIN_J2000';...
                          'RA','Dec';...
                          'XWIN_IMAGE','YWIN_IMAGE';...
                          'X_IMAGE','Y_IMAGE';...
                          'X','Y';
                          'xpos','ypos'}; 
DefV.RefColNames       = NaN;   % NaN means copy from RefColNames.
DefV.CatUnits          = {'deg'; 'deg'; 'pix'; 'pix'; 'pix'; 'pix'};
DefV.RefUnits          = NaN;
DefV.CooType           = {'sphere';'sphere';'plane';'plane';'plane';'plane'};
DefV.SearchRad         = 2;
DefV.SearchRadUnits    = {'arcsec';'arcsec';'pix';'pix';'pix';'pix'};
DefV.SortCat           = false;
DefV.SortRef           = false;
DefV.CheckSortCat      = true;
DefV.SearchMethod      = 'binms';
% Output columns
DefV.CatOutCol         = [];
DefV.RefOutCol         = NaN;
DefV.RefColPrefix      = 'REF_';
DefV.NmatchName        = 'N_Match';
DefV.DmatchName        = 'Dist_Match';
DefV.DistUnits         = 'arcsec';
DefV.AuxCol            = [true true];   % show Aux columns [Nfound, Dist]
DefV.EdgeCol           = [true true true true]; % show dis edge columns [Flag, DistEdge, DistCenter, DistPA]
DefV.Verbose           = false;

%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

NauxCol  = sum(InPar.AuxCol);
NedgeCol = sum(InPar.EdgeCol);

if (InPar.SkipWCS)
    InPar.CatColNames = {'XWIN_IMAGE','YWIN_IMAGE';...
                         'X_IMAGE','Y_IMAGE';...
                         'X','Y';
                         'xpos','ypos'}; 
    
    InPar.RefColNames = NaN;
    InPar.CatUnits    = {'pix'; 'pix'; 'pix'; 'pix'};
    InPar.RefUnits    = NaN;
    InPar.CooType     = {'plane';'plane';'plane';'plane'};
    InPar.SearchRadUnits = {'pix';'pix';'pix';'pix'};
end

if (~iscell(InPar.CooType))
    InPar.CooType = {InPar.CooType};
end

% If Ref cat is not provided use one of the elements of AstC
if (isempty(AstR))
    switch lower(InPar.RefSelectMethod)
        case 'first'
            AstR = AstC(1);
        case 'last'
            AstR = AstC(end);
        case 'max'
            [~,IndMax] = max(sizecat(AstC));
            AstR = AstC(IndMax);
        otherwise
            error('Unknown RefSelectMethod');
    end
end

if (numel(AstR)>1)
    error('Reference catalog should contain a single element');
end

% Make sure that AstR is an AstCat object.
% If needed convert to AstCat object
if (~AstCat.isastcat(AstR))
    Tmp  = AstR;
    AstR = AstCat;
    AstR.(CatField)     = Tmp;
    % assumes the columns are those in InPar.RefColNames{1,:}
    AstR.(ColCellField) = InPar.RefColNames(1,:);
    AstR                = colcell2col(AstR);
end

% if (~iscell(InPar.CatName)),
%     InPar.CatName = {InPar.CatName};
% end
% Ncn = numel(InPar.CatName);


% If RefColNames is NaN - copy content from CatColName
if (~iscell(InPar.RefColNames))
    if (isnan(InPar.RefColNames))
        InPar.RefColNames = InPar.CatColNames;
    end
end
if (~iscell(InPar.RefUnits))
    if (isnan(InPar.RefUnits))
        InPar.RefUnits = InPar.CatUnits;
    end
end



% Identify the coordinates column names in AstR:
% The user can supply a list of several pairs of columns (X/Y) names,
% the function searched for the first pair that appears in the table and
% use it.
[RefColNames,RefColInd,RefUseInd]=select_colnames(InPar.RefColNames,AstR);

% select the column units for AstR
RefUnits = get_colunits(AstR,RefColInd(1));
if (isempty(RefUnits))
    Nru      = numel(InPar.RefUnits);
    if (~iscell(InPar.RefUnits))
        InPar.RefUnits = {InPar.RefUnits};
    end
    RefUnits = InPar.RefUnits{min(Nru,RefUseInd)};
end

% Search radius
if (~iscell(InPar.SearchRadUnits))
    InPar.SearchRadUnits = {InPar.SearchRadUnits};
end
Nrsu = numel(InPar.SearchRadUnits);
InPar.SearchRadUnits = InPar.SearchRadUnits{min(Nrsu,RefUseInd)};

% convert the search radius to radians:
% If 'pix' - then will do nothing
ConvSearchRad = convert.angular(InPar.SearchRadUnits,'rad',1,'pix');
InPar.SearchRad = InPar.SearchRad.* ConvSearchRad;




% sort Ref
if (InPar.SortRef)
    AstR = sortrows(AstR,RefColInd(2));
end
[RefNrow, RefNcol] = size(AstR.(CatField));

% Set CooType according to the units of RefUnits or
% the default CooType
[RefFactor,InPar.CooType] = set_cootype(RefUnits,InPar.CooType{RefUseInd});

% read the coordinates in the Ref and convert them to radians
RefCoo    = AstR.(CatField)(:,RefColInd) .* RefFactor;   % [rad] - [RA, Dec]  or [X,Y]


% Set the columns of Ref that will appended into the columns of Cat
% Options are:
%   NaN - no columns
%   []  - all columns
%   Cell array of column names, or vector of column indices
RefOutCol = appended_columns_ind(InPar.RefOutCol,RefNcol,AstR);
NoutColRef = numel(RefOutCol);

switch lower(InPar.CooType)
    case 'sphere'
        % convert radians to DistUnits:
        ConvDistUnits = convert.angular('rad',InPar.DistUnits);
    case 'plane'
        ConvDistUnits = 1;
    otherwise
        error('Unknown CooType option');
end
       


% if (isempty(InPar.RefOutCol)),
%     % Append all the columns of Ref into the output Cat
%     RefOutCol = (1:1:RefNcol);
% else  
%     if (~iscell(InPar.RefOutCol)),
%         % RefOutCol is a vector
%         if (isnan(InPar.RefOutCol)),
%             % Don't append Columns of Ref to the output
%             RefOutCol = [];
%         else
%             % InPar.RefOutCol is already a vector of indices
%             RefOutCol = InPar.RefOutCol;
%         end
%     else
%         % Convert InPar.RefOutCol to a vector of indices
%         RefOutCol = colname2ind(AstR,InPar.RefOutCol);
%     end
% end


% Is input catalog also SIM objects
IsSIMref = SIM.issim(AstR);
IsSIMcat = SIM.issim(AstC);

% Auxilary columns [Nfound, Dist] in putput catalog
% which columns to show
AuxCol          = [1 2];
AuxCol          = AuxCol(InPar.AuxCol);


% % for each catalog
% Nc = numel(AstC);
% if (InPar.MatchIndiv),
%     AstOut = astcatdef(size(AstC));
% else
%     AstOut = AstCat;
% end

AstUM  = AstCat(size(AstC));
AstOut = AstCat(size(AstC));
Nc     = numel(AstC);
for Ic=1:1:Nc
    if (InPar.Verbose)
        fprintf('Match catalog %d\n',Ic);
    end
    
    % Identify the coordinates column names in AstC:
    % The user can supply a list of several pairs of columns (X/Y) names,
    % the function searched for the first pair that appears in the table and
    % use it.
    [CatColNames,CatColInd,CatUseInd]=select_colnames(InPar.CatColNames,AstC(Ic));

    % select the column units for AstC
    CatUnits = get_colunits(AstC(Ic),CatColInd(1));
    if (isempty(CatUnits))
        Ncu      = numel(InPar.CatUnits);
        CatUnits = InPar.CatUnits{min(Ncu,CatUseInd)};
    end
    
    % Check if CooType is consistent with units
    [CatFactor,CatCooType] = set_cootype(CatUnits,InPar.CooType);

    if (~strcmp(InPar.CooType,CatCooType))
        error('Ref and Cat Coo types are not consistent');
    end
        
    % read the coordinates in the Cat and convert them to radians
    CatCoo    = AstC(Ic).(CatField)(:,CatColInd) .* CatFactor;   % [rad] - [RA, Dec]  or [X,Y]

    % Cat size
    [CatNrow, CatNcol] = size(AstC(Ic).(CatField));
    
    % catalog may need conversion to radians
    % RefCoo is already in radians (if CooType='sphere')
    if (InPar.SortCat)
        AstC(Ic) = sortrows(AstC(Ic),CatColInd(2));
    else
        if (InPar.CheckSortCat)
            if (~issorted(CatCoo(:,2)))
                AstC(Ic) = sortrows(AstC(Ic),CatColInd(2));
            end
        end
    end
    clear CatCoo;
    
    % Match the two catalogs
    % RefCoo is a two column matrix (already in radians if needed)
    [Res,CatUM]     = VO.search.search_cat(AstC(Ic).(CatField)(:,CatColInd).*CatFactor,RefCoo,[],...
                                 'CooType',InPar.CooType,'IsRad',true,...
                                 'SearchMethod',InPar.SearchMethod,...
                                 'SearchRad',InPar.SearchRad);
    % CatUM is a flag indicating sources in AstC which were not matched
    IndRef = logical([Res.Nfound])';
                             
    % get the indices of the columns in AstC to keep
    CatOutCol  = appended_columns_ind(InPar.CatOutCol,CatNcol,AstC(Ic));
    NoutColCat = numel(CatOutCol);
    
    NcolOut = NoutColCat + NoutColRef + NauxCol + NedgeCol;
    
    IndCat = nan(RefNrow,1);
    Dist   = nan(RefNrow,1);
    for Ires=1:1:RefNrow
        if (~isempty(Res(Ires).IndCat))
            IndCat(Ires) = Res(Ires).IndCat(1);
            Dist(Ires)   = Res(Ires).DistRAD(1);
        end
    end 
    
    %----------------------------------------------
    %--- Number of lines in AstOut like in AstR ---
    %----------------------------------------------
    % The lines of Astout corresponds to the lines of AstR.
    % User requested that the output catalog will have the same number
    % of lines as in the Ref.

    % Initiate the AstOut with NaN
    AstOut(Ic).(CatField) = nan(RefNrow,NcolOut);
    AstOutCol             = cell(NcolOut,1);
    % copy the requested columns (CatOutCol) in AstC to the appropraite
    % place in AstOut.
    % The lines of AstCat will corresponds to the lines of AstR
    AstOut(Ic).(CatField)(IndRef,1:NoutColCat) = AstC(Ic).(CatField)(IndCat(IndRef),CatOutCol);
    AstOutCol(1:NoutColCat) = AstC(Ic).(ColCellField)(CatOutCol);

    % append the requested AstR columns
    AstOut(Ic).(CatField)(:,NoutColCat+1:NoutColCat+NoutColRef) = AstR.(CatField)(:,RefOutCol);
    AstOutCol(NoutColCat+1:NoutColCat+NoutColRef) = Util.cell.cellstr_prefix(AstR.(ColCellField)(RefOutCol), InPar.RefColPrefix);


    % append the additional meta adat to the AstOut combined catalog
    CurrentCol = NoutColCat+NoutColRef;
    if (InPar.AuxCol(1))
        CurrentCol = CurrentCol + 1;
        AstOut(Ic).(CatField)(:,CurrentCol) = [Res.Nfound].';
        AstOutCol{CurrentCol} = InPar.NmatchName;
    end
    if (InPar.AuxCol(1))
        CurrentCol = CurrentCol + 1;
        AstOut(Ic).(CatField)(:,CurrentCol) = Dist.*ConvDistUnits;
        AstOutCol{CurrentCol} = InPar.DmatchName;
    end

    if (IsSIMref && any(InPar.EdgeCol))        
        % Only for unmatched sources!
        % Flag sources in Cat which are outside the boundries of Ref
        [InAstC]=coo_inimage(AstR,AstC(Ic).(CatField)(CatUM,CatColInd),...
                             'CooType',InPar.CooType,'ColNames',CatColNames,'ColUnits',CatUnits);

        % a matrix containing the requested columns
        CatEdgeCol     = InAstC.(CatField)(:,InPar.EdgeCol);

        % TBD:
        % The following line doesn't work - so the CatEdgeCol is appended
        % only to CatUM...
        %---- AstOut(Ic).(CatField)(CatUM,CurrentCol+1:CurrentCol+NedgeCol) = CatEdgeCol;
        AstOutCol(CurrentCol+1:CurrentCol+NedgeCol) = InAstC.(ColCellField)(InPar.EdgeCol);
        EdgeCol = InPar.EdgeCol;
        EdgeColNames = InAstC.(ColCellField)(EdgeCol);
    else
        EdgeCol    = false(1,4);
        CatEdgeCol = [];
        EdgeColNames = [];
        % bug fix
        AstOutCol = AstOutCol(1:end-4);
    end        

    % populate the AstOut catalog column names        
    AstOut(Ic).(ColCellField) = AstOutCol;   
    AstOut(Ic) = colcell2col(AstOut(Ic));
    
    
    if (nargout>1)
        % Unmatched sources in the catalog
        AstUM(Ic).(CatField)   = [AstC(Ic).(CatField)(CatUM,CatOutCol), CatEdgeCol];
        AstUM(Ic).ColCell      = [AstC(Ic).ColCell(CatOutCol), EdgeColNames];
        AstUM(Ic)              = colcell2col(AstUM(Ic));
        AstUM(Ic).ColUnits     = AstC(Ic).ColUnits;
        AstUM(Ic).Source       = sprintf('Unmatched sources in Cat (%s) compared with Ref',AstC(Ic).Source,AstR.Source);
        AstUM(Ic).SortedBy     = AstC(Ic).SortedBy;
        AstUM(Ic).SortedByCol  = AstC(Ic).SortedByCol;
    end
    
        
       
    
end


end


%---------------------
%--- Aux functions ---
%---------------------


function [ColNames,ColInd,UseInd]=select_colnames(ColNames,AstC)
    %------------
    % Output : - ColNames is a two element cell array of column names
    %          - ColInd is a two element array of column indices
    %          - UseInd is the index of the line in the input ColNames
    %            that was selected.

    if (~iscell(ColNames))
        % Assume the user supplied the column indices in a two element matrix
        ColInd = ColNames;
        ColNames  = ind2colname(AstC,ColInd);
        UseInd    = 1;
    else
        % The user supplied a cell array of column names

        % check how many pairs of column names the user supplied
        if (size(ColNames,2)~=2)
            error('Cell array of X/Y column names must contain two columns');
        end
        NxyCol = size(ColNames,1);
        if (NxyCol==1)
            % The user supplied only one pairs of keyword names - use as is.
            % Use ColNames as is
            UseInd = 1;
        else
            % The user supplied multiple pairs of column names
            % select the first pair that appears in the AstC catalog
            Found    = is_colname(AstC,ColNames);
            IndF     = find(Found,1,'first');
            UseInd   = IndF;   % line of RefColNames in which the columns names were found
            ColNames = ColNames(UseInd,:);
        end

        ColInd   = colname2ind(AstC,ColNames);
        if (numel(ColInd)~=2)
            error('Number of ColNames columns found is not 2');
        end



    end

end


function [Factor,CooType]=set_cootype(Units,CooType)
    % Set the value of CooType according to the units | default parameters
    % Example: [RefFactor,InPar.CooType] = set_cootype(RefUnits,InPar.CooType);
    % Set the coordinate types (sphere | plane)
    if (isempty(Units))
        % Use default CooType
        % assume units are pix or radians
        Factor = 1;
    else
        if (iscell(Units))
            Units = Units{1};
        end
        switch lower(Units)
            case {'rad','deg','radian'}
                % override CooType
                CooType = 'sphere';
            case {'pix','pixel'}
                % override CooType
                CooType = 'plane';
            otherwise
                error('Unknown RefUnits option');
        end

        % Set the conversion factor for the Ref catalog units        
        switch lower(CooType)
            case 'plane'
                % no conversion
                Factor = 1;
            case 'sphere'
                % spherical coordinates - conversion may be needed

                % Ref Conversion factor (to radians)
                %Factor = convert_units(Units,'rad');
                Factor = convert.angular(Units,'rad');
                % Cat Conversion factor (to radians)
                %CatFactor = convert_units(InPar.CatUnits,'rad');

            otherwise
                error('Unknown CooType option');
        end
    end

end


function OutCol=appended_columns_ind(OutCol,Ncol,AstR)
    % Example: RefOutCol=appended_columns_ind(InPar.RefOutCol,RefNcol,AstR)
    % Set the columns of Ref that will appended into the columns of Cat
    % Options are:
    %   NaN - no columns
    %   []  - all columns
    %   Cell array of column names, or vector of column indices
    if (isempty(OutCol))
        % Append all the columns of Ref into the output Cat
        OutCol = (1:1:Ncol);
    else  
        if (~iscell(OutCol))
            % OutCol is a vector
            if (isnan(OutCol))
                % Don't append Columns of Ref to the output
                OutCol = [];
            else
                % InPar.RefOutCol is already a vector of indices
                % OutCol = OutCol;
            end
        else
            % Convert InPar.RefOutCol to a vector of indices
            OutCol = colname2ind(AstR,OutCol);
        end
    end


end
