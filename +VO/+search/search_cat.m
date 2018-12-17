function [Res,CatUM,Cat]=search_cat(Cat,X,Y,varargin)
%--------------------------------------------------------------------------
% search_cat function                                            Catalogue
% Description: Given a catalog with Long,Lat coordinates position,
%              search for lines near a list of reference positions.
%              This function can be used to search for a near(est) position
%              in a catalog or to match two catalogs.
%              This function replaces cat_search.m and cat_match.m
% Input  : - Input catalog to saerch.
%            The input catalog should contains at least two columns
%            containing the longitude (X) and latitude (Y) coordinates.
%            The catalog can be a matrix,
%            a structure (or SIM) containing an astronomical catalog (e.g.,
%            FIRST.mat), or a string containing the name of a mat file
%            containing a structure of an an astronomical catalog (e.g.,
%            'FIRST.mat').
%            The last two options are available only if the number of
%            input arguments is not 2 (e.g., set Y=[] if needed).
%          - A single column vector containing the reference catalog
%            longitude to search,
%            or a two column matrix [Long, Lat] to search.
%            If two columns then the units are either radians or deg
%            (according to the value of 'IsRefdeg', default is radians).
%            If single column and IsRefdeg=false then will use convertdms
%            to convert the longitude into radians.
%            Note that if there are only two input argument the program
%            works in a fast mode, in which all the initilaization is
%            skipped.
%          - Either an empty matrix, or a single column vector of
%            latitudes. If IsRefdeg=false then will use convertdms
%            to convert the longitude into radians.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SearchRad' - Search radius (radians). Default is 2./(RAD*3600) (2").
%            'CooType' - {'sphere'|'plane'}. Default is 'sphere'.
%                        Alternatively this can be a function handle name
%                        which calculate the distance between two points
%                        (e.g., @sphere_dist).
%            'IsCatdeg'- A flag indicating if the input catalog coordinates
%                        are in radians (false) or deg (true).
%                        Default is false.
%            'IsRefdeg'- A flag indicating if the input reference catalog
%                        coordinates are in radians (false) or deg (true).
%                        Default is false.
%            'IsCatSorted'- A flag indicating if the input catalog is
%                        sorted {true|false}. Default is true.
%            'IsRad'   - If true, then will assume that input reference
%                        (search) coordinates are already in radians,
%                        and no conversion is needed. Default is false.
%            'ColX'    - Column index or name of the longitude in the input
%                        catalog. Default is 1.
%            'ColY'    - Column index or name of the latitude in the input
%                        catalog. Default is 2.
%            'SortCol' - By which column the input catalog is sorted.
%                        Either 'X'|'Y' or column index. Default is 'Y'
%                        (i.e., latitude).
%            'RefSearchCol' - Corresponding sorted column index in the
%                        reference catalog. Note that the reference catalog
%                        not need to be sorted.
%                        Default is 2.
%            'SearchMethod' - Options are:
%                        'find' - use find slow search.
%                        'bin'  - one-d binary search on sorted column.
%                        'binm' - simultaneous binary search. 
%                        'binms'- simultaneous binary search, with fast
%                                 implemntation of spherical ang. distance
%                                 and no position angle calculation.
%                                 This is valid only for CooType='sphere'.
%                                 In this case CooType can't be a funtion
%                                 handle. Default.
%                        'binms1' - like 'binms', but return the index
%                                 and distance of the neasrest source only.
%                        'binmdup' - like 'binms', but in this case the
%                                 search is only for non-identical
%                                 duplicate entries in the same list.
%                                 This option should be used when you like
%                                 to remove non-identical duplicate entries
%                                 from a list.
%                                 In this case CooType can't be a funtion
%                                 handle.
%            'FunBinS' - Binary search function handle.
%                        Default is @bin_sear2.
% Output : - Structure array with element per each entry in the reference
%            [X,Y] catalog. The structure containing the following
%            information:
%            .IndCat    - Vector of indices of matched sources in the
%                         sorted input catalog.
%            .DistRAD   - Vector of angular distances [radians] per each
%                         match.
%            .PA        - Vector of position angles [radians] per each
%                         match.
%            .Nfound    - Scalar indicating the number of matches.
%          - A boolean vector indicating if an entry in the (sorted)
%            input catalog was matched (false) or not (true).
%          - The sorted input catalog. The indices in the output structure
%            reffers to this catalog.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: R=sortrows(rand(1000000,2),2);
%          [Res,CatUM,Cat]=search_cat(R(1:10000,1:2),R(1:10000,1:2));
%          [Res,CatUM,Cat]=search_cat(R,'02:10:10.0','+38:45:12','SearchRad',1./RAD);
%          [Res,CatUM,Cat]=search_cat('FIRST.mat','12:10:10.0','+38:45:12','SearchRad',1./RAD);
%          [Res,CatUM,Cat]=search_cat(R,10.0,+12.0,'SearchRad',0.1./RAD,'IsRefdeg',true);
%          % match catalog
%          load FIRST.mat
%          [Res,CatUM,Cat]=search_cat(FIRST.Cat,FIRST.Cat(:,1:2));
% Reliable: 2
%--------------------------------------------------------------------------
import Util.find.*

InvRAD = pi./180;


DefV.SearchRad       = 2.*InvRAD./3600;
DefV.CooType         = 'sphere';  % {'sphere'|'plane'} or handle (e.g., '@sphere_dist')
DefV.IsCatdeg        = false;
DefV.IsRefdeg        = false;
DefV.IsCatSorted     = true;
DefV.IsRad           = false;
DefV.ColX            = 1;
DefV.ColY            = 2;
DefV.SortCol         = 'Y';     % either column index or 'X' | 'Y'
DefV.RefSearchCol    = 2;
DefV.SearchMethod    = 'binms';   % {'bin'|'find'|'binm'|'binms'|'binmdup'}
DefV.FunBinS         = @bin_sear; % @binary_search, @bin_sear, @find_bin
%DefV.CalcPA          = false;

if (nargin==2)
    % fast set up
    InPar         = DefV;
    InPar.SortCol = InPar.ColY;
    RefCat        = X;
else
    % slow (full) setup
    %InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);

    % SortCol - convert X/Y char to column index
    if (ischar(InPar.SortCol))
        switch lower(InPar.SortCol)
            case 'x'
                InPar.SortCol = InPar.ColX;
            case 'y'
                InPar.SortCol = InPar.ColY;
            otherwise
                error('Unknown SortCol option');
        end
    end

    % Deal with various Cat input class
    if (isnumeric(Cat))
        % do nothing
        % Cat is a matrix
    else
        if (ischar(Cat))
            % load catalog from a MAT file
            Cat = Util.IO.load2(Cat);
        end
        
        if (isstruct(Cat) && ~AstCat.isastcat(Cat))
            % assume Cat is a matlab astronomical catalog structure
            InPar.IsCatSorted = true;
            InPar.ColX        = Cat.Col.RA;
            InPar.ColY        = Cat.Col.Dec;
            Cat               = Cat.Cat;
            
        end
    end
    
    
    % deal with the possibility that Y is not provided
    if (nargin.*0.5==floor(nargin.*0.5))
        % Y is not provided
        varargin(2:end+1) = varargin;
        varargin{1}       = Y;
        Y                 = [];
    end


    % deal with X/Y input
    if (isempty(Y))
        % Y is empty
        Y = X(:,2);
        X = X(:,1);
        if (InPar.IsRefdeg && strcmpi(InPar.CooType,'sphere'))
            % convert to radians
            X = X.*InvRAD;
            Y = Y.*InvRAD;
        end
    else
        switch lower(InPar.CooType)
            case 'sphere'
                if (InPar.IsRefdeg)
                    % convert to radians
                    X = X.*InvRAD;
                    Y = Y.*InvRAD;
                else
                    if (~InPar.IsRad)
                        X = celestial.coo.convertdms(X,'gH','r');
                        Y = celestial.coo.convertdms(Y,'gD','R');
                    end
                end
            case 'plane'
                % do nothing
            otherwis
                error('Unknown CooType option');
        end
    end
    RefCat = [X,Y];
    
    % sort input Cat
    if (~InPar.IsCatSorted)
        Cat = sortrows(Cat,InPar.SortCol);
    end
    
    if (InPar.IsCatdeg)
        % Input catalog is in degrees
        % convert to radians
        Cat(:,[InPar.ColX, InPar.ColY]) = Cat(:,[InPar.ColX, InPar.ColY]).*InvRAD;
    end
end

if (isa(InPar.CooType,'function_handle'))
    DistFun = InPar.CooType;
else
    switch lower(InPar.CooType)
        case 'sphere'
            DistFun = @celestial.coo.sphere_dist_fast;
        case 'plane'
            DistFun = @plane_dist;
        otherwise
            error('Unknown CooType option');
    end
end

Ncat = size(Cat,1);  % number of sources in Catalog to search
Nin  = size(RefCat,1);    % number of search coordinates

CatUM = true(Ncat,1);
Res   = struct('IndCat',cell(1,Nin),'DistRAD',cell(1,Nin),'PA',cell(1,Nin),'Nfound',cell(1,Nin)); %num2cell(zeros(1,Nin)));
if (~isempty(Cat))
    % If Cat is empty, will return empty arrays
    
    switch lower(InPar.SearchMethod)
        case 'bin'
            for Is=1:1:Nin
                % assume Cat is sorted
                %Ireg              = (find_bin(Cat(:,InPar.SortCol),RefCat(Is,InPar.RefSearchCol)-InPar.SearchRad):1:find_bin(Cat(:,InPar.SortCol),RefCat(Is,InPar.RefSearchCol)+InPar.SearchRad));
                %Ireg              = (InPar.FunBinS(Cat(:,InPar.SortCol),RefCat(Is,InPar.RefSearchCol)-InPar.SearchRad):1:InPar.FunBinS(Cat(:,InPar.SortCol),RefCat(Is,InPar.RefSearchCol)+InPar.SearchRad));
                Ireg1             = InPar.FunBinS(Cat(:,InPar.SortCol),RefCat(Is,InPar.RefSearchCol)-InPar.SearchRad);
                Ireg2             = InPar.FunBinS(Cat(:,InPar.SortCol),RefCat(Is,InPar.RefSearchCol)+InPar.SearchRad);
                Ireg              = (Ireg1:1:Ireg2);
                [Dist,PA]         = DistFun(RefCat(Is,1),RefCat(Is,2),Cat(Ireg,InPar.ColX),Cat(Ireg,InPar.ColY));
                Ifr               = find(Dist<InPar.SearchRad);
                Res(Is).IndCat    = Ireg(Ifr);
                Res(Is).DistRAD   = Dist(Ifr);
                Res(Is).PA        = PA(Ifr);
                Res(Is).Nfound    = numel(Res(Is).IndCat);
                CatUM(Res(Is).IndCat) = false;

            end

        case 'binm'
            Ind1 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'-InPar.SearchRad);
            Ind2 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'+InPar.SearchRad)+1;

            for Is=1:1:Nin
%                if (Is==1866)
%                    save Eran.mat RefCat Cat Is Ireg InPar
%                end

                Ireg = (Ind1(Is):1:Ind2(Is));
                %if (InPar.CalcPA)
                    [Dist,PA]         = DistFun(RefCat(Is,1),RefCat(Is,2),Cat(Ireg,InPar.ColX),Cat(Ireg,InPar.ColY));
                    Ifr               = find(Dist<InPar.SearchRad);
                    Res(Is).PA        = PA(Ifr);
                %else
                %    [Dist]            = DistFun(RefCat(Is,1),RefCat(Is,2),Cat(Ireg,InPar.ColX),Cat(Ireg,InPar.ColY));
                %    Ifr               = find(Dist<InPar.SearchRad);
                %end

                Res(Is).IndCat    = Ireg(Ifr);
                Res(Is).DistRAD   = Dist(Ifr);

                %[Dist]         = DistFun(RefCat(Is,1),RefCat(Is,2),Cat(Ireg,InPar.ColX),Cat(Ireg,InPar.ColY));

                Res(Is).Nfound    = numel(Res(Is).IndCat);
                CatUM(Res(Is).IndCat) = false;

            end

       case 'binms'
            switch lower(InPar.CooType)
                case 'sphere'
                    % for spherical coordinates
                    Ind1 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'-InPar.SearchRad);
                    Ind2 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'+InPar.SearchRad);
                    Ind2 = min(Ind2+1,Ncat);
                    SinR1 = sin(Cat(:,InPar.ColX));
                    SinD1 = sin(Cat(:,InPar.ColY));
                    CosR1 = cos(Cat(:,InPar.ColX)); %sqrt(1-SinR1.^2);
                    CosD1 = cos(Cat(:,InPar.ColY)); %sqrt(1-SinD1.^2);

                    SinR2 = sin(RefCat(:,1));    % note this is RefCat!
                    SinD2 = sin(RefCat(:,2));
                    CosR2 = cos(RefCat(:,1)); % sqrt(1-SinR2.^2);
                    CosD2 = cos(RefCat(:,2)); % sqrt(1-SinD2.^2);
                    CosSearchRad = cos(InPar.SearchRad);

                    % for each source in RefCat

                    for Is=1:1:Nin
                        Ireg = (Ind1(Is):1:Ind2(Is));   % indices in Cat
                        %Ireg = (Is+1:1:Ind2(Is));

                        CosDist = SinD1(Ireg).*SinD2(Is) + CosD1(Ireg).*CosD2(Is).*(CosR1(Ireg).*CosR2(Is)+SinR1(Ireg).*SinR2(Is));
                        %Ifr               = find(CosDist>CosSearchRad);
                        Ifr               = (CosDist>CosSearchRad);

                        Res(Is).IndCat    = Ireg(Ifr);
                        Res(Is).DistRAD   = acos(CosDist(Ifr));

                        Res(Is).Nfound    = numel(Res(Is).IndCat);
                        CatUM(Res(Is).IndCat) = false;

                    end
                case 'plane'
                    % for plane coordinates
                    Ind1 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'-InPar.SearchRad);
                    Ind2 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'+InPar.SearchRad)+1;

                    % for each source in RefCat
                    SearchRad2 = InPar.SearchRad.^2;
                    for Is=1:1:Nin
                        Ireg = (Ind1(Is):1:Ind2(Is));   % indices in Cat
                        %Ireg = (Is+1:1:Ind2(Is));

                        Dist2 = (RefCat(Is,1) - Cat(Ireg,InPar.ColX)).^2 + (RefCat(Is,2) - Cat(Ireg,InPar.ColY)).^2;
                        %Ifr               = find(Dist2<SearchRad2);
                        Ifr               = (Dist2<SearchRad2);

                        Res(Is).IndCat    = Ireg(Ifr);
                        Res(Is).DistRAD   = sqrt(Dist2(Ifr));

                        Res(Is).Nfound    = numel(Res(Is).IndCat);
                        CatUM(Res(Is).IndCat) = false;

                    end
                otherwise
                    error('In case SearchMethod=binms CooType cannot be a function handle');
            end

        case 'binms1'
            switch lower(InPar.CooType)
                case 'sphere'
                    % for spherical coordinates
                    Ind1 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'-InPar.SearchRad);
                    Ind2 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'+InPar.SearchRad)+1;
                    SinR1 = sin(Cat(:,InPar.ColX));
                    SinD1 = sin(Cat(:,InPar.ColY));
                    CosR1 = cos(Cat(:,InPar.ColX)); %sqrt(1-SinR1.^2);
                    CosD1 = cos(Cat(:,InPar.ColY)); %sqrt(1-SinD1.^2);

                    SinR2 = sin(RefCat(:,1));    % note this is RefCat!
                    SinD2 = sin(RefCat(:,2));
                    CosR2 = cos(RefCat(:,1)); %sqrt(1-SinR2.^2);
                    CosD2 = cos(RefCat(:,2)); %sqrt(1-SinD2.^2);
                    CosSearchRad = cos(InPar.SearchRad);

                    % for each source in RefCat
                    for Is=1:1:Nin
                        Ireg = (Ind1(Is):1:Ind2(Is));   % indices in Cat
                        %Ireg = (Is+1:1:Ind2(Is));

                        CosDist = SinD1(Ireg).*SinD2(Is) + CosD1(Ireg).*CosD2(Is).*(CosR1(Ireg).*CosR2(Is)+SinR1(Ireg).*SinR2(Is));
                        Ifr               = find(CosDist>CosSearchRad);
                        Res(Is).Nfound    = length(Ifr);
                        [~,MinI]          = max(CosDist(Ifr));

                        Res(Is).IndCat    = Ireg(Ifr(MinI));
                        Res(Is).DistRAD   = acos(CosDist(Ifr(MinI)));

                        CatUM(Res(Is).IndCat) = false;  

                    end
                case 'plane'
                    % for plane coordinates
                    Ind1 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'-InPar.SearchRad);
                    Ind2 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'+InPar.SearchRad)+1;

                    % for each source in RefCat
                    SearchRad2 = InPar.SearchRad.^2;
                    for Is=1:1:Nin
                        Ireg = (Ind1(Is):1:Ind2(Is));   % indices in Cat
                        %Ireg = (Is+1:1:Ind2(Is));

                        Dist2 = (RefCat(Is,1) - Cat(Ireg,InPar.ColX)).^2 + (RefCat(Is,2) - Cat(Ireg,InPar.ColY)).^2;
                        Ifr               = find(Dist2<SearchRad2);
                        Res(Is).Nfound    = length(Ifr);
                        [~,MinI]          = min(Dist2(Ifr));

                        Res(Is).IndCat    = Ireg(Ifr(MinI));
                        Res(Is).DistRAD   = sqrt(Dist2(Ifr(MinI)));

                        CatUM(Res(Is).IndCat) = false; 

                    end
                otherwise
                    error('In case SearchMethod=binms CooType cannot be a function handle');
            end

         case 'binmdup'
            % search only for non-identical duplicate entries
             switch lower(InPar.CooType)
                case 'sphere'
                    % for spherical coordinates
                    %Ind1 = mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'-InPar.SearchRad);
                    Ind2 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'+InPar.SearchRad)+1;
                    SinR1 = sin(Cat(:,InPar.ColX));
                    SinD1 = sin(Cat(:,InPar.ColY));
                    CosR1 = cos(Cat(:,InPar.ColX)); %sqrt(1-SinR1.^2);
                    CosD1 = cos(Cat(:,InPar.ColY)); %sqrt(1-SinD1.^2);

                    SinR2 = sin(RefCat(:,1));
                    SinD2 = sin(RefCat(:,2));
                    CosR2 = cos(RefCat(:,1)); %sqrt(1-SinR2.^2);
                    CosD2 = cos(RefCat(:,2)); %sqrt(1-SinD2.^2);
                    CosSearchRad = cos(InPar.SearchRad);

                    IndReg = (1:1:Nin)';
                    for Is=1:1:Nin
                        %-1,
                        %Ireg = (Ind1(Is):1:Ind2(Is));
                        %Ireg = (Is+1:1:Ind2(Is));
                        Ireg = Util.array.delete_ind(IndReg,Is);

                        CosDist = SinD1(Ireg).*SinD2(Is) + CosD1(Ireg).*CosD2(Is).*(CosR1(Ireg).*CosR2(Is)+SinR1(Ireg).*SinR2(Is));
                        Ifr               = find(CosDist>CosSearchRad);

                        Res(Is).IndCat    = Ireg(Ifr);
                        Res(Is).DistRAD   = acos(CosDist(Ifr));

                        Res(Is).Nfound    = numel(Res(Is).IndCat);
                        CatUM(Res(Is).IndCat) = false;

                    end    
                case 'plane'
                    % for planar coordinates
                    %Ind1 = mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'-InPar.SearchRad);
                    Ind2 = Util.find.mfind_bin(Cat(:,InPar.SortCol),RefCat(:,InPar.RefSearchCol).'+InPar.SearchRad)+1;
                    
                    SearchRad2 = InPar.SearchRad.^2;
                 
                    IndReg = (1:1:Nin)';
                    for Is=1:1:Nin
                        %-1,
                        %Ireg = (Ind1(Is):1:Ind2(Is));
                        %Ireg = (Is+1:1:Ind2(Is));
                        Ireg = Util.array.delete_ind(IndReg,Is);
                        
                        Dist2 = (RefCat(Is,1) - Cat(Ireg,InPar.ColX)).^2 + (RefCat(Is,2) - Cat(Ireg,InPar.ColY)).^2;
                        Ifr               = find(Dist2<SearchRad2);
                        
                        [~,MinI]          = min(Dist2(Ifr));
                        Res(Is).IndCat    = Ireg(Ifr(MinI));
                        Res(Is).DistRAD   = sqrt(Dist2(Ifr(MinI)));

                        Res(Is).Nfound    = length(Ifr);
                        CatUM(Res(Is).IndCat) = false;
                    end
                    
                otherwise
                    error('In case SearchMethod=binmdup CooType cannot be a function handle');
            end
        case 'find'
            for Is=1:1:Nin
                [Dist,PA]         = DistFun(RefCat(Is,1),RefCat(Is,2),Cat(:,InPar.ColX),Cat(:,InPar.ColY));
                Res(Is).IndCat    = find(Dist<InPar.SearchRad);
                Res(Is).DistRAD   = Dist(Res(Is).IndCat);
                Res(Is).PA        = PA(Res(Is).IndCat);
                Res(Is).Nfound    = numel(Res(Is).IndCat);
                CatUM(Res(Is).IndCat) = false;
            end

        case 'kd'
            error('kd tree search not supported yet');

        otherwise
            error('Unknown SearchMethod option');
    end
end











