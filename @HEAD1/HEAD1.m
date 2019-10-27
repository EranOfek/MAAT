%--------------------------------------------------------------------------
% HEAD class                                                         class
% Description: A class of for an image header (Header).
%              This class can be use to store header meta data information
%              associated with images or spectra, and to access and search
%              this information.
%              An header information consist of multiple keyword names.
%              Each keyword have a value and an optional comment.
%              The Header information can be stored in a 
%              three columns cell array {Keyword,Value,Comment}.
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef HEAD1 %< WorldCooSys
    properties (SetAccess = public)
        Header
        Keys
        WCS
    end
    
    % consider add hidden properties:
    % isaligned
    
%     properties (Hidden = true)
%         %mean1
%         header
%     end

    % constructor
    methods

        %-------------------------
        %--- Class constructor ---
        %-------------------------
        
        function H=HEAD1(varargin)
            % HEAD constructor
            % Package: @HEAD
            % Description: Constrct an HEAD object with specific number of
            %              elements, and with an empty Header field.
            % Input  : * Arbitrary number of number of elements in each
            %            dimension. Deafult is 1,1
            % Output : - An HEAD object
            % Example: HEAD(2,3,2)
            %          HEAD([2,3,2])
            %          HEAD(2)
            %          HEAD
            
            if (numel(varargin)==0)
                varargin{1} = 1;
                varargin{2} = 1;
            elseif (numel(varargin)==1)
                
                if numel(varargin{1})>1
                    varargin = num2cell(varargin{1});
                else
                    varargin{2} = 1;
                end
            else
                % do nothing
            end
            Dim = cell2mat(varargin);
            
            Nel = prod(Dim);
            
            for I=1:1:Nel
                H(I).Header = cell(0,3);
            end
            
            H=reshape(H,varargin{:});
            
        end
        
        
    end
    
    % getters/setters
    methods
        function H=get.Header(Head)
            % getter for cell Header from an HEAD object
            % Package: @HEAD
            % Input  : - An HEAD object
            % Output : - If the HEAD object contains a single element, than
            %            this is the content of the Header field (i.e., a 3
            %            column cell array).
            %            Otherwise this is a cell array of cell arrays.
            % By : Eran O. Ofek
            
            
            N = numel(Head);
            % first update Header according to Key structure
            % need to limit infinite recusion / BUG
            % Head = updateHeader_fromStruct(Head);
            
            if (N==1)
                H = Head.Header;
                
            else
                for I=1:1:N
                    H{I} = Head(I).Header;
                    
                end
            end
            
        end
        
        function S=get.Keys(H)
            % getter for struct Header from an HEAD object
            % Package: @HEAD
            % Description: Populate the Keys field in an HEAD object and
            %              return the structure array of header keywords
            %              and their values.
            % Input  : - An HEAD object
            % Output : - Structure array of keywords and their values.
            % By : Eran O. Ofek
            % Example: H=HEAD.basic; H(1).Keys, H(1).Keys.TYPE
            % Reliable: 2
            
            
            N = numel(H);
            for I=1:1:N
                if isempty(H(I).Keys)
                    H(I).Keys = cell2struct(H(I));
                end
                S(I) = H(I).Keys;
            end
            
        end
        
%         function H=set.Keys(H,Val)
%             %
%             Val
%             if isstruct(Val)
%                 % assign new values into the cell header
%                 Val
%                 
%                 H.Header = HEAD1.stvruct2cell(H.Keys);
%             else
%                 warning('You asigned a non structure object into Keys');
%             end
%             
%         end
        
    end
    
    
    
    % Static methods: isHEAD, HeaderField, basic
    methods (Static)
        
        function Ans=isHEAD(Obj)
            % Return true if object is HEAD
            % Description: Check if object is of HEAD class.
            % Input  : - Object
            % Output : - {true|false}.
            %     By : Eran O. Ofek                    Oct 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: S=1; HEAD.ishead(S);
            % Reliable: 2
            Ans = isa(Obj,'HEAD1');
        end
        
        function Name=HeaderField
            % Description: Return the header field name in HEAD
            % Input  : null
            % Output : Header field name
            Name = 'Header';
        end
        
        function H=basic
            % Construct a basic HEAD object populated with some fields
            % Package: @HEAD
            % Description: Construct a basic HEAD object populated with
            %              some fields. This is used mainly for testing.
            % Input  : null
            % Output : - An HEAD object with some Header fields.
            % Example: H=HEAD1.basic
            % Reliable: 2
           
            H = HEAD1;
            H.Header = {'SIMPLE',true','does file conform to the Standard?';
                        'BSCALE',1.0,'linear factor in scaling equation';
                        'BZERO',0.0,'zero point in scaling equation';
                        'BUNIT','DN','physical units of the array values';
                        'BITPIX',32,'';
                        'EXPTIME',1,'';
                        'TYPE','bias','';
                        'FILTER','r','';
                        'NAXIS',2,'';
                        'NAXIS1',1024,'';
                        'NAXIS2',2048,'';
                        'UTC-OBS','2000-01-01T00:00:00.0','date';
                        'END','',''};
        end
       
        function [Dict,Val]=dictSearch(Dict,Key,Exact)
            % Search string in a dictionary of strings
            % Package: @HEAD
            % Description: Search string in a dictionary of strings. The
            %              dictionary is a structure array in which the
            %              "Key" field indicate a string and the "Synonym"
            %              field contains a cell array of synonyms. The
            %              synonyms may be a regular expressions.
            % Input  : - Dictionary structure array.
            %            E.g., D(1).Key='TYPE';
            %                  D(2).Synonym={'TYPE','\w*[tT][yY][pP][eE]'};
            %          - Keyword name to search.
            %            The keyword is first searched in the .Key fields,
            %            and if not found it is searched in the .Synonym
            %            fields.
            %          - A logical indicating if to use exact search
            %            (true), or regular exprssions (false.
            %            Default is true.
            % Output : - A structure array of dictionary items that were 
            %            found.
            %          - The the Dictionary '.Key' string that was found.
            %            If more then one were found than this is the
            %            first. Return empty if nothing was found.
            % Example: D=HEAD1.dictKey; Dv=HEAD1.dictValType;
            %          HEAD1.dictSearch(D,'TYPE')
            %          HEAD1.dictSearch(D,'IMTYPE')
            %          HEAD1.dictSearch(Dv,'superbias',false)
            % Reliable: 2
            
            if nargin<3
                Exact = true;
            end
            
            if Exact
                Ikey = find(strcmp(Key,{Dict.Key}));
            else
                Tmp = regexp({Dict.Key},Key,'match');
                Ikey = find(~Util.cell.isempty_cell(Tmp));
            end

            if isempty(Ikey)
                % Search Key in Dict.Synonym
                Ndict = numel(Dict);
                Ikey = [];
                for Idict=1:1:Ndict
                    if Exact
                        if any(strcmp(Key,Dict(Idict).Synonym))
                            Ikey = [Ikey; Idict];
                        end
                    else
                        %Tmp = regexp(Dict(Idict).Synonym,Key,'match');
                        Tmp = regexp(Key,Dict(Idict).Synonym,'match');
                        if any(~Util.cell.isempty_cell(Tmp))
                            Ikey = [Ikey; Idict];
                        end
                    end
                end
                if isempty(Ikey)
                    Dict = [];
                else
                    Dict = Dict(Ikey);
                end
            else
                Dict = Dict(Ikey);
            end
            if isempty(Dict)
                Val = [];
            else
                Val = Dict(1).Key;
            end
            
            
        end % end dictSearch function
        
        function [Dict,Val]=dictKey(Key,Exact)
            % Keywords synomyms dictionaty
            % Package: @HEAD
            % Description: Dictionary of synonyms for keyword names.
            %              Each keyword may have a list of synomyms.
            %              In case the keyword doesn't exist, the synonym
            %              may replace the keyword. The order of synonyms
            %              have meaning, where the first one is more
            %              likely.
            % Input  : - Optional Keyweord name to search in dictionary.
            %            Default is empty. If empty, then return entire
            %            dictionary.
            %          - Flag indicating if to perform exact search (true)
            %            or regexp search (false).
            %            Default is true.
            % Output : - Structure array of 'Key' and 'Synonym'.
            %          - The the Dictionary '.Key' string that was found.
            %            If more then one were found than this is the
            %            first. Return empty if nothing was found.
            % Example: Dict = HEAD1.dictKey;
            %          Dict = HEAD1.dictKey('EXPTIM',false);
            %          Dict = HEAD1.dictKey('AEXPTI',false);
            %          Dict = HEAD1.dictKey('EXPTIME')
            % Reliable: 2
            
            if nargin<2
                Exact = true;
                if nargin<1
                    Key = [];
                end
            end
            
            I = 0;
            I = I + 1;
            Dict(I).Key     = 'EXPTIME';
            Dict(I).Synonym = {'AEXPTIME','EXPTIME','EXPOSURE'};
            I = I + 1;
            Dict(I).Key     = 'RA';
            Dict(I).Synonym = {'RA','RADEG','RIGHTASC','JRA','RASEX','TELRA','OBJRA'};
            I = I + 1;
            Dict(I).Key     = 'DEC';
            Dict(I).Synonym = {'DEC','DECDEG','DECLINAT','JDEC','TELDEC','OBJDEC'};
            I = I + 1;
            Dict(I).Key     = 'HA';
            Dict(I).Synonym = {'HA','HOURANG'};
            I = I + 1;
            Dict(I).Key     = 'TYPE';
            Dict(I).Synonym = {'TYPE','OBSTYPE','IMTYPE','IMGTYPE','IMGTYP'};
            I = I + 1;
            Dict(I).Key     = 'OBSLON';
            Dict(I).Synonym = {'OBSLON','OBSLONG','LONG','LON'};
            I = I + 1;
            Dict(I).Key     = 'OBSLAT';
            Dict(I).Synonym = {'OBSLAT','OBSLAT','LAT'};
            I = I + 1;
            Dict(I).Key     = 'HEIGHT';
            Dict(I).Synonym = {'HEIGHT','OBSHEIGH'};
            I = I + 1;
            Dict(I).Key     = 'EQUINOX';
            Dict(I).Synonym = {'EQUINOX','EQUIN','EPOCH'};
            I = I + 1;
            Dict(I).Key     = 'READNOI';
            Dict(I).Synonym = {'READNOI','READNOIS','READNO','RN'};
            I = I + 1;
            Dict(I).Key     = 'FILTER';
            Dict(I).Synonym = {'FILTER','FILTERID','FILTERSL','FILTER1','FILT'};
            I = I + 1;
            Dict(I).Key     = 'FILTER2';
            Dict(I).Synonym = {'FILTER2','FILTERID2','FILTERSL2','FILTER1','FILT2'};
            I = I + 1;
            Dict(I).Key     = 'JD';
            Dict(I).Synonym = {'JD','OBSJD'};
            I = I + 1;
            Dict(I).Key     = 'MIDJD';
            Dict(I).Synonym = {'MIDJD','JDMID'};
            I = I + 1;
            Dict(I).Key     = 'MJD';
            Dict(I).Synonym = {'MJD','OBSMJD'};
            I = I + 1;
            Dict(I).Key     = 'MIDMJD';
            Dict(I).Synonym = {'MIDMJD','MJDMID'};
            I = I + 1;
            Dict(I).Key     = 'UTC-OBS';
            Dict(I).Synonym = {'UTC-OBS','UTC_OBS','UTCOBS',...
                               'UTC',...
                               'DATE-OBS','DATE_OBS','DATEOBS',...
                               'OBS-DATE','OBS_DATE','OBSDATE',...
                               'TIME-OBS','TIME_OBS','TIMEOBS',...
                               'DATE','ISODATE','DATEISO'};
            
            
            % Searck Key in Dict.Key
            if ~isempty(Key)
                [Dict,Val]=HEAD1.dictSearch(Dict,Key,Exact);
            else
                if nargout>1
                    error('Can not return searched value when keyword search is not provided');
                end
            end
            

        end % end dictKey function
        
        function [Dict,Val]=dictValType(Key,Exact)
            % dictionary of TYPE keyword values (e.g., bias, flat,...)
            % Package: @HEAD
            % Description: Dictionary of TYPE keyword values. The
            %              dictionary is a structure array in which the
            %              "Key" field indicate a string and the "Synonym"
            %              field contains a cell array of synonyms. The
            %              synonyms may be a regular expressions.
            % Input  : - Optional Keyweord name to search in dictionary.
            %            Default is empty. If empty, then return entire
            %            dictionary.
            %          - Flag indicating if to perform exact search (true)
            %            or regexp search (false).
            %            Default is false.
            % Output : - A dictionary strcture array there were found.
            %            The
            %            dictionary is a structure array in which the
            %            "Key" field indicate a string and the "Synonym"
            %            field contains a cell array of synonyms. The
            %            synonyms may be a regular expressions.
            %          - The the Dictionary '.Key' string that was found.
            %            If more then one were found than this is the
            %            first. Return empty if nothing was found.
            % Example: Dict=HEAD1.dictValType
            %          Dict=HEAD1.dictValType('superflat',false)
            %          Dict=HEAD1.dictValType('superflat',true); % return empty
            %          Dict=HEAD1.dictValType('Flat');
            % Reliable: 2
            
            if nargin<2
                Exact = false;
                if nargin<1
                    Key = [];
                end
            end
            
            I = 0;
            I = I + 1;
            Dict(I).Key     = 'bias';
            Dict(I).Synonym = {'\w*[bB][iI][aA][sS]\w*'};
            I = I + 1;
            Dict(I).Key     = 'flat';
            Dict(I).Synonym = {'\w*[fF][lL][aA][tT]\w*'};
            I = I + 1;
            Dict(I).Key     = 'dark';
            Dict(I).Synonym = {'\w*[dD][aA][rR][kK]\w*'};
            I = I + 1;
            Dict(I).Key     = 'science';
            Dict(I).Synonym = {'\w*[sC][cC][iI][eE][nN][cC][eE]\w*','\w*[oO][bB][jJ][eE][cC][tT]\w*'};
            I = I + 1;
            Dict(I).Key     = 'weight';
            Dict(I).Synonym = {'\w*[wW][eE][iI][gG][hH][tT]\w*'};
            I = I + 1;
            Dict(I).Key     = 'fringe';
            Dict(I).Synonym = {'\w*[fF][rR][iI][nN][gG][eE]\w*'};
            I = I + 1;
            Dict(I).Key     = 'coadd';
            Dict(I).Synonym = {'\w*[cC][oO][aA][dD][dD]\w*'};
            
            % Searck Key in Dict.Key
            if ~isempty(Key)
                [Dict,Val]=HEAD1.dictSearch(Dict,Key,Exact);
             else
                if nargout>1
                    error('Can not return searched value when keyword search is not provided');
                end
            end
            
            
        end
    end
    
    % static methods for cell header
    methods (Static)
        
        % get key value from a header in a cell format
        function [SubCell,FlagI]=getKeyCell(Cell,Key,Exact,Col,ReturnN)
            % get keyword value from an header in a cell format
            % Package: @HEAD
            % Description: Given a 3 column cell array [Key, Val, Comment]
            %              search for a keyword name and return its value
            %              and comment.
            % Input  : - A 3 columnn cell array [Key, Val, Comment].
            %          - Keyword name, or a cell array of keyword names.
            %          - A flag indicating if to perform an exact name
            %            search (true). Default is true.
            %            If false, then use regexp.
            %          - Column in the cell array in which to perform the
            %            search. Default is 1 (i.e., the keyword column).
            %            Note that if 2 is used, than this may fail if the
            %            some values are numeric.
            %          - If more than one keyword was matched, this is the
            %            number of first keywords/values to return.
            %            If Inf, then retun all lines containing the
            %            keyword name.
            %            Default is Inf.
            % Output : - A 3 column cell array with only the lines that
            %            satisfy the search.
            %          - Line indices in the input cell that satisfy the
            %            search criteria.
            % By: Eran O. Ofek
            % Example: H=HEAD1.basic; [SubCell]=HEAD1.getKeyCell(H.Header,{'TYPE','B'})
            %          H=HEAD1.basic; [SubCell]=HEAD1.getKeyCell(H.Header,{'TYP','EXP'},false)
            %          H=HEAD1.basic;
            %          [SubCell]=HEAD1.getKeyCell([H.Header;H.Header],{'TYP','EXP'},false,1,1)
            % Reliable: 2
            
            if (nargin<5)
                ReturnN = Inf;
                if (nargin<4)
                    Col = 1;
                    if (nargin<3)
                        Exact     = true;
                    end
                end
            end
            
            if (ischar(Key))
                Key = {Key};
            end
            Nkey = numel(Key);
            % for each keyword name
            FlagI = [];
            for Ikey=1:1:Nkey
                if (Exact)
                    % use exact name search
                    Flag = strcmp(Cell(:,Col),Key{Ikey});
                else
                    % use regular expression
                    Flag = regexp(Cell(:,Col),Key{Ikey},'match');
                    Flag = ~Util.cell.isempty_cell(Flag);
                end
                %Nfound{Ikey} = 
                FlagI = [FlagI; find(Flag,ReturnN,'first')];
                
            end
            SubCell = Cell(FlagI,:);
                    
        end  % end getKeyCell function
        
        function [SubCell]=getFilledKeyCell(Cell,Key,FillVal,ReturnI,Col)
            % get key from an header in a cell format, NaN if not exist
            % Package: @HEAD
            % Description: Given a 3 column cell array [Key, Val, Comment]
            %              search for exact keyword name and return one
            %              line for each searched keyword. Value is NaN is
            %              keyword doesn't exist.
            %              If more than one keyword exist then select only
            %              one.
            %              For non exact matching without fill values use:
            %              HEAD.getKeyCell.m
            % Input  : - A 3 columnn cell array [Key, Val, Comment].
            %          - Keyword name, or a cell array of keyword names.
            %          - Value to populate the keyword value in case the
            %            keyword doesn't exist. Default is NaN.
            %          - If more than one keyword was matched, this is the
            %            index of the keywords/values to return.
            %            (in anycase only one line will be returned).
            %            Inf will return the last keyword.
            %            Default is 1.
            %          - Column in the cell array in which to perform the
            %            search. Default is 1 (i.e., the keyword column).
            %            Note that if 2 is used, than this may fail if the
            %            some values are numeric.
            % Output : - A 3 column cell array with only the lines that
            %            satisfy the search.
            % By: Eran O. Ofek
            % Example: H=HEAD1.basic;
            %          [SubCell]=HEAD1.getFilledKeyCell(H.Header,{'TYPE','B'},NaN)
            %          [SubCell]=HEAD1.getFilledKeyCell([H.Header;H.Header],{'TYPE','B'},NaN,2)
            % Reliable: 2
            
            Def.ReturnI = 1;
            if (nargin<5)
                Col = 1;
                if (nargin<4)
                    ReturnI = Def.ReturnI;
                    if (nargin<3)
                        FillVal = NaN;
                    end
                end
            end
            
            if (ischar(Key))
                Key = {Key};
            end
            Nkey = numel(Key);
            % for each keyword name
            SubCell = cell(Nkey,3);
            SubCell(:,1) = Key(:);
            for Ikey=1:1:Nkey
                
                % use exact name search
                Flag  = strcmp(Cell(:,Col),Key{Ikey});
                FlagI = find(Flag);
                if (numel(FlagI)<ReturnI)
                    ReturnI = numel(FlagI);
                else
                    ReturnI = Def.ReturnI;
                end
                if isempty(FlagI)
                    SubCell{Ikey,2} = FillVal;
                    SubCell{Ikey,3} = '';
                else
                    FlagI = FlagI(ReturnI);
                    SubCell(Ikey,2:3) = Cell(FlagI,2:3);
                end
            end
                   
        end
        
        function [IsExist,Counter]=isKeyExistCell(Cell,Key,Exact,Col)
            % Check (and count) if a keyword name exist in a cell header
            % Package: @HEAD
            % Description: Check if keyword name exist in a cell header,
            %              and also count the number of appearnces for each
            %              keyword. Can either use exact match or regular
            %              exprssions.
            % Input  : - A cell header {Keyword, Value, Comment}
            %          - A keyword string or a cell array of keyword
            %            strings to search.
            %          - A flag indicating if to perform exact match (true)
            %            or regular exprssion matching (false).
            %          - Column index in the cell header on which to
            %            perfoprm the search. Default is 1.
            % Output : - A row vector of logical indicating if each keyword
            %            exist in the header.
            %          - A row vector of counters of the number of
            %            appearnces of each keyword in the header.
            % Example: H=HEAD1.basic;
            %          [IsExist,Counter]=HEAD1.isKeyExistCell([H.Header;H.Header],{'TYPE','N','EXPTI'},false)
            
            if nargin<4
                Col = 1;
                if nargin<3
                    Exact = true;
                end
            end
            
            
            if ischar(Key)
                Key = {Key};
            end
            Nkey = numel(Key);
            IsExist = false(1,Nkey);
            Counter = zeros(1,Nkey);
            for Ikey=1:1:Nkey
                if Exact
                    Flag = strcmp(Cell(:,Col),Key{Ikey});
                else
                    Flag = regexp(Cell(:,Col),Key{Ikey},'match');
                    Flag = ~Util.cell.isempty_cell(Flag);
                end
                IsExist(Ikey) = any(Flag);
                Counter(Ikey) = sum(Flag);
            end
                    
            
        end
        
        % delete key from cell header
        function [SubCell]=delKeyCell(Cell,Key,Exact,Col)
            % delete specific keyword name from an header in a cell format
            % Package: @HEAD
            % Description: Given a 3 column cell array [Key, Val, Comment]
            %              search for a keyword name and return a cell
            %              array without the lines containing this keyword.
            % Input  : - A 3 columnn cell array [Key, Val, Comment].
            %          - Keyword name, or a cell array of keyword names.
            %          - A flag indicating if to perform an exact name
            %            search (true). Default is true.
            %            If false, then use regexp.
            %          - Column in the cell array in which to perform the
            %            search. Default is 1 (i.e., the keyword column).
            %            Note that if 2 is used, than this may fail if the
            %            some values are numeric.
            % Output : - A 3 column cell array with only the lines that
            %            satisfy the search.
            % By: Eran O. Ofek
            % Reliable: 2
            
            
            if (nargin<4)
                Col = 1;
                if (nargin<3)
                    Exact     = true;
                end
            end
            
            if (ischar(Key))
                Key = {Key};
            end
            Nkey = numel(Key);
            % for each keyword name
            FlagA = false(size(Cell,1),1);
            for Ikey=1:1:Nkey
                if (Exact)
                    % use exact name search
                    Flag = strcmp(Cell(:,Col),Key{Ikey});
                else
                    % use regular expression
                    Flag = regexp(Cell(:,Col),Key{Ikey},'match');
                    Flag = ~Util.cell.isempty_cell(Flag);
                end
                
                FlagA = FlagA | Flag;
                
            end
              
            SubCell = Cell(~FlagA,:);
            
        end
        
        % remove duplicates from cell header
        function [SubCell,Ia,Ic]=uniqueKeyCell(Cell,Occur)
            % delete duplicate keywords from a Cell array header
            % Package: @HEAD
            % Description: Given a 3 column cell array [Key, Val, Comment]
            %              search for duplicate keywords and leave only the
            %              first or last keyword.
            % Input  : - A 3 columnn cell array [Key, Val, Comment].
            %          - Which occurance of duplicate keywords to save.
            %            'first' | 'last. Default is 'first.
            % Output : - A 3 column cell array with unique keywords.
            %          - Index Ia, for SubCell = Cell(Ia), where SubCell is
            %            the output and Cell is the input.
            %          - Index Ic, for Cell = SubCell(Ic).
            % By: Eran O. Ofek
            % Reliable: 2
            
            
            if (nargin<2)
                Occur = 'first';
            end
            
            [~,Ia,Ic] = unique(Cell(:,1),Occur);
            SubCell = Cell(Ia,:);
        end
        
        % concat cell headers
        function [Cell]=addKeyCell(Cell1,Cell2,Location)
            % Add a cell header to another cell array header
            % Package: @HEAD
            % Description: Given a 2/3 column cell array [Key, Val, Comment]
            %              add to it another 3 column or 2 coluymn cell
            %              array in a given location.
            % Input  : - A 2/3 columnn cell array [Key, Val, Comment].
            %          - A 2/3 columnn cell array [Key, Val, Comment].
            %          - Location in which to add the 2nd cell to the 1st
            %            cell. Inf for adding at the end (but before the
            %            END keyword), 0 for adding
            %            at the begining. Default is Inf.
            % Output : - A 3 column cell array of the combined two arrays.
            % By: Eran O. Ofek
            % Example: HEAD.addKeyCell({'EXPTIME',2;'TYPE','bias'},{'FILTER','r',''})
            % Reliable: 2
            
            if (nargin<3)
                Location = Inf;
            end
            
            [Nrow1,Ncol1] = size(Cell1);
            [Nrow2,Ncol2] = size(Cell2);
            
            if Ncol1==3
                % do nothing
            elseif Ncol1==2
                [Tmp{1:Nrow1}]=deal('');
                Cell1 = [Cell1, Tmp.'];
            else
                error('Input cell array must have 2 or 3 columns');
            end
                
            if Ncol2==3
                % do nothing
            elseif Ncol2==2
                [Tmp{1:Nrow2}]=deal('');
                Cell2 = [Cell2, Tmp.'];
            else
                error('Input cell array must have 2 or 3 columns');
            end
            
            if Location==Inf || Location>Nrow1
                % concat Cell2 at end of Cell1
                if strcmp(Cell1{end,1},'END')
                    Cell = [Cell1(1:end-1,:); Cell2; Cell1(end,:)];
                else
                    Cell = [Cell1;Cell2];
                end
            elseif Location==0
                Cell = [Cell2;Cell1];
            else
                Cell = [Cell1(1:Location,:); Cell2; Cell1(Location+1:end,:)];
            end
            
        end
        
        % replace keyword name in cell header
        function Cell=replaceKeyNameCell(Cell,OldKeyName,NewKeyName)
            % Replace keyword name in cell header
            % Package: @HEAD
            % Description: Replace keyword name in cell header.
            %              If keyword name is not found than skip the
            %              replecment.
            % Input  : - A 3 column cell array {keyname, value, comment}
            %          - A string or a cell array of strings containing
            %            the old keyword names to search in the header.
            %          - Astring or a cell array of strings containing
            %            new keyword names that will replace the old names.
            %            This input should have the same 
            %            number of keywords as the second input arguments.
            % Output : - A cell array in which the keyword names were
            %            replaced.
            
            Col = 1;
            if ischar(OldKeyName)
                OldKeyName = {OldKeyName};
            end
             if ischar(NewKeyName)
                NewKeyName = {NewKeyName};
             end
            if numel(OldKeyName)~=numel(NewKeyName)
                error('Number of new keywords names should be equal to old keyword names');
            end
            
            Nkey = numel(NewKeyName);
            for Ikey=1:1:Nkey
                Ipos = find(strcmp(Cell(:,Col),OldKeyName{Ikey}));
                if (~isempty(Ipos))
                    Cell{Ipos,Col} = NewKeyName{Ikey};
                end
            end
            
            
        end
        
        % replace keyword value in cell header
        function Cell=replaceKeyValCell(Cell,KeyName,KeyVal)
            % Replace keyword value in cell header
            % Package: @HEAD
            % Description: Replace keyword value in cell header.
            % Input  : - A 3 column cell array {keyname, value, comment}
            %          - A string or a cell array of strings containing
            %            the keyword names to search in the header,
            %            and for which to replace the value.
            %          - Keyword values. The same number of values as
            %            keyword names.
            %            This can be a scalar, a cell array or a vector.
            % Output : - A cell array in which the keyword values were
            %            replaced.
            
            Col = 1;
            ColVal = 2;
            if ischar(KeyName)
                KeyName = {KeyName};
            end
            if ischar(KeyVal)
                KeyVal = {KeyVal};
            elseif isnumeric(KeyVal)
                KeyVal = num2cell(KeyVal);
            else
                % do nothing
            end
            if numel(KeyName)~=numel(KeyVal)
                error('Number of keywords names should be equal to number of values');
            end
            
            Nkey = numel(KeyName);
            for Ikey=1:1:Nkey
                Ipos = find(strcmp(Cell(:,Col),KeyName{Ikey}));
                if (~isempty(Ipos))
                    Cell{Ipos,ColVal} = KeyVal{Ikey};
                end
            end
            
            
        end
            
        function [Cell]=fixColCell(Cell)
            % If a cell array contains two columns, add a third with ''
            % Package: @HEAD
            % Description: If a cell array contains two columns, add a
            %              third column with empty comments ('').
            % Input  : - A cell header with 2 or 3 columns.
            % Output : - A cell header with 3 columns.
            % Reliable: 2
            
            [Nrow,Ncol] = size(Cell);
            
            if Ncol==3
                % do nothing
            elseif Ncol==2
                [Tmp{1:Nrow}]=deal('');
                Cell = [Cell, Tmp.'];
            else
                error('Input cell array must have 2 or 3 columns');
            end
            
                
        end
        
        function Cell=struct2cell(S)
            % Convert a structure into a 3 column cell header
            % Package: @HEAD
            % Description: Convert a single-element structure into a 3
            %              column cell header.
            %              For conversion of a structure array, use the
            %              HEAD.struct2head static method.
            % Input  : - A structure.
            % Output : - A 3 column cell array {Keyword, Value, Comment},
            %            where the Keywords are taken from the structure
            %            field names, and the values are taken from the
            %            field values.
            % Example: HEAD1.struct2cell(S)
            % Reliable: 2
        
            if numel(S)>1
                error('HEAD.struct2cell works on a single element structure - use HEAD.struct2head');
            end
            FN = fieldnames(S);
            Cell = HEAD1.fixColCell([FN(:), struct2cell(S)]);
            
        end
        
        function H=struct2head(S)
            % Convert a structure array into an HEAD object
            % Package: @HEAD
            % Description: Convert a structure array into an HEAD object
            % Input  : - A structure array.
            % Output : - An HEAD object. Each element in the structure
            %            array corresponds to an element in the HEAD object.
            %            Each field/value in the structure corresponds to
            %            a line in the cell header.
            % Example: H2=HEAD1.struct2head([S;S])
            % Reliable: 2
            
            H = HEAD1(size(S));
            N = numel(S);
            for I=1:1:N
                H(I).Header = HEAD1.struct2cell(S(I));
            end
            
            
        end
        
        function C=conv2numCell(C,KeepAsStr)
            % Attempt to convert all elements in a cell array to numeric
            % Package: @HEAD
            % Description: Attempt to convert all elements in a cell array
            %              to numeric.
            % Input  : - A cell array.
            %          - A logical flag indicating if to keep strings as
            %            strings (true), or to replace them with NaNs
            %            (false). Default is true.
            % Output : - A cell array, with numeric values.
            % Example: HEAD1.conv2numCell({'1','NaN',NaN,'-Inf',11,'aa'})
            
            if nargin<2
                KeepAsStr = true;
            end
           
            N = numel(C);
            for I=1:1:N
                if ischar(C{I})
                    if strcmpi(C{I},'nan')
                        C{I} = NaN;
                    else
                        Tmp = str2double(C{I});
                        if isnan(Tmp)
                            if ~KeepAsStr
                                C = NaN;
                            end
                        else
                            C{I} = Tmp;
                        end
                    end
                end
            end
            
        end
        
        function C=replace_iligal_charCell(C,Col)
            % replace iligal characters in keyword names (e.g. '-').
            % Package: @HEAD
            % Description: Given a 3 column cell array. Replace iligal
            %              characters in the first column.
            %              List of illigal characters:
            %              '-' replace with '_'
            % Input  : - A 3 columns cell array
            %          - Column index in which to replace characaters.
            %            Default is 1.
            % Output : - A 3 column cell array.
            % Example: H=HEAD1.basic; H=H.addKey({'UTC-OBS','2009'});
            %          HEAD1.replace_iligal_charCell(H.Header)
            % Reliable: 2
            
            if nargin<2
                Col = 1;
            end
            
            C(:,Col) = regexprep(C(:,Col),'-','_');
            
        end
        
        function C=updateCell_fromStruct(S,C)
            % update cell header by adding elements in struct not in header
            % Package: @HEAD
            % Description: Given a 3 column cell array (i.e., cell header)
            %              and a structure with keywords and values.
            %              First make sure all the header keywords have no
            %              iligal characters. Then look for fields in the
            %              structure that are not in the cell header and
            %              add them to the cell header.
            % Input  : - A structure.
            %          - A cell header.
            % Output : - An updated cell header.
            % Example: C=HEAD1.updateCell_fromStruct(S,C)
            
            Col = 1;
            FN   = fieldnames(S);
            Nfs  = numel(FN);  % number of fields in struct
            Nkey = size(C,1);
            if (Nfs==Nkey)
                % assume there is no need to update
            elseif (Nfs<Nkey)
                error('Bug: Nfs<Nkey');
            else
                % Nfs>Nkey
                % Look for new Key values that donot appear in Header
                Cl = HEAD1.replace_iligal_charCell(C);
                [NewKey, NewKeyInd] = setdiff(FN,Cl(:,Col));
                SC        = struct2cell(S);
                NewHeader = [NewKey(:), SC(NewKeyInd)];
                NewHeader = HEAD1.fixColCell(NewHeader);
                C         = [C; NewHeader];
                
            end
                    
            
            
        end
        
    end
    
    % dictionaries
%     methods (Static)
%         function Dict=dictionary
%             % get keywords dictionary
%         end
%     end  % end methods
    
    
    % override methods: isfield, isstruct, isempty
    methods 
        
        %--------------------------
        %--- Structre functions ---
        %--------------------------
        function obj=isfield(Head,Field)
            % isfield 
            obj = any(strcmp(fieldnames(Head),Field));
        end

        function obj=isstruct(Head)
            % isstruct
            obj = true;  %isstruct(Sim) || isa(Sim,'SIM');
        end
        
        function Res=isempty(H)
            % isempty override on HEAD object
            % Package: @HEAD
            % Description: Check if HEAD object isempty.
            %              Return a logical for each HEAD object element.
            %              The logical indicates if the Header field is
            %              empty.
            % Input  : - An HEAD object.
            % Output : - A logical for each HEAD object element.
            %            The logical indicates if the Header field is
            %            empty.
            % Example: R=isempty(H);
            % Reliable: 2
            
            Res = false(size(H));
            N = numel(H);
            for I=1:1:N
                Res(I) = isempty(H(I).Header);
            end
            
        end

    end
    
    
    
    % set/ get keyword data - read data from header
    % disp
    % copy
    % updateHeader_fromStruct
    % delKey, uniqueKey, addKey, replaceKeyName, eplaceKeyVal
    % getKey, getVal_best (previously getkey_fromlist), mgetkey (like old)
    % regexp (like old), regexprep (like old)
    % lowerKey (previously lower_key), upperKey (previously upper_key)
    % numKey (previously numkey)
    % spacedel (new), val2num (new)
    % cell2struct (new)
    % isKeyExist (new), isKeyVal (previoisly iskeyval)
    % getkey  - obsolete
    % iskeyval - obsolete
    % add_end
    methods
        
        function disp(H,UseDisp)
            % display HEAD object Header keywords, values, and comments.
            % Package: @HEAD1
            % Description: Display HEAD object Header keywords, values,
            %              and comments. If multiple elements, then display
            %              all headers.
            % Input  : - An HEAD object.
            %          - A logical flag indicating if to use the disp
            %            command (true), or just the variable name without
            %            ";" (false).
            %            Default is true.
            % Output * null
            % Example: H=HEAD1.basic; disp(H), disp(H,false)
            % Relauble: 2
            
            if nargin<2
                UseDisp = true;
            end
            
            Nh = numel(H);
            for Ih=1:1:Nh
                if UseDisp
                    disp(H(Ih).Header)
                else
                    H(Ih).Header
                end
            end
            
        end
        
        function H=copy(H,C)
            % Populate an HEAD object elements with a Cell or HEAD object.
            % Package: @HEAD
            % Description: Copy a cell header or a single-element HEAD
            %              object content to all the elements in an HEAD
            %              object.
            % Input  : - An HEAD object.
            %          - A cell header (i.e., {key, val, comment}), or a
            %            single-element HEAD object.
            % Output : - An HEAD object in which the second input argument
            %            is copied into all its elements.
            % Example: copy(H,{'TYPE',3,''})
            %          H=HEAD(2,2); copy(H,H(1))
            
            if (HEAD1.isHEAD(C))
                IsHEAD = true;
                if numel(C)>1
                    error('Second input arg must be a cell or a single-element HEAD object');
                end
            else
                IsHEAD = false;
            end
            
            N = numel(H);
            for I=1:1:N
                if IsHEAD
                    H(I) = C;
                else
                    H(I).Header = HEAD1.fixColCell(C);
                end
            end
            
            
        end
        
        function H=updateHeader_fromStruct(H)
            % update Header in HEAD object from Keys structure
            % Package: @HEAD
            % Description: Given an HEAD object, in each element look for
            %              keywords in the "Keys" property structure that
            %              do not appear in the "Header" property cell, and
            %              add them to the "Header" property cell array.
            % Input  : - An HEAD object.
            % Output : - An updated HEAD object.
            % Example: H=HEAD1.basic; H.Keys.A = 1; H1=[H,H];
            %          H2=updateHeader_fromStruct(H1)
            % Reliable: 2
            
            Nh = numel(H);
            for Ih=1:1:Nh
                H(Ih).Header = HEAD1.updateCell_fromStruct(H(Ih).Keys,H(Ih).Header);
            end
        
        end
        
        function H=delKey(H,KeyName,Exact)
            % Delete keywords from Header object
            % Package: @HEAD
            % Description: Delete keywords from all the elements of an
            %              Header object. The search can be exact or
            %              non-exact.
            % Input  : - HEAD object.
            %          - String or cell array of strings of keyword names
            %            to search and delete.
            %          - Exact keyword name matching. Default is true.
            % Output : - HEAD object with the the deleted lines.
            % Example: H=delKey(H,'TYPE');
            % Reliable: 2
            
            if (nargin<3)
                Exact = true;
            end
            
            N = numel(H);
            for I=1:1:N
                H(I).Header = HEAD1.delKeyCell(H(I).Header,KeyName,Exact);
            end
        end
        
        function H=uniqueKey(H,varargin)
            % Select only uniuqe keywords in an HEAD object
            % Package: @HEAD
            % Description: Given an HEAD object, in each element save only
            %              the unique keyword names.
            % Input  : - An HEAD object.
            %          - In case of duplicate keywords, which cccurance 
            %            to keep: 'first' | 'last'. Default is 'first'.
            % Example: H=uniqueKey(H)
            
            N = numel(H);
            for I=1:1:N
                H(I).Header = HEAD1.uniqueKeyCell(H(I).Header,varargin{:});
            end
            
        end
        
        function H=addKey(H,New,varargin)
            % Add keywords to an HEAD object
            % Package: @HEAD
            % Description: Add keywords to all the elements in an HEAD
            %              object. The keywords can be a single-element
            %              HEAD object or a 2 or 3 columns cell array of
            %              {keyword, value, comment}.
            % Input  : - An HEAD object.
            %          - A single-element HEAD object or a 2 or 3 columns
            %            cell array of {keyword, value, comment}.
            %          - Location in which to add the keywords.
            %            Inf for adding at the end (but before the END
            %            keyword), 0 for adding
            %            at the begining. Default is Inf.
            % Output : - AN HEAD object.
            % Example: addKey(H,H(1))
            %          addKey(H,{'A',2},2)
            % Reliable: 2
            
            if (HEAD1.isHEAD(New))
                if (numel(New)>1)
                    error('Second argument must be a cell header or a single-element HEAD object');
                end
                Cell = New.Header;
            else
                Cell = HEAD1.fixColCell(New);
            end
                
            N = numel(H);
            for I=1:1:N
                H(I).Header = HEAD1.addKeyCell(H(I).Header,Cell,varargin{:});
            end
            
            
        end
        
        function H=replaceKeyName(H,OldKeyName,NewKeyName)
            % Replace keyword name in HEAD object
            % Package: @HEAD
            % Description: Replace keyword name in all elements in an
            %              HEAD object.
            %              If keyword name is not found than skip the
            %              replecment.
            %              The replacment is applied for all HEAD elements.
            % Input  : - AN HEAD object.
            %          - A string or a cell array of strings containing
            %            the old keyword names to search in the header.
            %          - Astring or a cell array of strings containing
            %            new keyword names that will replace the old names.
            %            This input should have the same 
            %            number of keywords as the second input arguments.
            % Output : - An HEAD object with the new keyword names.
            % Example: H1=replaceKeyName(H,{'a','RA'},{'DECDEG','RADEG'})
            % REliable: 2
            
            Nh = numel(H);
            for Ih=1:1:Nh
                H(Ih).Header = HEAD1.replaceKeyNameCell(H(Ih).Header,OldKeyName,NewKeyName);
            end
        end
            
        function H=replaceKeyVal(H,KeyName,KeyVal)
            % Replace keyword value in HEAD object
            % Package: @HEAD
            % Description: Replace keyword value in all elements in an
            %              HEAD object.
            %              If keyword name is not found than skip the
            %              replecment.
            %              The replacment is applied for all HEAD elements.
            % Input  : - An HEAD object.
            %          - A string or a cell array of strings containing
            %            the keyword names to search in the header,
            %            and for which to replace the values.
            %          - Keyword values. The same number of values as
            %            keyword names.
            %            This can be a scalar, a cell array or a vector.
            % Output : - An HEAD object with the new keyword names.
            % Example: H1=replaceKeyVal(H,{'a','RA'},{'best',102.12})
            % REliable: 2
            
            Nh = numel(H);
            for Ih=1:1:Nh
                H(Ih).Header = HEAD1.replaceKeyValCell(H(Ih).Header,KeyName,KeyVal);
            end
            
            
        end
        
        function Out=getKey(Head,KeyName,varargin)
            % Get keyword value from object header
            % Package: HEAD
            % Description: Search keyeod names in an HEAD object and return
            %              a new HEAD object with the requested keywords.
            %              The search can be exact or by using regular
            %              expression. Use getFilledKey in order to return
            %              filled values for non-existing keywords.
            % Input  : - HEAD object
            %          - Keyword name, or a cell array of keyword names.
            %          - A flag indicating if to perform an exact name
            %            search (true). Default is true.
            %            If false, then use regexp.
            %          - Column in the cell array in which to perform the
            %            search. Default is 1 (i.e., the keyword column).
            %            Note that if 2 is used, than this may fail if the
            %            some values are numeric.
            %          - If more than one keyword was matched, this is the
            %            number of first keywords/values to return.
            %            If Inf, then retun all keywords.
            %            Default is Inf.
            % Output : - An HEAD object containing only lines containing
            %            the requested keywords.
            % Example: H=HEAD1.basic; K=getKey(H,'TYPE')
            
            N = numel(Head);
            % for each element in object 
            Out = HEAD1(size(Head));
            for I=1:1:N
                Out(I).Header = HEAD1.getKeyCell(Head(I).Header,KeyName,varargin{:});
            end
            
        end  % end getKey function
        
        function Out=getFilledKey(Head,KeyName,varargin)
            % Get filled keyword value from object header
            % Package: HEAD
            % Description: Perform exact keyeod name search in an HEAD
            %              object and return
            %              a new HEAD object with the requested keywords,
            %              and filled values in case the keyword doesn't
            %              exist. Use getKey for not filled/non exact
            %              search.
            % Input  : - HEAD object
            %          - Keyword name, or a cell array of keyword names.
            %          - Value to populate the keyword value in case the
            %            keyword doesn't exist. Default is NaN.
            %          - If more than one keyword was matched, this is the
            %            index of the keywords/values to return.
            %            (in anycase only one line will be returned).
            %            Inf will return the last keyword.
            %            Default is 1.
            %          - Column in the cell array in which to perform the
            %            search. Default is 1 (i.e., the keyword column).
            %            Note that if 2 is used, than this may fail if the
            %            some values are numeric.
            % Output : - An HEAD object with the lines containing the
            %            searched keywords. The number of lines always
            %            equal to the number of searched keywords and
            %            appears in the same order. If the keywords doesn't
            %            exist in the header than a filled value is
            %            returned.
            % Example: H=HEAD1.basic; K=getFilledKey(H,'TYPE')
            %          K=getFilledKey(H,{'TYPE','B'})
            
            N = numel(Head);
            % for each element in object 
            Out = HEAD1(size(Head));
            for I=1:1:N
                Out(I).Header = HEAD1.getFilledKeyCell(Head(I).Header,KeyName,varargin{:});
            end
            
        end
      
        function [CellVal,SubHead] = mgetkey(Head,KeyName,SpaceDel,Conv2Num)
            % Get the value of multiple keywords from a multi-elemet HEAD
            % Package: @HEAD
            % Description: Get the value of a multiple keywords from a
            %              mult-element HEAD object. Keywords that doesn't
            %              exist will be replaced with NaN. The output is a
            %              cell arraey of size Number HEAD elements by
            %              Number of keywords.
            %              Optionaly delete spaces from strings and convert
            %              string to numbers if possible.
            % Input  : - An HEAD object.
            %          - A string or a cell array of string of keyword
            %            names.
            %          - A logical indicating if to delete leading and
            %            trailing spaces. Default is true.
            %            implemented using HEAD/spacedel.
            %          - A logical indicating if to convert strings to
            %            numbers if possible. Default is false.
            %            Implemented using HEAD/val2num
            % Output : - A cell array of size Number HEAD elements by
            %            Number of keywords, containing the requested
            %            keyword values. NaN if keyword doesn't exist.
            %          - An Haed  object with the requested keywords only,
            %            and missing keywords filled by NaN values.
            % Example: H=HEAD1.basic; mgetkey([H;H],{'TYPE','B'})
            % Reliable: 2
            
            Col = 2;
            
            if nargin<4
                Conv2Num = false;
                if nargin<3
                    SpaceDel = true;
                end
            end
            
            if ischar(KeyName)
                KeyName = {KeyName};
            end
            
            SubHead = getFilledKey(Head,KeyName,NaN,1,1);
            if SpaceDel
                SubHead = spacedel(SubHead,'leadtrail');
            end
            
            if (Conv2Num)
                SubHead = val2num(SubHead);
            end
            
            Nkey = numel(KeyName);
            Nh   = numel(SubHead);
            CellVal = cell(Nh,Nkey);
            for Ih=1:1:Nh
                CellVal(Ih,:) = SubHead(Ih).Header(:,Col).';
            end
            
            
        end
        
        function Out=regexp(H,Col,varargin)
            % Execture regexp on single element HEAD object
            % Package: @HEAD
            % Description: Execute regexp.m on one of the columns of an HEAD object
            %                with a single element.
            % Input  : - A a single element HEAD object.
            %          - Column index on which to execute regexp.m (1|2|3).
            %          * Additional arguments to pass to regexp.m, ...,EXPRESSION,...
            % Output : - The output of regexp.m
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Apr 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: regexp(Head,1,'A_\d+_\d+')
            % Reliable: 2

            if numel(H)>1
                error('Input must be a single element HEAD object');
            end
                        
            Out = regexp(H.Header(:,Col),varargin{:});
        end
        
        function H=regexprep(H,Col,varargin)
            % Execture regexprep on single element HEAD object
            % Package: @HEAD
            % Description: Execute regexprep.m on one of the columns of an
            %              HEAD object.
            % Input  : - An HEAD object.
            %          - Column index on which to execute regexp.m (1|2|3).
            %          * Additional arguments to pass to regexprep.m,
            %            ...,EXPRESSION,...
            % Output : - An HEAD object with the replaced strings.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Apr 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: regexprep(Head,1,'A','B')
            % Reliable: 2

            N = numel(H);
            for I=1:1:N
                H(I).Header(:,Col) = regexprep(H(I).Header(:,Col),varargin{:});
            end
        end
        
        function H=lowerKey(H)
            % Convert all keyword names to lower case
            % Package: @HEAD
            % Description: Convert all keyword names to lower case
            % Input  : - An HEAD object
            % Output : - An HEAD object in which the key names are lower
            %            case.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Head=lowerKey(Head)
            % Reliable: 2
            
            Col = 1;
            for I=1:1:numel(H)
                H(I).Header(:,Col) = lower(H(I).Header(:,Col));
            end
            
        end
        
        function H=upperKey(H)
            % Convert all keyword names to upper case
            % Package: @HEAD
            % Description: Convert all keyword names to upper case
            % Input  : - An HEAD object
            % Output : - An HEAD object in which the key names are upper
            %            case.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Head=upperKey(Head)
            % Reliable: 2
            
            Col = 1;
            for I=1:1:numel(H)
                H(I).Header(:,Col) = upper(H(I).Header(:,Col));
            end
            
        end
            
        function Num=numKey(H)
            % Count the number of keywords in HAED object
            % Package: @HEAD
            % Description: Count the number of keywords in HAED object
            % Input  : - An HEAD object
            % Output : - An array in which each element contains the number
            %            of keywords in each HEAD element.
            % Example: numKey(H)
            % Reliable: 2
            Num = zeros(size(H));
            N = numel(H);
            for I=1:1:N
                Num(I) = size(H.Header,1);
            end
            
        end
        
        function H=spacedel(H,Type,Col)
            % Deleted spaces from all character values
            % Package: @HEAD
            % Description: Delete spaces from all values in an HEAD object.
            %              Can deal with all spaces, leading of trailing.
            % Input  : - An HEAD object.
            %          - Space removal option:
            %            'all' - remove all spaces.
            %            'trail' - remove trailing spaces.
            %            'leadtrial' - remove leading and trailing spaces.
            %                          Default.
            %            'lead' - Remove leading spaces.
            %          - Header column on which to oerate.
            %            Default is 2 (i.e., Values column).
            % Output : - An HEAD object.
            % Example: H=spacedel(H)
            % Reliable: 2
           
            if nargin<3
                Col = 2;
                if nargin<2
                    Type = 'leadtrail';
                end
            end
            
            N = numel(H);
            for I=1:1:N
                Flag = cellfun(@ischar,H(I).Header(:,Col)); % char values
                
                switch lower(Type)
                    case 'all'
                        H(I).Header(Flag,Col) = Util.string.spacedel(H(I).Header(Flag,Col));
                    case 'trail'
                        H(I).Header(Flag,Col) = deblank(H(I).Header(Flag,Col));
                    case 'leadtrail'
                        H(I).Header(Flag,Col) = strip(H(I).Header(Flag,Col),'both');
                    case 'lead'
                        H(I).Header(Flag,Col) = strip(H(I).Header(Flag,Col),'left');
                    otherwise
                        error('Unknown Typr option');
                end
            end
            
        end
        
        function H=val2num(H)
            % Convert values in HEAD object to numeric if possible
            % Package: @HEAD
            % Description: Convert values in HEAD object from strings to
            %              numeric if possible.
            % Input  : - An HEAD object.
            % Output : - An HEAD object in which a numeric values in
            %            strings are converted to numbers.
            %            A string like 'NaN' will not be converted to NaN.
            % Example: H=HEAD1.basic; H=val2num(H);
            % Reliable: 2 
            
            
            Nh = numel(H);
            for Ih=1:1:Nh
                FlagStr = cellfun(@ischar,H(Ih).Header(:,2));
                NumVal  = cellfun(@str2double,H(Ih).Header(:,2));
                Flag    = ~isnan(NumVal) & FlagStr;
                H(Ih).Header{Flag,2} = NumVal(Flag);
            end
                
            
            
        end
        
        function H=replace_iligal_char(H,Col)
            % replace iligal characters in keyword names (e.g. '-').
            % Package: @HEAD
            % Description: Given an HEAD object. Replace iligal
            %              characters in the first column.
            %              List of illigal characters:
            %              '-' replace with '_'
            % Input  : - An HEAD object.
            %          - Column index in which to replace characaters.
            %            Default is 1.
            % Output : - An HEAD object.
            % Example: H=HEAD1.basic; H=H.addKey({'UTC-OBS','2009'});
            %          H=replace_iligal_char(H)
            % Reliable: 2
            
            if nargin<2
                Col = 1;
            end
            
            Nh = numel(H);
            for Ih=1:1:Nh
                H(Ih).Header = HEAD1.replace_iligal_charCell(H(Ih).Header,Col);
            end            
            
        end % end replace_iligal_char function
        
        function S=cell2struct(H,Unique)
            % Convert the Header fields in HEAD object to structure
            % Package: @HEAD
            % Description: Convert the Header fields in HEAD object to
            %              structure array in which the field name is the
            %              keyword name and the field value is the keyword
            %              value. Optionally run uniqueKey on the HEAD
            %              object.
            %              Before the conversion the function replace
            %              iligal characters in the keyword names.
            % Input  : - An HEAD object.
            %          - A logical falg indicating if to run uniqueKey on
            %            the HEAD object prior to conversion to structure.
            %            Default is true.
            % Output : - A structure in which the field names are the
            %            keywords.
            % Example: cell2struct(HEAD1.basic)
            % Relible: 2
            
            
            if (nargin<2)
                Unique = true;
            end
            
            if Unique
                H = uniqueKey(H);
            end
           
            N = numel(H);
            for I=1:1:N
                H(I).Header = HEAD1.replace_iligal_charCell(H(I).Header);
                S(I) = cell2struct(H(I).Header(:,2),H(I).Header(:,1));
            end
            
        end
        
        function [Val,Out] = getVal_best(H,KeyList,SelectMethod,ExactDict,FillVal)
            % Select the (first) existing header line out of many keywords
            % Package: @HEAD
            % Description: Given an HEAD object and a list of keywords
            %              select only one available keyword (e.g., the
            %              first available keyword). This function is
            %              useful when the exact name of keywprd is
            %              unknown. The list of keywords is a cell array of
            %              strings. However, if a single string is given
            %              then the function search for this string in the
            %              dictKey dictionary function that return a cell
            %              array of synonsyms.
            % Input  : - An HEAD object.
            %          - A cell array (or a single string) or keywords to
            %            test. E.g., {'TYPE','IMTYPE','IMGTYPE'}.
            %            Alternatively this can be a string (e.g., 'TYPE').
            %            In this case, then will use the HEAD.dictKey
            %            command to retrieve the list of synonyms to this
            %            keywords, and these synonyms will be used as a
            %            KeyList.
            %          - Keyword selsection option:
            %            'first' - select the first keyword available out
            %                      of the keyword list (second input 
            %                      argument), by their order in the list
            %                      (rather than order in the header).
            %            'last'  - Like first but for the 'last' keyword.
            %          - A logical indicating if to use exact search (true)
            %            when the dictionary is used. Default is true.
            %          - Fill value in the cell array output.
            %            Default is NaN.
            % Output : - A cell array of keyword values. The size of the
            %            cell array is equal to the size of the HEAD
            %            object. Each cell element contains a keyword value
            %            (string or numeric). If the keyword doesn't exist,
            %            then the corresponding cell element contains the
            %            fill value.
            %          - An HEAD object with only one line in the Header of
            %            each HEAD element (i.e., the selected keyword).
            %            If Keyword not found then return an empty header.
            % Example: H=HEAD1.basic;
            %          [V,O]=getVal_best([H,H],{'YY','IMTYPE','TYPE'});
            %          [V,O]=getVal_best([H,H],'TYPE'); % in this option the
            %                                   dictionary replaces 'TYPE' with a cell array
            %                                   {'TYPE','OBSTYPE',...}.
            %          [V,O]=getVal_best(H,'RA');
            %          [V,O]=getVal_best(H,{'JRA'}); % return empty
            %          [V,O]=getVal_best(H,'RA');
            %          [V,O]=getVal_best(H,'JRA'); % returns the RA keyword
            % Reliable: 2
            
            ColVal = 2;
            %if nargin<6
            %    Conv2Num = true;
            if nargin<5
                FillVal = NaN;
                if nargin<4
                    ExactDict = true;
                    if nargin<3
                        SelectMethod = 'first';
                    end
                end
            end
            %end
            
            
            if ischar(KeyList)
                Dict = HEAD1.dictKey(KeyList,ExactDict);
                if isempty(Dict)
                    error('Key: %s was not found in dictionary',KeyList);
                end
                
                if isempty(KeyList)
                    error('Key: %s was not found in dictionary',KeyList);
                end
                KeyList = Dict(1).Synonym;
            end
            
            N  = numel(H);
            Out = getKey(H,KeyList);
            Val = cell(size(H));
            for I=1:1:N
                if ~isempty(Out(I).Header)
                    switch lower(SelectMethod)
                        case 'first'
                            Out(I).Header = Out(I).Header(1,:);
                        case 'last'
                            Out(I).Header = Out(I).Header(end,:);
                        otherwise
                    end
                    Val{I} = Out(I).Header{1,ColVal};
                    
                else
                    % fill value only for Val cell array
                    Val{I} = FillVal;
%                     if Conv2Num && ischar(Val{I})
%                         Val{I} = str2double(Val{I});
%                     end
                        
                end
            end
            
            
        end
        
        function [IsExist,Counter] = isKeyExist(H,Key,varargin)
            % Check (and count) if a keyword name exist in an HEAD object
            % Package: @HEAD
            % Description: Check if keyword name exist in an HEAD object,
            %              and also count the number of appearnces for each
            %              keyword. Can either use exact match or regular
            %              exprssions.
            % Input  : - A cell header {Keyword, Value, Comment}
            %          - A keyword string or a cell array of keyword
            %            strings to search.
            %          - A flag indicating if to perform exact match (true)
            %            or regular exprssion matching (false).
            %          - Column index in the cell header on which to
            %            perfoprm the search. Default is 1.
            % Output : - A matrix of logicals indicating if each keyword
            %            exist in the header.
            %            Rows represents different HEAD elements, while
            %            columns repesents different keywords.
            %          - A matrix of counters of the number of
            %            appearnces of each keyword in the header.
            %            Rows represents different HEAD elements, while
            %            columns repesents different keywords.
            % Example: H=HEAD1.basic;
            %          [IsExist,Counter]=isKeyExist([H;H],{'TYPE','N','EXPTI'},false)
            
            if ischar(Key)
                Key = {Key};
            end
            Nkey = numel(Key);
            
            Nh = numel(H);
            IsExist = false(Nh,Nkey);
            Counter = zeros(Nh,Nkey);
            for Ih=1:1:Nh
                [IsExist(Ih,:),Counter(Ih,:)] = HEAD1.isKeyExistCell(H(Ih).Header,Key,varargin{:});
            end
            
        end
        
        function IsEqual=isKeyVal(H,Key,Val,FillVal)
            % Compare the keyword values in HEAD object with some vcalues
            % Package: @HEAD
            % Descrption: Given an HEAD object, keyword names, and vlaues,
            %             search for the keword names in the headers and
            %             compare their values with the function input
            %             values (3rd input argument). Return a matrix of
            %             logicals. The number of rows equal to the number
            %             of elements in the HEAD object and the number of
            %             columns equal to the number of keyword names. The
            %             matrix is poulated with true, where the keyword
            %             value equal to the input value.
            %             Note the comparison is extended to NaNs and
            %             empties.
            % Input  : - An HEAD object.
            %          - String of a keyword name, or a cell array of
            %            strings of keyword names.
            %          - Cell array of string/numeric values, or vector of
            %            numeric values. The number of values corresponds
            %            to the number of keyword names. Each value is
            %            compared with the keyword value in the header.
            %          - Fill value in case the keyword name doesn't exist.
            %            Default is NaN.
            % Output : - A matrix of logicals.
            %            The number of rows equal to the number
            %             of elements in the HEAD object and the number of
            %             columns equal to the number of keyword names. The
            %             matrix is poulated with true, where the keyword
            %             value equal to the input value.
            % Example: H=HEAD1.basic; IE=isKeyVal(H,'TYPE','bias')
            %          IE=isKeyVal([H;H],{'TYPE','EXPTIME'},{'bias','a'})
            %          IE=isKeyVal([H;H],{'TYPE','EXPTIME'},[1 1])      
            
            if nargin<4
                FillVal = NaN;
            end
            
            Col = 2;
            if ischar(Key)
                Key = {Key};
            end
            
            if iscell(Val)
                % do nothingg
            elseif ischar(Val)
                Val = {Val};
            elseif isnumeric(Val)
                Val = num2cell(Val);
            else
                error('Val (3rd) argument must be numeric, char, or cell of char/numeric');
            end
            
            H = getFilledKey(H,Key,FillVal);
            
            Nkey = numel(Key);
            Nh   = numel(H);
            IsEqual = false(Nh,Nkey);
            for Ih=1:1:Nh
                for Ikey=1:1:Nkey
                    if ischar(H(Ih).Header{Ikey,Col}) && ischar(Val{Ikey})
                        IsEqual(Ih,Ikey) = strcmp(H(Ih).Header{Ikey,Col},Val{Ikey});
                    elseif isnumeric(H(Ih).Header{Ikey,Col}) && isnumeric(Val{Ikey})
                        if isempty(H(Ih).Header{Ikey,Col}) || isempty(Val{Ikey})
                            if isempty(H(Ih).Header{Ikey,Col}) && isempty(Val{Ikey})
                                % treat empties
                                IsEqual(Ih,Ikey) = true;
                            else
                                IsEqual(Ih,Ikey) = false;
                            end
                        elseif isnan(H(Ih).Header{Ikey,Col}) && isnan(Val{Ikey})
                            % treat NaNs
                            IsEqual(Ih,Ikey) = true;
                        else
                            IsEqual(Ih,Ikey) = H(Ih).Header{Ikey,Col} == Val{Ikey};
                        end
                    else
                        % KeyVal and Val are not of the same type
                        IsEqual(Ih,Ikey) = false;
                    end
                end
            end
                         
             
        end
        
        function H=add_end(H)
            % Add END keyword at the end of heaer
            % Package: @HEAD
            % Description: Add END keyword at the end of each header in an
            %              HEAD object. If the END keyword appears in the
            %              middle of an header it will be removed.
            % Input  : - An HEAD object.
            % Output : - An HEAD object.
            % Example: H=HEAD1.basic; H=add_end(H)
            % Reliable: 2
            
            Nh = numel(H);
            for Ih=1:1:Nh
                switch H(Ih).Header{end,1}
                    case 'END'
                        % end appear in end of header do nothing
                    otherwise
                        % add end
                        H(Ih).Header{end+1,1} = 'END';
                end
                
                % verify there is no end in the middle of the header
                Flag = [strcmp(H(Ih).Header(1:end-1,1),'END'); false];
                if any(Flag)
                    H(Ih).Header = H(Ih).Header(~Flag,:);
                end
                
            end
            
        end

    end  % methods
    
    
    % get/set special keywords
    % naxis, getCoo, getObsCoo, julday
    % find_groups
    % getType, isType, isbias, isflat, isdark, isarc
    methods
        function Naxis=naxis(H)
            % Get naxis realted keywords from HEAD object.
            % Package: @HEAD
            % Description: Get the value of the NAXIS1, NAXIS2,... keywords
            %              from an HEAD object.
            % Input  : - An HEAD object.
            % Output : - A matrix of NAXIS value keywords. each line corresponds to
            %            on HEAD element (or e.g., a SIM image). The number of columns
            %            is equal to the maximum NAXIS value, and the columns
            %            corresponds to NAXIS1, NAXIS2, etc. NaN for
            %            non-existing values.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: ImSize=naxis(S)
            % Reliable: 2

            ColName = 1;
            ColVal  = 2;
            SubH = getKey(H,'NAXIS\w+',false); % non exact keyword search
            Nkey = numKey(SubH);  % number of keywords in each element
            
            Nh       = numel(SubH);
            Naxis    = nan(Nh,max(Nkey));
            for Ih=1:1:Nh
                for Ikey=1:1:Nkey(Ih)
                    DimIndex = str2double(SubH(Ih).Header{Ikey,ColName}(6:end));
                    Naxis(Ih,DimIndex) = SubH(Ih).Header{Ikey,ColVal};
                end
            end
            

        end
        
        % find_groups
        function Groups=find_groups(Head,Keys,DelSpace)
            % group headers by unique keyword values
            % Package: @HEAD
            % Description: Find groups of headers which have keywords with
            %              identical values. Given an HEAD object with
            %              multiple elements and a list of keywords, find
            %              unique groups defined as having the same keyword
            %              values. For example, images which were taken
            %              with the same filter and same exposure time.
            %              This function use the
            %              Util.cell.cell_find_groups.m function.
            % Input  : - An HEAD object.
            %          - Cell array of keywords by which to find groups.
            %          - Delete spaces from strings prior to comparison
            %            {true|false}.
            %            Default is true.
            % Output : - Structure of groups. Each structure represent a group of rows
            %            in the input cell array which have equal values.
            %            The structure contains the following fields:
            %            .Conntent  - Cell array of values define the group.
            %            .ptr       - The indices of the rows that belong to the
            %                         group.
            % See also: cell_find_groups.m
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Groups=find_groups(Sim,{'EXPTIME','FILTER'});
            %          Groups=find_groups(lower_key(Sim),{'exptime','filter'});
            % Reliable: 2

            Cell = mgetkey(Head,Keys);

            % find groups:
            Groups = Util.cell.cell_find_groups(Cell,DelSpace);

            
        end
        
        function [RA,Dec,Equinox]=getCoo(H,varargin)
            % get RA/Dec/Equinox from HEAD object and convert to deg
            % Package: @HEAD
            % Description: Search for RA/DEC/EQUINOX keywords in HEAD
            %              object. Of they does not exist then return NaN.
            %              If they are numeric then return as is.
            %              If RA/DEC are strings then assume they are in
            %              sexagesimal format and convert to degrees.
            %              If EQUINOX is string then it checks if the first
            %              character is J (Julian) or B (Besselian) and
            %              convert to Julian years. If J or B are not
            %              indicated then assume "J".
            % Input  : - An HEAD object.
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'KeyRA' - R.A. keyword name. If string
            %                      then will use the dictKey and getVal_best
            %                      functions to find possible aliases.
            %                      If cell array of strings then look for
            %                      the first keyword that appaers in the
            %                      Header.
            %                      Default is 'RA'.
            %            'KeyDec' - Like 'KeyRA', but for the J2000.0
            %                      Declination. Default is 'DEC'.
            %            'KeyEquinox' - Like 'KeyRA', but for the
            %                      Equinox. Default is 'EQUINOX'.
            % Output : - A numeric array of R.A. (likely in
            %            degrees). The size of the array is as the size of
            %            the HEAD object.
            %          - A numeric array of Dec.
            %          - A numeric array of Equinox. Likely in Julian
            %            years.
            % Example: [RA,Dec,Equinox]=getCoo(H)
            % Reliable: 2
            
            ColVal = 2;
            
            DefV.KeyRA                = 'RA';
            DefV.KeyDec               = 'DEC';
            DefV.KeyEquinox           = 'EQUINOX';
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            
            Val_RA      = getVal_best(H,InPar.KeyRA);
            Val_Dec     = getVal_best(H,InPar.KeyDec);
            Val_Equinox = getVal_best(H,InPar.KeyEquinox);
            
            % convert RA to requested units
            Nh = numel(H);
            RA      = nan(size(H));
            Dec     = nan(size(H));
            Equinox = nan(size(H));
            for Ih=1:1:Nh
                if isnumeric(Val_RA{Ih})
                    RA(Ih) = Val_RA{Ih};
                elseif ischar(Val_RA{Ih})
                    RA(Ih)  = celestial.coo.convertdms(Val_RA{Ih},'gH','d');
                else
                    error('Unknown RA value type');
                end
                
                if isnumeric(Val_Dec{Ih})
                    Dec(Ih) = Val_Dec{Ih};
                elseif ischar(Val_Dec{Ih})
                    Dec(Ih)  = celestial.coo.convertdms(Val_Dec{Ih},'gD','d');
                else
                    error('Unknown RA value type');
                end
            
                if isnumeric(Val_Equinox{Ih})
                    Equinox(Ih) = Val_Equinox{Ih};
                elseif ischar(Val_Equinox{Ih})
                    switch lower(Val_Equinox{Ih}(1))
                        case 'j'
                            Equinox(Ih) = str2double(Val_Equinox{Ih}(2:end));
                        case 'b'
                            Equinox(Ih) = str2double(Val_Equinox{Ih}(2:end));
                            % convert to julian years
                            Equinox(Ih) = convert.time(Val_Equinox{Ih}, upper(Val_Equinox{Ih}(1)),'J');
                        otherwise
                            % assuming Equnox is in Julian years
                            Equinox(Ih) = str2double(Val_Equinox{Ih});
                    end
                else
                    error('Unknown RA value type');
                end
            end
                
        end
          
        function [Long,Lat,Height]=getObsCoo(H,varargin)
            % Get observatory geodetic coordinates from HEAD object
            % Package: @HEAD
            % Description: Get observatory geodetic coordinates from HEAD
            %              object. Observatory long/lat/height are searched
            %              using getVal_best and dictKey dictionary.
            % Input  : - An HEAD object.
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'KeyLong' - Observatory longitude dictionary
            %                       leyword. Default is 'OBSLON'.
            %            'KeyLat'  - Observatory latitude dictionary
            %                       leyword. Default is 'OBSLAT'.
            %            'KeyHeight'- Observatory height dictionary
            %                       leyword. Default is 'HEIGHT'.
            % Output : - Matrix of observatory longitude. Each element
            %            corresponds to HEAD object element.
            %            NaN if value not found, or the value is a string
            %            and can not converted to numeric.
            %          - Matrix of observatory latitude.
            %          - Matrix of observatory height.
            % Example: getObsCoo(H)
             
            ColVal = 2;
            
            DefV.KeyLong              = 'OBSLON';
            DefV.KeyLat               = 'OBSLAT';
            DefV.KeyHeight            = 'HEIGHT';
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            
            Long      = getVal_best(H,InPar.KeyLong);
            Lat       = getVal_best(H,InPar.KeyLat);
            Height    = getVal_best(H,InPar.KeyHeight);
            
            % convert to numeric values - strings to NaN
            Long   = HEAD1.conv2numCell(Long,false);
            Lat    = HEAD1.conv2numCell(Lat,false);
            Height = HEAD1.conv2numCell(Height,false);
            
        end
        
        function [MidJD,ExpTime]=julday(H)
            % Calculate mid exposure JD and ExpTime for HEAD object
            % Package: @HEAD
            % Description: Given an HEAD object calculate the mid exposure
            %              Julian day (JD) and exposure time for each
            %              header. The JD is calculated by attempting
            %              reading various time related keywords starting
            %              with JD, MJD, UTC-OBS etc.
            % Input  : - An HEAD object.
            % Outout : - Matrix of mid exposure time JD, for each HEAD
            %            element. NaN if can not calvulate JD.
            %          - Matrix of exposure times (usually seconds), for
            %            each HEAD element. NaN if not found.
            % Example: H=HEAD1.basic; [JD,E]=julday(H)
            SEC_IN_DAY = 86400;
            
            
            
            ExpTime      = getVal_best(H,'EXPTIME');  % s
            JD           = getVal_best(H,'JD');  % day
            MIDJD        = getVal_best(H,'MIDJD');  % day
            MJD          = getVal_best(H,'MJD');  % day
            MIDMJD       = getVal_best(H,'MIDMJD');  % day
            Date         = getVal_best(H,'UTC-OBS');  % day
            
            
            MidJD = nan(size(Date));
            for I=1:1:numel(Date)

                if ~isnan(MIDJD{I})
                    MidJD(I) = MIDJD{I};
                else
                    if ~isnan(MIDMJD{I})
                        MidJD(I) = convert.time(MIDMJD{I},'MJD','JD');
                    else
                       if ~isnan(JD{I})
                           MidJD(I) = JD{I} + 0.5.*ExpTime{I}./SEC_IN_DAY;
                       else
                           if ~isnan(MJD{I})
                               MidJD(I) = convert.time(MJD{I},'MJD','JD') + 0.5.*ExpTime{I}./SEC_IN_DAY;
                           else
                               if ~isnan(Date{I})
                                   % convert date to JD
                                   MidJD(I) = celestial.time.julday(Date{I}) + 0.5.*ExpTime{I}./SEC_IN_DAY;
                               else
                                   % error - can not read date/JD...
                                   MidJD(I) = MIDJD{I};  % NaN
                               end
                           end
                       end
                    end
                end
            end
            
            
        end
        
        function [ValTrans,ValType]=getType(H)
            % Read synonymous TYPE keyword and translate using dictionary
            % Package: @HEAD
            % Description: Read the TYPE keyword (or synnymous keyword)
            %              and translate using dictValType dictionary.
            % Input  : - An HEAD object
            % Output : - A cell array of the translated TYPE values.
            %            The size of the cell array is identical to the
            %            size of the HEAD object.
            %            Each element contains the TYPE keyword value
            %            translated using HEAD1.dictValType.
            %            For example: 'bias10' value will be translated to
            %            'bias'.
            %          - A cell array of the TYPE values (not translated).
            % Example: H=HEAD1.basic; ImType=getType([H;H])
            % Reliable: 2
            
            ValType      = getVal_best(H,'TYPE');
            ValTrans     = cell(size(ValType));
            Nh = numel(H);
            for Ih=1:1:Nh
                [~,ValTrans{Ih}] = HEAD1.dictValType(ValType{Ih},false);
            end
                 
             
            
        end  % end getType function
        
        function [IsType,TypeValTrans]=isType(H,TypeVal)
            % Is TYPE or synonymous keyword in HEAD object equal value
            % Package: @HEAD
            % Description: Check if TYPE or one of its synonymous keyword 
            %              (as defined in in HEAD1.dictKey) in an HEAD
            %              object equal some value. The comparison to the
            %              value is done against synonyms in the
            %              HEAD.dictValType function.
            % Inpput  : - An HEAD object.
            %           - User input TYPE value (e.g., 'bias').
            %             Legal type values are
            %             listed in HEAD1.dictValType.Key.
            %             The function will return an error if the TYPE
            %             value does not appear in HEAD1.dictValType.
            % Output : - A matrix of logical indicating the the TYPE value
            %            equal the use input TYPE value.
            %          - A cell array of the translated TYPE values.
            %            The size of the cell array is identical to the
            %            size of the HEAD object.
            %            Each element contains the TYPE keyword value
            %            translated using HEAD1.dictValType.
            %            For example: 'bias10' value will be translated to
            %            'bias'.
            % Example: H=HEAD1.basic; isType(H,'bias')
            % Reliable:
            
            TypeValTrans = getType(H);
            
            % check if user TypeVal (TYPE value) is in the dictionary
            DV = HEAD1.dictValType;
            if ~any(strcmp({DV.Key},TypeVal))
                error('TypeVal:%s is not found in dictionary',TypeVal);
            end
            IsType = strcmp(TypeValTrans,TypeVal);
            
            
        end  % end isType function
        
        function Ans=isbias(H,Val)
            % Check if HEAD object TYPE keyword is bias or synonymous
            % Package: @HEAD
            % Description: Check if HEAD object TYPE keyword is bias or
            %              synonymo of bias. The TYPE keyword is selected
            %              using the dictKey function, while the bias value
            %              is selected using the dictValType function.
            % Input  : - An HEAD object.
            %          - The keyword value in the HEAD1.dictValType
            %            function indicating a bias image.
            %            Default is 'bias'.
            % Output : - A matrix of logical indicating if each one of the
            %            HEAD object elements TYPE is consistent with bias.
            % Example: H=HEAD1.basic; Ans=isbias(H)
            % Reliable: 2
            
            if nargin<2
                Val = 'bias';
            end
            
            Ans = isType(H,Val);
            
        end
        
        function Ans=isflat(H,Val)
            % Check if HEAD object TYPE keyword is flat or synonymous
            % Package: @HEAD
            % Description: Check if HEAD object TYPE keyword is flat or
            %              synonymous of flat. The TYPE keyword is selected
            %              using the dictKey function, while the flat value
            %              is selected using the dictValType function.
            % Input  : - An HEAD object.
            %          - The keyword value in the HEAD1.dictValType
            %            function indicating a flat image.
            %            Default is 'flat'.
            % Output : - A matrix of logical indicating if each one of the
            %            HEAD object elements TYPE is consistent with flat.
            % Example: H=HEAD1.basic; Ans=isflat(H)
            % Reliable: 2
            
            if nargin<2
                Val = 'flat';
            end
            
            Ans = isType(H,Val);
            
        end
            
        function Ans=isdark(H,Val)
            % Check if HEAD object TYPE keyword is dark or synonymous
            % Package: @HEAD
            % Description: Check if HEAD object TYPE keyword is dark or
            %              synonymous of dark. The TYPE keyword is selected
            %              using the dictKey function, while the dark value
            %              is selected using the dictValType function.
            % Input  : - An HEAD object.
            %          - The keyword value in the HEAD1.dictValType
            %            function indicating a dark image.
            %            Default is 'dark'.
            % Output : - A matrix of logical indicating if each one of the
            %            HEAD object elements TYPE is consistent with dark.
            % Example: H=HEAD1.basic; Ans=isdark(H)
            % Reliable: 2
            
            if nargin<2
                Val = 'dark';
            end
            
            Ans = isType(H,Val);
            
        end
        
        function Ans=isarc(H,Val)
            % Check if HEAD object TYPE keyword is arc or synonymous
            % Package: @HEAD
            % Description: Check if HEAD object TYPE keyword is arc or
            %              synonymous of arc. The TYPE keyword is selected
            %              using the dictKey function, while the arc value
            %              is selected using the dictValType function.
            % Input  : - An HEAD object.
            %          - The keyword value in the HEAD1.dictValType
            %            function indicating a arc image.
            %            Default is 'arc'.
            % Output : - A matrix of logical indicating if each one of the
            %            HEAD object elements TYPE is consistent with arc.
            % Example: H=HEAD1.basic; Ans=isarc(H)
            % Reliable: 2
            
            if nargin<2
                Val = 'arc';
            end
            
            Ans = isType(H,Val);
            
        end
        
    end  % methods
        
    
    % WCS related methods
    methods
        
    end % end methods
end

            
