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
            
            if (numel(varargin)==0)
                varargin{1} = 1;
                varargin{2} = 1;
            elseif (numel(varargin)==1)
                varargin{2} = 1;
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
    
    % getters
    methods
        function H=get.Header(Head)
            % getter for Header from HEAD object
            % Package: @HEAD
            % Input  : - An HEAD object
            % Output : - If the HEAD object contains a single element, than
            %            this is the content of the Header field (i.e., a 3
            %            column cell array).
            %            Otherwise this is a cell array of cell arrays.
            % By : Eran O. Ofek
            
            
            N = numel(Head);
            if (N==1)
                H = Head.Header;
            else
                for I=1:1:N
                    H{I} = Head(I).Header;
                end
            end
            
        end
        
      
        
    end
    
    
    
    % Static methods: isHEAD, HeaderField
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
            H.Header = {'EXPTIME',1,'';
                        'TYPE','bias','';
                        'FILTER','r',''};
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
            %          - Keyword name.
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
            % Output : - A 3 column cell array with only the lines that
            %            satisfy the search.
            %          - Line indices in the input cell that satisfy the
            %            search criteria.
            % By: Eran O. Ofek
            % Example: [SubCell]=getKeyCell(Cell)
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
            %            cell. Inf for adding at the end, 0 for adding
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
                Cell = [Cell1;Cell2];
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
        
        
     
    end
    
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
    methods
        % copy
        % delKey
        % uniqueKey
        % addKey
        % replaceKeyName
        % replaceKeyVal
        % getKey
        % getVal
        % regexp
        % regexprep
        % lowerKey (previously lower_key)
        % upperKey (previously upper_key)
        % numKey (previously numkey)
        % spacecel
        % cell2struct
        
        % getkey
        % mgetkey
        %
     
        % find_groups
        % iskeyval
        % disp
        % istype
        % julday
        % geodpos
        % coo
        % isarc
        % isbias
        % isdark
        % isflat
        % naxis
        
        
        
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
            %            Inf for adding at the end, 0 for adding
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
            
            N = numel(Head);
            % for each element in object 
            for I=1:1:N
                Out(I).Header = HEAD1.getKeyCell(Head(I).Header,KeyName,varargin{:});
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
        
        function S=cell2struct(H,Unique)
            % Convert the Header fields in HEAD object to structure
            % Package: @HEAD
            % Description: Convert the Header fields in HEAD object to
            %              structure array in which the field name is the
            %              keyword name and the field value is the keyword
            %              value. Optionally run uniqueKey on the HEAD
            %              object.
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
                S(I) = cell2struct(H(I).Header(:,2),H(I).Header(:,1));
            end
            
        end
        
    end
    
    
end

            
