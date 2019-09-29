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
    methods

        %-------------------------
        %--- Class constructor ---
        %-------------------------
        
        function obj=HEAD1(varargin)
            obj(1).Header = cell(0,3);
            %obj(1).UserData = [];
            %obj = struct_def({'Header','UserData'},varargin{:});
            
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
            
            
            N = numel(H);
            if (N==1)
                H = Head.Header;
            else
                for I=1:1:N
                    H{I} = Head(I).Header;
                end
            end
            
        end
        
    end
    
    
    
    %--- Static methods ---
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
        
    end
    
    % static methods for cell header
    methods (Static)
        
        % get key value from a header in a cell format
        function [Val,Comment,KeyName]=getkeyCell(Cell,Key,Exact,ReturnN)
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
            %          - A flag indicating what to do if multiple keywords
            %            are found. If true return all, if false, only the
            %            first. Default is flase.
            % Output : - A cell array of keyword values.
            %          - A cell array of keyword comment.
            %          - A cell array of keyword name.
            % By: Eran O. Ofek
            
            if (nargin<4)
                ReturnN = 1;
                if (nargin<3)
                    Exact     = true;
                end
            end
            
            if (ischar(Key))
                Key = {Key};
            end
            Nkey = numel(Key);
            % for each keyword name
            Val  = cell(1,Nkey);
            Comment = cell(1,Nkey);
            KeyName = cell(1,Nkey);
            for Ikey=1:1:Nkey
                if (Exact)
                    % use exact name search
                    Flag = strcmp(Cell(:,1),Key{Ikey});
                else
                    % use regular expression
                    Flag = regexp(Cell(:,1),Key{Ikey},'match');
                    Flag = ~Util.cell.isempty_cell(Flag);
                end
                Nfound{Ikey} = 
                FlagI = find(Flag,ReturnN,'first');
                

                Val(:,Ikey)     = Cell(Flag,2);
                Comment(:,Ikey) = Cell(Flag,3);
                KeyName(:,Ikey) = Cell(Flag,1);
%                 if (~ReturnAll && any(Flag))
%                     Val = Val(1);
%                     Comment = Comment(1);
%                     KeyName = KeyName(1);
%                 end
            end
                
        end
    end
    
    % isfield, isstruct
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

    end
    
    % read data from header
    methods
        function getkey(Head,KeyName)
            % Get keyword value from object header
            % Package: HEAD
            % Input  : - HEAD object
            %          - A string containing a single keyword name.
            
            N = numel(Head);
            % for each element in object 
            for I=1:1:N
                Head(I).Header
            end
            
        end
    end
    
    
end

            
