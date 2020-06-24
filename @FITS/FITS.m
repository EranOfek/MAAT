%--------------------------------------------------------------------------
% FITS class                                                         class
% Description: A static class for astronomical FITS images.
%              Include functions to read, write and manipulate FITS images,
%              FITS tables, and FITS headers.
%              Type "FITS." followed by <tab> to see the full list of
%              functions.
%              Full manual is available in manual_FITS.pdf
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef FITS 
    
    % FITS HDU utilities
    methods (Static)
        
        % get number of HDUs
        function Nhdu=num_hdu(Image)
            % get number of HDUs
            % Package: @FITS
            % Description: Get number of HDU in FITS file
            % Input  : - FITS image name.
            % Output : - Number of HDUs
            % Example: N=FITS.num_hdu('Image.fits');
            
            import matlab.io.*
            Fptr = fits.openFile(Image);
            Nhdu = fits.getNumHDUs(Fptr);
            
            fits.closeFile(Fptr);
        end
    end % methods
    
    % FITS header
    methods (Static)
        
        % get header from a sinle FITS file
        function [HeadCell,Nhdu]=read_header(Image,HDUnum)
            % Read a single header into a three column cell array
            % Package: @FITS (Static)
            % Input  : - A single FITS file name.
            %          - HDU number. Default is 1.
            % Output : - A three column cell array with the entire header.
            %          - Number of HDUs in file.
            % Example:
            % [HeadCell,Nhdu]=FITS.read_header('PTF_201211224002_i_p_scie_t093619_u014676207_f02_p100037_c02.fits')
            
            if nargin<2
                HDUnum = 1;
            end
            
            import matlab.io.*
            

            KeyPos = 9;
            ComPos = 32;

            Fptr = fits.openFile(Image);
            Nhdu = fits.getNumHDUs(Fptr);
            if (Nhdu>=HDUnum)
                Htype = fits.movAbsHDU(Fptr,HDUnum);

                Nkey = fits.getHdrSpace(Fptr);
                HeadCell = cell(Nkey,3);
                for Ikey = 1:1:Nkey
                   Card     = fits.readRecord(Fptr,Ikey);
                   LenCard = length(Card);
                   if (LenCard>=9)

                       if (strcmpi(Card(KeyPos),'='))
                           HeadCell{Ikey,1}  = Util.string.spacedel(Card(1:KeyPos-1));
                           % update comment position due to over flow
                           Islash = strfind(Card(ComPos:end),'/');
                           if (isempty(Islash))
                               UpdatedComPos = ComPos;
                           else
                               UpdatedComPos = ComPos + Islash(1)-1;
                           end
                           Value = Card(KeyPos+1:min(LenCard,UpdatedComPos-1));
                           PosAp = strfind(Value,'''');

                           if (isempty(PosAp))
                               if contains('TF',upper(strtrim(Value)))
                                   % a boolean
                                   Value=upper(strtrim(Value))=='T';
                               else
                                   % possible number
                                   Value = str2double(Value);
                               end
                           else
                               if (length(PosAp)>=2)
                                   % a string
                                   Value = strtrim(Value(PosAp(1)+1:PosAp(2)-1));
                               else
                                   Value = Card(PosAp(1)+10:end);
                               end
                           end

                           HeadCell{Ikey,2}  = Value; %Card(KeyPos+1:min(LenCard,ComPos-1));
                           if (LenCard>UpdatedComPos)
                               HeadCell{Ikey,3}  = Card(UpdatedComPos+1:end);    
                           else
                               HeadCell{Ikey,3}  = '';
                           end

                       end
                   end
                   
                   % look for history and comment keywords
                   if (strcmpi(Card(1:7),'HISTORY'))
                       HeadCell{Ikey,1} = 'HISTORY';
                       HeadCell{Ikey,2} = Card(KeyPos:end);
                       HeadCell{Ikey,3} = '';
                   end
                   if (strcmpi(Card(1:7),'COMMENT'))
                       HeadCell{Ikey,1} = 'COMMENT';
                       HeadCell{Ikey,2} = Card(KeyPos:end);
                       HeadCell{Ikey,3} = '';
                   end
                   
                end
            end
            
            fits.closeFile(Fptr);
        
        end
        
        % get header into an headCl object
        function [H,Nhdu]=header2headCl(ImageList,HDUnum)
            % Read a list of FITS headers into an headCl object
            % Package: @FITS (Static)
            % Description: Read a specific Header Data Unit (HDU) in a FITS file
            %              into a an headCl object.
            % Input  : - FITS file name.
            %          - Index of HDU. Default is 1.
            % Output : - An Head object containing the header information:
            %            {Keyword, Value, comment}.
            %          - Number of HDU identified in FITS image.
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Jul 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [H,Nhdu]=FITS.header2headCl('*.fits');
            % Reliable: 2
            
            if nargin<2
                HDUnum = 1;
            end
            
            [~,List] = Util.files.create_list(ImageList,NaN);
            Nlist    = numel(List);
            H        = headCl(Nlist,1);
            Nhdu     = zeros(Nlist,1);

            
            for Ilist=1:1:Nlist
                %Ilist
                Image = List{Ilist};

                [HeadCell,Nhdu(Ilist)]=FITS.read_header(Image,HDUnum);
                
                H(Ilist).Header = HeadCell;
                
            end % end for Ilist...

            
            
        end
        
        % get header from FITS file into an HEAD object
        function [Head,Nhdu]=get_head(ImageList,HDUnum,PopWCS)
            % Read FITS header into a cell array
            % Package: @FITS
            % Description: Read a specific Header Data Unit (HDU) in a FITS file
            %              into a cell array of {Keyword, Value, comment}.
            % Input  : - FITS file name.
            %          - Index of HDU. Default is 1.
            %          - Populate WCS. Default is true.
            % Output : - An Head object containing the header information:
            %            {Keyword, Value, comment}.
            %          - Number of HDU identified in FITS image.
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Jul 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [HeadCell,Nhdu]=FITS.get_head(Image,3);
            % Reliable: 2
            %--------------------------------------------------------------------------
            HeaderField = HEAD.HeaderField;

            Def.HDUnum = 1;
            Def.PopWCS = true;
            if (nargin==1)
                HDUnum   = Def.HDUnum;
                PopWCS   = Def.PopWCS;
            elseif (nargin==2)
                PopWCS   = Def.PopWCS;
            elseif (nargin==3)
                % do nothing
            else
                error('Illegal number of input arguments: get_head(ImageList,[HDUnum,PopWCS])');
            end

            if (isnan(HDUnum))
                HDUnum = 1;
            end

            [~,List] = Util.files.create_list(ImageList,NaN);
            Nlist    = numel(List);
            Head     = HEAD(Nlist,1);
            Nhdu     = zeros(Nlist,1);

            for Ilist=1:1:Nlist
                %Ilist
                Image = List{Ilist};

                [HeadCell,Nhdu(Ilist)]=FITS.read_header(Image,HDUnum);
                
                Head(Ilist).(HeaderField) = HeadCell;
                
            end % end for Ilist...

            
            if (PopWCS)
                Head = populate_wcs(Head);
            end
        end  % end function
        
        % get keyword values from a single FITS file
        function [KeysVal,KeysComment,Struct]=get_keys(Image,Keys,HDUnum,Str)
            % Get single FITS header keywords value
            % Package: @FITS
            % Description: Get the values of specific keywords from a single
            %              FITS file header. Use only for existing keywords.
            % Input  : - FITS image name.
            %          - Cell array of keys to retrieve from the image header.
            %          - HDU number. Default is 1.
            %            If NaN, then set to 1.
            %          - Check if the keyword value is char and try to convert to
            %            a number {false|true}. Default is false.
            % Output : - Cell array of keyword values.
            %          - Cell array of keyword comments.
            %          - Structure containing the keyword names (as fields)
            %            and their values.
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Jul 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [KeysVal,KeysComment,Struct]=FITS.get_keys('A.fits',{'NAXIS1','NAXIS2'});
            % Reliable: 2 
            %--------------------------------------------------------------------------

            if (nargin==2)
                HDUnum = 1;
                Str    = false;
            elseif (nargin==3)
                Str    = false;
            else
                % do nothing
            end

            if (isnan(HDUnum))
                HDUnum = 1;
            end

            import matlab.io.*
            Fptr = fits.openFile(Image);
            N = fits.getNumHDUs(Fptr);
            if (HDUnum>N)
                fits.closeFile(Fptr);
                error('requested HDUnum does not exist');
            end
            fits.movAbsHDU(Fptr,HDUnum);

            if (ischar(Keys))
                Keys = {Keys};
            end

            Nkey = numel(Keys);


            KeysVal     = cell(size(Keys));
            KeysComment = cell(size(Keys));
            for Ikey=1:1:Nkey
                [KeysVal{Ikey},KeysComment{Ikey}] = fits.readKey(Fptr,Keys{Ikey});
                if (ischar(KeysVal{Ikey}) && Str)
                    Tmp = str2double(KeysVal{Ikey});
                    if (isnan(Tmp))
                        % do nothing - keep as a string
                    else
                        KeysVal{Ikey} = Tmp;
                    end
                end

            end
            fits.closeFile(Fptr);

            if (nargout>2)
               Struct = cell2struct(KeysVal,Keys,2);
            end
        end
        
        % get keywords from multiple FITS files
        function [KeysVal,KeysComment,Struct,List]=mget_keys(Images,Keys,HDUnum,Str)
            % Get header keywords value from multiple FITS
            % Package: @FITS
            % Description: Get the values of specific keywords from a list of
            %              FITS files header.
            % Input  : - List of FITS image names. See Util.files.create_list.m for options.
            %          - Cell array of keys to retrieve from the image header.
            %          - HDU number. Default is 1.
            %          - Check if the keyword value is char and try to convert to
            %            a number {false|true}. Default is false.
            % Output : - Cell array (per image) of cell array of keyword values.
            %          - Cell array (per image) of cell array of keyword comments.
            %          - Structure array (element per image) containing the keyword
            %            names (as fields) and their values.
            %          - Cell array containing the list of images.
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Jul 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example:
            % [KeysVal,KeysComment,Struct,List]=FITS.mget_keys('PTF_201202*.fits',{'NAXIS1','NAXIS2'});
            % Reliable: 2
            %--------------------------------------------------------------------------


            if (nargin==2)
                HDUnum = 1;
                Str    = false;
            elseif (nargin==3)
                Str    = false;
            else
                % do nothing
            end

            [~,List] = Util.files.create_list(Images,NaN);
            Nim = numel(List);
            KeysVal     = cell(size(List));
            KeysComment = cell(size(List));
            for Iim=1:1:Nim
               [KeysVal{Iim},KeysComment{Iim},Struct(Iim)]=FITS.get_keys(List{Iim},Keys,HDUnum,Str);
            end

        end
        
        % delete keywirds
        function delete_keys(ImageName,Keywords)
            % Delete a lits of keywords from a list of FITS headers
            % Package: @FITS
            % Description: Delete a list of header keywords from a list of
            %              FITS images.
            % Input  : - List of FITS image names to read. See Util.files.create_list.m for
            %            options.
            %          - Cell array of keyword names to delete.
            % Output : null
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Jun 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: FITS.delete_keys('A.fits',{'PTFPID','OBJECT'})
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (~iscell(Keywords))
                Keywords = {Keywords};
            end
            Nkey = numel(Keywords);

            [~,List] = Util.files.create_list(ImageName,NaN);
            Nim = numel(List);

            import matlab.io.*
            for Iim=1:1:Nim
                Fptr = fits.openFile(List{Iim},'readwrite');
                for Ikey=1:1:Nkey
                    fits.deleteKey(Fptr,Keywords{Ikey});
                end
                fits.closeFile(Fptr);
            end
        end
        
        % write keywords
        function write_keys(ImageName,KeyCell)
            % Insert or update FITS header keywords
            % Package: @FITS
            % Description: Insert new, or update existing FITS header keywords in
            %              a list of FITS images.
            % Input  : - List of FITS image names to edit. See Util.files.create_list.m for
            %            options.
            %          - A cell array of two or three columns of key/cal/comments to
            %            add to FITS header.
            % Output : null
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Jun 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: FITS.write_keys('A.fits',{'try','A','comm';'try2',6,'what'});
            % Reliable: 2
            %--------------------------------------------------------------------------


            Nkey = size(KeyCell,1);

            [~,List] = Util.files.create_list(ImageName,NaN);
            Nim = numel(List);

            import matlab.io.*
            for Iim=1:1:Nim
                Fptr = fits.openFile(List{Iim},'readwrite');
                for Ikey=1:1:Nkey
                    %KeyCell{Ikey,:}
                    if isempty(KeyCell{Ikey,3})
                        KeyCell{Ikey,3} = ' ';
                    end
                    if (strcmp(KeyCell{Ikey,1},'SIMPLE'))
                        KeyCell{Ikey,2} = true;
                    end
                    if (~isempty(KeyCell{Ikey,1}))
                        fits.writeKey(Fptr,KeyCell{Ikey,:});
                    end
                end    
                fits.closeFile(Fptr);
            end

        end
        
        % show header
        function HeadStr=head(File,varargin)
            % Print FITS header to screen
            % Package: @FITS
            % Description: Read FITS file header, convert it to a string, and
            %              print it to screen. Supress empty header lines.
            % Input  : - A fits image file name, or a structure array containing
            %            an .Header field.
            %          * Additional arguments that will be passed to fitsinfo.m.
            % Output : - A string containing the header, with carridge return
            %            between lines.
            % See also: disphead.m
            % Tested : Matlab R2011b
            %     By : Eran O. Ofek                    Sep 2013
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: fitshead PTF201001131475_2_o_14235_00.w.fits
            %          HeadStr=FITS.head('PTF201001131475_2_o_14235_00.w.fits');
            % Reliable: 2
            %--------------------------------------------------------------------------

            HeaderField = 'Header';

            if (isstruct(File))
               HeaderCell = File.(HeaderField);
            else
               H = fitsinfo(File);
               HeaderCell = H.PrimaryData.Keywords;
            end
            N = size(HeaderCell,1);
            HeadStr = '';
            for I=1:1:N
               if (isempty(HeaderCell{I,1}) && isempty(HeaderCell{I,2}) && isempty(HeaderCell{I,3}))
                  % empty line - supress
               else
                  if (ischar(HeaderCell{I,2}))
                      HeadStr = sprintf('%s%8s%20s%s\n',HeadStr,HeaderCell{I,1},HeaderCell{I,2},HeaderCell{I,3});
                  else
                      HeadStr = sprintf('%s%8s%20.10f%s\n',HeadStr,HeaderCell{I,1},HeaderCell{I,2},HeaderCell{I,3});
                  end
               end
            end
        end

        
        
    end  % end methods
    
    % Cell header
    % treating header in cell
    % More recomended to use the HEAD methods
    methods (Static)
        
        function [NewHeader]=cellhead_addcomment(Header,Type,Comments)
            %--------------------------------------------------------------------------
            % FITS.cellhead_addcomment function                            class/@FITS
            % Description: Add comment to a FITS header in a cell array
            % Input  : - An Nx3 cell array of FITS header.
            %          - Either 'COMMENT' or 'HISTORY'.
            %            If empty then 'COMMENT'.
            %          - Cell array of comments or history comments.
            % Output : - New header.
            % Tested : Matlab R2013a
            %     By : Eran O. Ofek                    Mar 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [NewHeader]=FITS.cellhead_addcomment([],[],{'This is a comment'});
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (isempty(Type))
                Type = 'COMMENT';
            end

            if (~iscell(Comments))
                Comments = {Comments};
            end

            Ncom = numel(Comments);
            AddHead = cell(Ncom,3);
            for Icom=1:1:Ncom
                AddHead{Icom,1} = Type;
                AddHead{Icom,2} = '';
                AddHead{Icom,3} = Comments{Icom};
            end

            NewHeader = FITS.cellhead_addkey(Header,AddHead);
        end
        
        function [NewCellHead]=cellhead_addkey(CellHead,varargin)
            %--------------------------------------------------------------------------
            % FITS.cellhead_addkey function                                class/@FITS
            % Description: A utility program to add new keywords, values and
            %              comments to a cell array containing A FITS header
            %              information.
            %              The FITS header cell array contains an arbitrary number
            %              of rows and 3 columns, where the columns are:
            %              {keyword_name, keyword_val, comment}.
            %              The comment column is optional.
            % Input  : - Cell array containing the FITS header information.
            %          * An header cell array to add at the end of the existing
            %            array (e.g., FITS.cellhead_addkey(CellHead,AddHead);).
            %            Alternatively, vecor of positions in the existing header
            %            array in which to add the new keys
            %            (e.g., FITS.cellhead_addkey(CellHead,NewPos,AddHead);).
            %            Alternatively, triplets of :...,key,value,comment,...
            %            to add at the end of the existing header
            %            (e.g., FITS.cellhead_addkey(CellHead,'NAXIS',1024,'number of axes');)
            %            Alternatively, quadraplets of :...,pos,key,value,comment,...
            %            to add at a given position in the existing header specified by
            %            the integer pos.
            %            (e.g., FITS.cellhead_addkey(CellHead,Pos,'NAXIS',1024,'number of axes');)
            % Output : - The new header cell array.
            % Tested : Matlab 7.10
            %     By : Eran O. Ofek                    Jun 2010
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Reliable: 2
            %--------------------------------------------------------------------------
            if (isempty(CellHead))
               CellHead = cell(0,3);
            end

            Narg = length(varargin);
            if (Narg==0)
               % Do nothing
               AddHead = [];
            elseif (Narg==1)
               AddHead = varargin{1};
               VecPos  = zeros(size(AddHead,1),1)+Inf;
            elseif (Narg==2)
               VecPos     = varargin{1};
               AddHead = varargin{2};
            else
               if (Narg./3==floor(Narg./3) && ischar(varargin{1})==1)
                  % assume triplets: ...,key,value,comment,...
                  Counter = 0;
                  AddHead = cell(Narg./3,3);
                  VecPos  = zeros(Narg./3,1) + Inf;
                  for Iarg=1:3:Narg
                     Counter = Counter + 1; 
                     [AddHead{Counter,1:3}] = deal(varargin{Iarg:1:Iarg+2});
                  end
               elseif (Narg./4==floor(Narg./4) && isnumeric(varargin{1})==1)
                  % assume quadraplets: ...,pos,key,value,comment,...
                  Counter = 0;
                  AddHead = cell(Narg./4,3);
                  VecPos  = zeros(Narg./4,1);
                  for Iarg=1:4:Narg
                     Counter = Counter + 1; 
                     VecPos(Counter) = varargin{Iarg};
                     [AddHead{Counter,1:3}] = deal(varargin{Iarg+1:1:Iarg+3});
                  end
               else
                  error('Unknown format of additional parameters');
               end
            end

            [~,Ncr] = size(CellHead);
            [~,Nar] = size(AddHead);

            if (Ncr==Nar)
               % do nothing
            else
               if (Ncr==2 && Nar==3)
                  CellHead{1,3} = [];
               elseif (Ncr==3 && Nar==2)
                  AddHead{1,3} = [];
               else
                  error('Illegal number of columns');
               end
            end

            % sort AddHead by VecPos
            NewCellHead = CellHead;

            [SortedVecPos,SortedInd] = sort(VecPos);
            SortedAddHead = AddHead(SortedInd,:);
            Nadd = length(SortedVecPos);
            for Iadd=1:1:Nadd
               NewCellHead = Util.array.insert_ind(NewCellHead,SortedVecPos(Iadd)+(Iadd-1),SortedAddHead(Iadd,:));
            end

            NewCellHead = FITS.cellhead_fix(NewCellHead);
        end
        
        function [NewCellHead]=cellhead_delkey(CellHead,Keys)
            %--------------------------------------------------------------------------
            % FITS.cellhead_delkey function                                class/@FITS
            % Description: A utility program to delete a keywords, values and
            %              comments from a cell array containing A FITS header
            %              information.
            %              The FITS header cell array contains an arbitrary number
            %              of rows and 3 columns, where the columns are:
            %              {keyword_name, keyword_val, comment}.
            %              The comment column is optional.
            % Input  : - Cell array containing the FITS header information.
            %          - A string containing a keyword name to delete,
            %            or a cell of strings containing keyword names to delete.
            % Output : - The new header cell array.
            % Tested : Matlab 7.10
            %     By : Eran O. Ofek                    Jun 2010
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (ischar(Keys))
               Keys = {Keys};
            elseif (iscell(Keys))
               % do nothing
            else
               error('Unknown Keys DataType');
            end

            NewCellHead = CellHead;
            Nkeys = length(Keys);
            for Ikeys=1:1:Nkeys
               [~,Lines]=FITS.cellhead_getkey(NewCellHead,Keys{Ikeys});
               NewCellHead = Util.array.delete_ind(NewCellHead,Lines);
               Nl = length(Lines);
               for I=1:1:Nl-1
                  [~,Lines]=FITS.cellhead_getkey(NewCellHead,Keys{Ikeys});
                  NewCellHead = Util.array.delete_ind(NewCellHead,Lines);
               end
            end
        end
        
        function NewHeader=cellhead_fix(Header)
            %--------------------------------------------------------------------------
            % FITS.cellhead_fix function                                   class/@FITS
            % Description: Given an Nx3 cell array of FITS header. Remove blank lines
            %              and make sure the END keyword is at the end of the header.
            % Input  : - An Nx3 cell array of FITS header.
            % Output : - A fixed header.
            % Tested : Matlab R2013a
            %     By : Eran O. Ofek                    Mar 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: NewHeader=FITS.cellhead_fix(Header);
            % Reliable: 2
            %--------------------------------------------------------------------------

            FlagEmpty = strcmp(Header(:,1),'') & strcmpi(Header(:,2),'') & strcmpi(Header(:,3),'');
            NewHeader = Header(~FlagEmpty,:);

            % remove END
            FlagEnd   = strcmp(NewHeader(:,1),'END');
            NewHeader = NewHeader(~FlagEnd,:);

            % add END
            NewHeader = [NewHeader; {'END','',''}];
        end
        
        function [NewCellHead,Lines]=cellhead_getkey(CellHead,Keys,NotExist,Multiple)
            %--------------------------------------------------------------------------
            % FITS.cellhead_getkey function                                class/@FITS
            % Description: A utility program to get a specific keywords, values and
            %              comments from a cell array containing A FITS header
            %              information.
            %              The FITS header cell array contains an arbitrary number
            %              of rows and 3 columns, where the columns are:
            %              {keyword_name, keyword_val, comment}.
            %              The comment column is optional.
            % Input  : - Cell array containing the FITS header information.
            %          - A string containing a keyword to look in the header
            %            or a cell array of keywords (case insensitive).
            %          - A parameter to control the behaviour when a specific keyword
            %            is not found.
            %            'Ignore' - ignore missing parameters and in that case the
            %                       length of NewCellHead will be shorter than the
            %                       length of (CellHead). Default.
            %            'NaN'    - Replace missing keywords by NaNs.
            %                       In that case the Lines vector and the NewCellHead
            %                       will contain NaNs
            %          - A flag indicating what to do if there are multiple apperances
            %            of a specific keyword. Options are: {'all','first','last'}.
            %            Default is 'last'.
            % Output : - The header cell array with only the requested lines.
            %          - Indices of requested lines in header.
            % Tested : Matlab 7.10
            %     By : Eran O. Ofek                    Jun 2010
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [NewCellHead,Lines]=FITS.cellhead_getkey(CellHead,Keys,NotExist);
            % Reliable: 2
            %--------------------------------------------------------------------------

            Def.NotExist = 'Ignore';
            Def.Multiple = 'last';
            if (nargin==2)
               NotExist = Def.NotExist;
               Multiple = Def.Multiple;
            elseif (nargin==3)
                Multiple = Def.Multiple;
            elseif (nargin==4)
               % do nothing
            else
               error('Illegal number of input arguments');
            end

            if (ischar(Keys))
               Keys = {Keys};
            end

            Nkeys = length(Keys);
            Lines = zeros(0,1);
            for Ikeys=1:1:Nkeys
                switch lower(Multiple)
                    case 'all'
                        Ifound = find(strcmpi(CellHead(:,1),Keys{Ikeys})==1);
                    otherwise
                        Ifound = find(strcmpi(CellHead(:,1),Keys{Ikeys})==1,1,Multiple);
                end


               Lines  = [Lines; Ifound];

               if (isempty(Ifound))
                  switch lower(NotExist)
                   case 'nan'
                      Lines = [Lines; NaN];
                   otherwise
                      % do nothing
                  end
               end
            end

            if (sum(isnan(Lines))>0)
               Inan = find(isnan(Lines));
               Lines(Inan) = 1;
               NewCellHead = CellHead(Lines,:);
               for In=1:1:length(Inan)
                  NewCellHead(Inan(In),:) = {NaN NaN NaN};
               end
            else
               NewCellHead = CellHead(Lines,:);
            end
        end
        
        function SubHead=cellhead_search(Head,String,varargin)
            %--------------------------------------------------------------------------
            % FITS.cellhead_search function                                class/@FITS
            % Description: Search for substring in FITS header stored as a cell
            %              array of 3 columns. Return and display all instances of
            %              the substring.
            % Input  : - Cell array containing 3 columns FITS header
            %            {key, val, comment}.
            %          - Sub string to search, or a regular expression patern to match.
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'CaseSens' - case sensative {false|true}. Default is false.
            %            'Col'      - Which column in the header to search:
            %                         1 - keywords column.
            %                         3 - comment column.
            %                         You can specificy one or many options
            %                         (e.g., [1 3] or 1). Default is [1 3].
            %            'Verbose'  - Show search result {true|false}. Default is true.
            % Output : - Cell array containing the lines which contain the requested
            %            sub string.
            % License: GNU general public license version 3
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Apr 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: SubHead=FITS.cellhead_search(Sim(1).Header,'FWHM');
            % Reliable: 3
            %--------------------------------------------------------------------------

            DefV.CaseSens          = false;
            DefV.Col               = [1 3];
            DefV.Verbose           = true;
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            if (InPar.CaseSens)
                FindFun = @regexp;
            else
                FindFun = @regexpi;
            end

            Ncol = numel(InPar.Col);
            AllLines = [];
            for Icol=1:1:Ncol
                ColInd = InPar.Col(Icol);
                AllLines = [AllLines; find(~Util.cell.isempty_cell(FindFun(Head(:,ColInd),String,'match')))];
            end
            AllLines = unique(AllLines);
            SubHead = Head(AllLines,:);

            if (InPar.Verbose)
                Nl = numel(AllLines);
                if (Nl>0)
                    fprintf('\n');
                    fprintf('Keyword      Value                   Comment\n');
                    fprintf('----------   --------------------    ----------------------------------------\n');
                else
                    fprintf('\n Substring not found \n');
                end
                for Il=1:1:Nl
                    if (isnumeric(SubHead{Il,2}))
                        fprintf('%-10s   %20.9f    %-40s\n',SubHead{Il,1},SubHead{Il,2},SubHead{Il,3});
                    else
                        fprintf('%-10s   %-20s    %-40s\n',SubHead{Il,1},SubHead{Il,2},SubHead{Il,3});
                    end
                end
            end

        end
        
        function CellHead=cellhead_update(CellHead,varargin)
            %--------------------------------------------------------------------------
            % FITS.cellhead_update function                                class/@FITS
            % Description: Update keywords in fits header cell.
            % Input  : - Cell array containing the FITS header information.
            %          * Arbitrary number of arguments:
            %            ...,{key,val,comment},{key,val},...
            %            ...,key,val,comment,key,val,comment,...
            % Output : - New cell arry of fits header.
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Jan 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: CellHead=FITS.cellhead_update(CellHead,'EXPTIME',60,'');
            %          CellHead=FITS.cellhead_update(CellHead,{'EXPTIME',60,'Exposure time'});
            % Reliable: 2
            %--------------------------------------------------------------------------


            Narg = numel(varargin);
            if (Narg>0)
                if (iscell(varargin{1}))
                    for Iarg=1:1:Narg
                        CellHead = FITS.cellhead_update1(CellHead,varargin{Iarg}{:});
                    end
                else

                    for Iarg=1:3:Narg-2
                        CellHead = FITS.cellhead_update1(CellHead,varargin{Iarg:Iarg+2});
                    end
                end
            end
            CellHead = FITS.cellhead_fix(CellHead);
        end
       
        function CellHead=cellhead_update1(CellHead,Key,Val,Comment)
            %--------------------------------------------------------------------------
            % FITS.cellhead_update1 function                               class/@FITS
            % Description: Update a single keyword in fits header cell.
            % Input  : - Cell array containing the FITS header information.
            %          * Arbitrary number of arguments:
            %            ...,{key,val,comment},{key,val},...
            %            ...,key,val,comment,key,val,comment,...
            % Output : - New cell arry of fits header.
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Jan 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: CellHead=FITS.cellhead_update1(CellHead,'EXPTIME',60,'');
            % Reliable: 2
            %--------------------------------------------------------------------------

            Nkey = size(CellHead,1);
            if (nargin==3)
                Comment = '';
            end
            Ikey = find(strcmp(CellHead(:,1),Key));
            if (isempty(Ikey))
                % key doesn't exist - add
                % add key val at the end of the 
                CellHead(end+1,:) = {Key,Val,Comment};
            else
                % update
                CellHead(Ikey(1),:) = {Key,Val,Comment};
                if (numel(Ikey)>1)
                    CellHead = CellHead(setdiff((1:1:Nkey),Ikey(2:end)),:);
                end
            end
        end
        
    end % methods
    
    % FITS WCS
    methods (Static)
       
        % get WCS/SIP info
        function SIP_Struct=get_sip(HeadCell)
            % Get the SIP(WCS) keywords information from a FITS image.
            % Package: @FITS
            % Description: Get the SIP(WCS) keywords information from a SIM image
            %              structure array or FITS images.
            % Input  : - An header cell array, a SIM containing a single image,
            %            or a single FITS file name.
            % Output : - A structure array containing the SIP(WCS) header keywords.
            % See also: HEAD/get_wcs
            % Tested : Matlab R2015b
            %     By : Tali Engel                      Dec 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: SIP=FITS.get_sip('A.fits');
            %          S = image2sim('A.fits'); SIP=FITS.get_sip(S);
            %          S = image2sim('A.fits'); SIP=FITS.get_sip(S.Header);
            % Reliable: 2
            %--------------------------------------------------------------------------

            Keywordsip{1} = 'A_ORDER';      % sip polynomial order, axis 1, detector to sky
            Keywordsip{2} = 'B_ORDER';      % sip polynonial order, axis 2, detector to sky  
            Keywordsip{3} = 'A_DMAX';       % maximum correction, axis 1 (bound on the maximum distortion over the array)
            Keywordsip{4} = 'B_DMAX';       % maximum correction, axis 2 (bound on the maximum distortion over the array)
                                            % A_DMAX and B_DMAX can be used to estimate the maximum error that would
                                            % result from not evaluating the distortion polynomial)

            if (ischar(HeadCell))
                % HeadCell is a single FITS file name
                % read fits
                Tmp = image2sim(HeadCell);
                HeadCell = Tmp.Header;
            elseif (SIM.issim(HeadCell))
                if (numel(HeadCell)==1)
                    % HeadCell is a single SIM image
                    HeadCell = HeadCell(1).Header;
                else
                    error('Input SIM contains multiple images');
                end
            elseif (iscell(HeadCell))
                % do nothing
            else
                error('Illegal input HeadCell type');
            end

            OrdHeadCell = FITS.cellhead_getkey(HeadCell,Keywordsip,'NaN','first');
            SIP_Struct = cell2struct(OrdHeadCell(:,2).',Keywordsip,2);

            if (~isnan(SIP_Struct.A_ORDER))
                Aord = SIP_Struct.A_ORDER;
                Bord = SIP_Struct.B_ORDER;

                Asip_len = sum(1:(Aord+1));    
                Bsip_len = sum(1:(Bord+1));    

                % initialization
                Asip_KeywordCell = cell(1,Asip_len);
                Bsip_KeywordCell = cell(1,Bsip_len);
                Ap = zeros(1,Asip_len);
                Aq = zeros(1,Asip_len);
                %Acoeff = zeros(1,Asip_len);
                Bp = zeros(1,Bsip_len);
                Bq = zeros(1,Bsip_len);
                %Bcoeff = zeros(1,Bsip_len);
                % Possible P,Q values (since Pa+Qa<=A_ORDER and Pb+Qb<=B_ORDER)
                % => KeyWords cells for the sip polynomials
                I = 0;                      %initialization
                for P = 0:Aord
                    for Q = 0:(Aord-P)
                        I = I+1;
                        Asip_KeywordCell{I} = sprintf('A_%d_%d',P,Q);
                        Ap(I) = P;
                        Aq(I) = Q;
                    end
                end
                A_SIP_HeadCell=FITS.cellhead_getkey(HeadCell,Asip_KeywordCell,'NaN','first');
                Acoeff = [A_SIP_HeadCell{:,2}];
                Acoeff(isnan(Acoeff)) = 0.0;

                I = 0;                      %initialization
                for P = 0:Bord
                    for Q = 0:(Bord-P)
                        I = I+1;
                        Bsip_KeywordCell{I} = sprintf('B_%d_%d',P,Q);     
                        Bp(I) = P;
                        Bq(I) = Q;
                    end
                end
                B_SIP_HeadCell=FITS.cellhead_getkey(HeadCell,Bsip_KeywordCell,'NaN','first');
                Bcoeff = [B_SIP_HeadCell{:,2}];
                Bcoeff(isnan(Bcoeff)) = 0.0;

                %SIP_Struct.Aorder = Aord;
                SIP_Struct.Ap = Ap;
                SIP_Struct.Aq = Aq;
                SIP_Struct.Acoeff = Acoeff;
                %SIP_Struct.Admax = Amax;
                %SIP_Struct.Border = Bord;
                SIP_Struct.Bp = Bp;
                SIP_Struct.Bq = Bq;
                SIP_Struct.Bcoeff = Bcoeff;
                %SIP_Struct.Bdmax = Bmax;
                %SIP_Struct = cell2struct(SIP_HeadCell(:,2).',SIPKeywordCell,2);

                %---
                % coefficients for inverse transformation (sky to detector), if provided in
                % header.
                %---

                Keywordsipinv{1} = 'AP_ORDER';     % inverse sip polynomial order, axis 1, sky to detector
                Keywordsipinv{2} = 'BP_ORDER';     % inverse sip polynomial order, axis 2, sky to detector
                OrdHeadCell_inv = FITS.cellhead_getkey(HeadCell,Keywordsipinv,'NaN','first');
                SIP_Struct_inv = cell2struct(OrdHeadCell_inv(:,2).',Keywordsipinv,2);

                Aord_inv = SIP_Struct_inv.AP_ORDER;
                Bord_inv = SIP_Struct_inv.BP_ORDER;
                if ((~isnan(Aord_inv))&&(~isnan(Bord_inv)))
                    SIP_Struct.inv = SIP_Struct_inv;
                    Asip_len_inv = sum(1:(Aord_inv+1));    
                    Bsip_len_inv = sum(1:(Bord_inv+1));    

                    % initialization
                    Asip_KeywordCell_inv = cell(1,Asip_len_inv);
                    Bsip_KeywordCell_inv = cell(1,Bsip_len_inv);
                    Ap_inv = zeros(1,Asip_len_inv);
                    Aq_inv = zeros(1,Asip_len_inv);
                    %Acoeff = zeros(1,Asip_len);
                    Bp_inv = zeros(1,Bsip_len_inv);
                    Bq_inv = zeros(1,Bsip_len_inv);
                    %Bcoeff = zeros(1,Bsip_len);
                    % Possible P,Q values (since Painv+Qainv<=AP_ORDER and Pbinv+Qbinv<=BP_ORDER)
                    % => KeyWords cells for the (inverse) sip polynomials
                    I = 0;                      %initialization
                    for P = 0:Aord_inv
                        for Q = 0:(Aord_inv-P)
                            I = I+1;
                            Asip_KeywordCell_inv{I} = sprintf('AP_%d_%d',P,Q);
                            Ap_inv(I) = P;
                            Aq_inv(I) = Q;
                        end
                    end
                    Ainv_SIP_HeadCell=FITS.cellhead_getkey(HeadCell,Asip_KeywordCell_inv,'NaN','first');
                    Acoeff_inv = [Ainv_SIP_HeadCell{:,2}];
                    Acoeff_inv(isnan(Acoeff_inv)) = 0.0;

                    I = 0;                      %initialization
                    for P = 0:Bord_inv
                        for Q = 0:(Bord_inv-P)
                            I = I+1;
                            Bsip_KeywordCell_inv{I} = sprintf('BP_%d_%d',P,Q);     
                            Bp_inv(I) = P;
                            Bq_inv(I) = Q;
                        end
                    end
                    Binv_SIP_HeadCell=FITS.cellhead_getkey(HeadCell,Bsip_KeywordCell_inv,'NaN','first');
                    Bcoeff_inv = [Binv_SIP_HeadCell{:,2}];
                    Bcoeff_inv(isnan(Bcoeff_inv)) = 0.0;

                    %SIP_Struct.Aorder = Aord;
                    SIP_Struct.inv.Ap = Ap_inv;
                    SIP_Struct.inv.Aq = Aq_inv;
                    SIP_Struct.inv.Acoeff = Acoeff_inv;
                    %SIP_Struct.Admax = Amax;
                    %SIP_Struct.Border = Bord;
                    SIP_Struct.inv.Bp = Bp_inv;
                    SIP_Struct.inv.Bq = Bq_inv;
                    SIP_Struct.inv.Bcoeff = Bcoeff_inv;
                    %SIP_Struct.Bdmax = Bmax;
                    %SIP_Struct = cell2struct(SIP_HeadCell(:,2).',SIPKeywordCell,2);

                end

            else
                SIP_Struct = [];
            end

        end
        
        % get WCS/PV info
        function TPV_Struct=get_tpv(HeadCell)
            % Get the PV(WCS) keywords information from a FITS image.
            % Package: @FITS
            % Description: Get the PV(WCS) keywords information from a SIM image
            %              structure array or FITS images.
            % Input  : - An header cell array, a SIM containing a single image,
            %            or a single FITS file name.
            % Output : - A structure array containing the PV(WCS) header keywords.
            % See also: HEAD/get_wcs
            % Tested : Matlab R2015b
            %     By : Tali Engel                      Dec 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: PV=FITS.get_tpv('A.fits');
            %          S = image2sim('A.fits'); PV=FITS.get_tpv(S);
            %          S = image2sim('A.fits'); PV=FITS.get_tpv(S.Header);
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (ischar(HeadCell))
                % HeadCell is a single FITS file name
                % read fits
                Tmp = image2sim(HeadCell);
                HeadCell = Tmp.Header;
            elseif (SIM.issim(HeadCell))
                if (numel(HeadCell)==1)
                    % HeadCell is a single SIM image
                    HeadCell = HeadCell(1).Header;
                else
                    error('Input SIM contains multiple images');
                end
            elseif (iscell(HeadCell))
                % do nothing
            else
                error('Illegal input HeadCell type');
            end

            % initialization
            PV1_KeywordCell = cell(1,40);
            PV2_KeywordCell = cell(1,40);
            PV11 = 0;
            PV21 = 0;

            % the distortion polynomials are with up to 40 coefficients (0-39)
            % (corresponding to a max polynomial order of 7)
            for Nc = 0:39
                PV1_KeywordCell{Nc+1} = sprintf('PV1_%d',Nc);
                PV2_KeywordCell{Nc+1} = sprintf('PV2_%d',Nc);
            end

            PV1_HeadCell=FITS.cellhead_getkey(HeadCell,PV1_KeywordCell,'NaN','first');
            PV2_HeadCell=FITS.cellhead_getkey(HeadCell,PV2_KeywordCell,'NaN','first');

            PV1 = [PV1_HeadCell{:,2}];
            PV2 = [PV2_HeadCell{:,2}];

            % missing PV keywords default to 0 except for PV1_1 and PV2_1 which default
            % to 1
            if (any(find(isnan(PV1))==2))
                PV11 = 1;
            end
            if (any(find(isnan(PV2))==2))
                PV21 = 1;
            end
            PV1(isnan(PV1)) = 0.0;
            PV2(isnan(PV2)) = 0.0;
            if (PV11)
                PV1(2) = 1;
            end
            if (PV21)
                PV2(2) = 1;
            end

            TPV_Struct.pv1 = PV1;
            TPV_Struct.pv2 = PV2;

        end
        
    end % methods
    
    
    % FITS image read
    methods (Static)
        
        % fitsread (matlab builtin/slow)
        function Data=fitsread(varargin)
            % Description: Call the fitsread.m function
            % Input  : - Single file name
            %          * Additional parameters to pass to fitsread.m
            %            See fitsread.m for options.
            % Output : - Data in FITS file.
            % Example: Data=FITS.fitsread('File.fits');
            
            Data=fitsread(varargin{:});
          
        end
        
        % read/ read section (faster than fitsread)
        function Im=read(ImageName,HDUnum,StartPix,EndPix)
            % Read a rectangular region of interest from a single FITS image.
            % Package: @FITS
            % Description: Read a rectangular region of interest from a single
            %              FITS image.
            % Input  : - String containing FITS image name.
            %          - HDU number. Default is 1. If empty use default.
            %          - Start pixels position [x, y].
            %            Alternatively, this can be [xmin, xmax, ymin, ymax].
            %            In this case the third input argument should not provided.
            %            If empty then read the entire image. Default is
            %            empty.
            %          - End pixels position [x, y].
            % Output : - FITS image sub section.
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Jun 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Im=FITS.read('Image.fits',[],[11 11],[70 70]);
            %          Im=FITS.read('Image.fits',[],[11 70 11 70]);
            % Reliable: 2
            %--------------------------------------------------------------------------

            Def.HDUnum   = 1;
            Def.StartPix = [];
            Def.EndPix   = [];
            
            if (nargin==1)
                HDUnum   = Def.HDUnum;
                StartPix = [];
                EndPix   = [];
            elseif (nargin==2)
                StartPix = [];
                EndPix   = [];
            elseif (nargin==3)
                if (isempty(StartPix))
                    % do nothing
                else
                    EndPix   = StartPix([2,4]);
                    StartPix = StartPix([1,3]);
                end
            else
                % do nothing
            end
            if (isempty(HDUnum))
                HDUnum = Def.HDUnum;
            end
            
            import matlab.io.*
            Fptr = fits.openFile(ImageName);
            fits.movAbsHDU(Fptr,HDUnum);
            
            if (isempty(StartPix))
                % read full image
                Im = fits.readImg(Fptr);
            else
                % read image section
                EndPix   = fliplr(EndPix);
                StartPix = fliplr(StartPix);
            
                Im = fits.readImg(Fptr,StartPix,EndPix);
            end
            fits.closeFile(Fptr);

        end

        % read to SIM
        function Sim=read2sim(Images,varargin)
            % Description: Read FITS images into SIM object.
            %              Can read N-dimensional images.
            %              Can also read multi extension files.
            % Input  : - List of images to read. See Util.files.create_list.m for
            %            details.
            %          * Arbitrary number of pairs of ...,key,val,...
            %            input arguments. Available keywords are:
            %            'HDUnum' - HDU to read. If multiple numbers are
            %                       given then will attempt to read multiple
            %                       HDU each one into different SIM
            %                       element.
            %                       If empty, then will attemp to read all
            %                       extensions.
            %                       Default is 1.
            %            'CCDSEC' - A four column matrix of CCDSEC
            %                       [xmin xmax ymin ymax] to read.
            %                       Either line per image or a single line.
            %            'Sim'    - An existing SIM into to write the FITS
            %                       images. If empty, then create a new
            %                       SIM object. Default is empty.
            %            'ExecField'- SIM field into which to write the
            %                       FITS images. Default is 'Im'.
            %            'ReadHead'- Read header into SIM. Default is true.
            %            'HDUnum' - Index of HDU. Default is 1.
            %            'PopWCS' - Populate WCS. Default is true.
            % Output: - A SIM object with the FITS images.
            % Example: S=FITS.read2sim('Image*.fits');
            %          S=FITS.read2sim('Image6[15-28].fits');
            %          S=FITS.read2sim('@list');
            %          S=FITS.read2sim('Image,fits','CCDSEC',[1 10 1 100]);
            % Reliable: 2
            
            HeaderField = HEAD.HeaderField;
            FileField   = SIM.FileNameField;
            WCSField    = 'WCS';
            
            DefV.HDUnum               = 1;
            DefV.CCDSEC               = [];  % section to read
            DefV.Sim                  = [];  % read into existing SIM
            DefV.ExecField            = SIM.ImageField;   % read into field
            DefV.ReadHead             = true;
            DefV.HDUnum               = 1;
            DefV.PopWCS               = true;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            [~,ListIm] = Util.files.create_list(Images,NaN);
            Nim = numel(ListIm);
            
            if (isempty(InPar.Sim))
                % allocate SIM
                Sim = SIM(Nim,1);
            else
                % Use existing SIM
                Sim = InPar.Sim;
            end
            
            % number of lines in CCDSEC
            Nccdsec = size(InPar.CCDSEC,1);
            
            % for each image
            Isim = 0; % image (including extensions) index
            for Iim=1:1:Nim
                % get number of HDU in FITS image
                if (isempty(InPar.HDUnum))
                    Nhdu   = FITS.num_hdu(ListIm{Iim});
                    HDUnum = (1:1:Nhdu);
                else
                    HDUnum = InPar.HDUnum;
                    Nhdu   = numel(HDUnum);
                end
                % for each HDU
                for Ihdu=1:1:Nhdu
                    Isim = Isim + 1;
                    % Read image to SIM
                    if (isempty(InPar.CCDSEC))
                        Sim(Isim).(InPar.ExecField) = FITS.read(ListIm{Iim},HDUnum(Ihdu));
                    else
                        Sim(Isim).(InPar.ExecField) = FITS.read(ListIm{Iim},HDUnum(Ihdu),InPar.CCDSEC(min(Iim,Nccdsec),:));
                    end
                    Sim(Isim).(FileField) = ListIm{Iim};

                    % read header
                    if (InPar.ReadHead)
                        H = FITS.get_head(ListIm{Iim},HDUnum(Ihdu),InPar.PopWCS);
                        Sim(Isim).(HeaderField) = H.(HeaderField);
                        Sim(Isim).(WCSField)    = H.(WCSField);
                    end
                end
            end
            
        end
        
        % read to cube (without header)
        function Cube=read2cube(Images,varargin)
            % Description: Read FITS images into a cube.
            %              Assume all the images have the same size.
            %              Can also read multi extension files.
            % Input  : - List of images to read. See Util.files.create_list.m for
            %            details.
            %          * Arbitrary number of pairs of ...,key,val,...
            %            input arguments. Available keywords are:
            %            'HDUnum' - HDU to read. If multiple numbers are
            %                       given then will attempt to read multiple
            %                       HDU each one into different SIM
            %                       element.
            %                       If empty, then will attemp to read all
            %                       extensions.
            %                       Default is 1.
            %            'CCDSEC' - A four column matrix of CCDSEC
            %                       [xmin xmax ymin ymax] to read.
            %                       Either line per image or a single line.
            % Output: - Cube of images in which the first dimension is the
            %           image index.
            % Example: S=FITS.read2cube('Image*.fits');
            %          S=FITS.read2cube('Image6[15-28].fits');
            %          S=FITS.read2cube('@list');
            %          S=FITS.read2cube('Image,fits','CCDSEC',[1 10 1 100]);
            % Reliable: 2
            
            DefV.HDUnum               = 1;
            DefV.CCDSEC               = [];  % section to read
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            [~,ListIm] = Util.files.create_list(Images,NaN);
            Nim = numel(ListIm);
            
            % number of lines in CCDSEC
            Nccdsec = size(InPar.CCDSEC,1);
            
            % for each image
            Isim = 0; % index of image (including extensions)
            for Iim=1:1:Nim
                % get number of HDU in FITS image
                if (isempty(InPar.HDUnum))
                    Nhdu   = FITS.num_hdu(ListIm{Iim});
                    HDUnum = (1:1:Nhdu);
                else
                    HDUnum = InPar.HDUnum;
                    Nhdu   = numel(HDUnum);
                end
                % for each HDU
                for Ihdu=1:1:Nhdu
                    Isim = Isim + 1;
                
                    % Read image to SIM
                    if (isempty(InPar.CCDSEC))
                        Image = FITS.read(ListIm{Iim},HDUnum(Ihdu));
                    else
                        Image = FITS.read(ListIm{Iim},HDUnum(Ihdu),InPar.CCDSEC(min(Iim,Nccdsec),:));
                    end

                    if (Isim==1)
                        Size = size(Image);
                        Cube = zeros(Nim,Size(1),Size(2));
                    end
                    if ~all(size(Image)==Size)
                        error('Image %s has different size compared with image %s',ListIm{Iim},ListIm{1});
                    end
                    
                    Cube(Isim,:,:) = Image;
                end
            end       
                
            
        end
        
    end % Static
    
    
    % FITS image write
    methods (Static)
        
        % fitswrite (matlab builtin)
        function fitswrite(varargin)
            % Description: Call the matlab builtin fitswrite function
            % Input  : - Data 
            %          - Filename
            %          * Additional parameters to pass to fitswrite.
            % Output : null
            % Example: FITS.fitswrite(rand(100,100),'Image.fits');
            % Reliable: 2
            
            fitswrite(varargin{:});
        end
        
        % write FITS file
        function Flag=write(Image,FileName,varargin)
            % Description: Write or append an image into FITS file.
            %              The image may have N-dimensions.
            %              Append will write multi extension FITS image.
            % Input  : - Array to save as FITS image.
            %          - FITS file name to save.
            %          * Arbitrary number of ...,key,val,... pairs.
            %            Following keywords are available:
            %            'Header' - Cell array of {key,val,comment} header
            %                       or an HEAD object to write into the
            %                       FITS file.
            %            'DataType' - Data type - default is 'single'
            %                       precision.
            %            'Append' - Append image as a multi extension to an
            %                       existing FITS file. Default is false.
            %            'OverWrite'- Overwrite an existing image. Default
            %                       is false.
            %            'WriteTime'- Add creation time to image header.
            %                       Default is false.
            % Example: Flag=FITS.write(rand(100,100),'Try.fits');
            %          Flag=FITS.write(rand(10,10,3),'Try.fits');
            %
            
            
            HeaderField = HEAD.HeaderField;
            
            DefV.Header             = [];
            DefV.DataType           = -32;
            DefV.Append             = false;  % true for multi extension
            DefV.OverWrite          = false;
            DefV.WriteTime          = false;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            % Set FITS DataType
            switch InPar.DataType
                 case {'int8',8}
                    DataType = 'uint8';
                 case {'int16',16}
                     % apparently uint16 is not supported?! in 2017b
                    DataType = 'int16';
                 case {'int32',32}
                     % apparently uint16 is not supported?! in 2017b
                    DataType = 'int32';
                 case {'int64',64}
                    DataType = 'int64';
                 case {'single','float32',-32}
                    DataType = 'single';
                 case {'double','float64',-64}
                    DataType = 'double';
                 otherwise
                    error('Unknown DataType option');
            end
            
            % Overwrite existing FITS file
            if (InPar.OverWrite)
                % delete existing FileName if exist
                if (exist(FileName,'file')~=0)
                    delete(FileName);
                end
            end
            
            % Prepare header
            if (HEAD.ishead(InPar.Header))
                % already HEAD object
                Header = InPar.Header;
            else
                % convert to HEAD object
                Header = HEAD;
                Header.(HeaderField) = InPar.Header;
            end
            %--- Set the FITS "mandatory" keywords ---
            %--- add BSCALE and BZERO ---
            % check if BZERO and BSCALE are already in HeaderInfo 
            Header = replace_key(Header,'BZERO',   single(0),  'offset data range to that of unsigned short',...
                                        'BSCALE',  single(1),  'default scaling factor');
                               
            %--- Write creation date to header ---
            if (InPar.WriteTime)
                Time = celestial.time.get_atime([],0,0); % Na'ama, 20180516
                %Header = replace_key(Header,'CRDATE',  Time.ISO,'Creation date of FITS file',...
                %                            'COMMENT', '',      'File Created by MATLAB FITS.write.m written by E. Ofek');
                Header = replace_key(Header,'CRDATE',  Time.ISO,'Creation date of FITS file',...
                                            'COMMENT', 'File Created by MATLAB FITS.write.m written by E. Ofek', ''); % Na'ama, 20180518
                                        
%=======
%                                            'COMMENT', 'File Created by MATLAB FITS.write.m written by E. Ofek', '');
%>>>>>>> 35a2a05383ff2610cc108265c1d7c311b103e2ac
            end
            Nline = size(Header.(HeaderField),1);

            import matlab.io.*
            if (InPar.Append)
                % append to existing FITS file
                Fptr = fits.openFile(FileName,'READWRITE');
            else
                % Create new FITS file
                Fptr = fits.createFile(FileName);
            end
            
            

            % create Image
            fits.createImg(Fptr,DataType,size(Image));
            % write Image
            fits.writeImg(Fptr,Image); %,Fpixels);
            
            % write Header
            for Inl=1:1:Nline
                if (~isempty(Header.(HeaderField){Inl,1}))
                    switch lower(Header.(HeaderField){Inl,1})
                        case 'comment'
                            fits.writeComment(Fptr,Header.(HeaderField){Inl,2});
                        case 'history'
                            fits.writeHistory(Fptr,Header.(HeaderField){Inl,2});
                        case {'extname', 'xtension'} % Na'ama, 20180905
                            % do nothing
                        case {'simple','bitpix','naxis','naxis1','naxis2','naxis3','naxis4'}
                            % do nothing - these keywords are written by
                            % the FITS creator
                         case {'bscale', 'bzero', 'datamax', 'datamin', 'epoch', 'equinox'} % Na'ama, 2018-06-06
                            % convert to floating point, required by matlab.io.fits.writeKey
                            if (ischar(Header.(HeaderField){Inl,2}))
                                Header.(HeaderField){Inl,2} = double(eval(Header.(HeaderField){Inl,2}));
                            end
                            if (isempty(Header.(HeaderField){Inl,3}))
                                Header.(HeaderField){Inl,3} = ' ';
                            end
                            fits.writeKey(Fptr,Header.(HeaderField){Inl,1},...
                                               Header.(HeaderField){Inl,2},...
                                               Header.(HeaderField){Inl,3});   
                        case 'end'
                            % do nothing
                        otherwise
                            if (isnan(Header.(HeaderField){Inl,2}))
                                Header.(HeaderField){Inl,2} = ' ';
                            end
                            if (isempty(Header.(HeaderField){Inl,3}))
                                Header.(HeaderField){Inl,3} = ' ';
                            end
%                             Header.(HeaderField){Inl,:}
%                             string(Header.(HeaderField){Inl,3})
%                             regexprep(Header.(HeaderField){Inl,3},'\W&\S','')
                            % dela with non-standard keywords
                            if isempty(Header.(HeaderField){Inl,2})
                                Header.(HeaderField){Inl,2} = ' ';
                            end
                            fits.writeKey(Fptr,Header.(HeaderField){Inl,1},...
                                               Header.(HeaderField){Inl,2},...
                                               Header.(HeaderField){Inl,3});
                    end
                end
            end
            % Close FITS file
            fits.closeFile(Fptr);
            Flag = sign(Fptr);


            
        end
        
        % write_old (previously fitswrite_my)
        function [Flag,HeaderInfo]=write_old(Image,FileName,HeaderInfo,DataType,varargin)
            % Write a simple 2D FITS image. OBSOLETE: Use FITS.write instead.
            % Package: @FITS
            % Description: Write a simple 2D FITS image.
            %              OBSOLETE: Use FITS.write instead.
            % Input  : - A 2D matrix to write as FITS file.
            %          - String containing the image file name to write.
            %          - Cell array containing header information to write to image.
            %            The keywords: SIMPLE, BITPIX, NAXIS, NAXIS1, NAXIS2, EXTEND
            %            will be re-written to header by default.
            %            The keywords BSCALE and BZERO wil be written to header if
            %            not specified in the header information cell array.
            %            Alternatively this could be a character array (Nx80)
            %            containing the header (no changes will be applied).
            %            If not given, or if empty matrix (i.e., []) than write a
            %            minimal header.
            %          - DataType in which to write the image, supported options are:
            %            'int8',8            
            %            'int16',16
            %            'int32',32
            %            'int64',64
            %            'single','float32',-32    (default)
            %            'double','float64',-64            
            %          * Arbitrary number of pairs of input arguments: 
            %            ...,keyword,value,... - possible keys are:
            %            'IdentifyInt'  - attempt to identify integers in header
            %                             automaticaly and print them as integer.
            %                             {'y' | 'n'}, default is 'y'.
            %            'ResetDefKey'  - Reset default keys {'y' | 'n'},
            %                             default is 'y'.
            %                             Default keys are:
            %                             {'SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2',
            %                              'EXTEND','BSCALE',BZERO'}.
            %            'OverWrite'    - Over write existing image {'y' | 'n'},
            %                             default is 'y'.
            % Output : - Flag indicating if image was written to disk (1) or not (0).
            %          - Actual header information written to file.
            % See also: fitswrite_my.m
            % Bugs   : Don't write standard FITS file
            % Tested : Matlab 7.10
            %     By : Eran O. Ofek                      June 2010
            %    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
            % Examples: [Flag,HeaderInfo]=FITS.write_old(rand(2048,1024),'Example.fits');
            % Reliable: 2
            %--------------------------------------------------------------------------


            Def.HeaderInfo = [];
            Def.DataType   = -32;
            if (nargin==2)
               HeaderInfo = Def.HeaderInfo;
               DataType   = Def.DataType;
            elseif (nargin==3)
               DataType   = Def.DataType;
            end

            % set default for additional keywords:
            DefV.IdentifyInt = 'y';
            DefV.ResetDefKey = 'y';
            DefV.OverWrite   = 'y';

            InPar = InArg.populate_keyval(DefV,varargin,mfilename);


            Flag = 1;

            switch DataType
             case {'int8',8}
                DataType = 8;
             case {'int16',16}
                DataType = 16;
             case {'int32',32}
                DataType = 32;
             case {'int64',64}
                DataType = 64;
             case {'single','float32',-32}
                DataType = -32;
             case {'double','float64',-64}
                DataType = -64;
             otherwise
                error('Unknown DataType option');
            end



            if (ischar(HeaderInfo)==0)

            %--- Set the FITS "mandatory" keywords ---
            ImSize = size(Image);
            switch lower(InPar.ResetDefKey)
               case 'n'
                  % do nothing
               case 'y'
                  if (isempty(HeaderInfo))
                     % do nothing
                  else
                     % delete "default" keys:
                     HeaderInfo = FITS.cellhead_delkey(HeaderInfo,{'SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND','BSCALE','BZERO'});
                  end
                  HeaderInfo = FITS.cellhead_addkey(HeaderInfo,...
                                                    0,'SIMPLE',true(1,1),'file does conform to FITS standard',...
                                                    0,'BITPIX',int32(DataType),'number of bits per data pixel',...
                                                    0,'NAXIS' ,int32(length(ImSize)),'number of data axes',...
                                                    0,'NAXIS1',int32(ImSize(2)),'length of data axis 1',...
                                                    0,'NAXIS2',int32(ImSize(1)),'length of data axis 2',...
                                                    0,'EXTEND',false(1,1),'FITS dataset may contain extensions',...
                                                    0,'BZERO' ,single(0),'offset data range to that of unsigned short',...
                                                    0,'BSCALE',single(1),'default scaling factor');
               otherwise
                  error('Unknown ResetDefKey option');
            end

            % check if last keyword is END - delete if needed.
            if (strcmp(HeaderInfo{end,1},'END'))
                HeaderInfo = HeaderInfo(1:end-1,:);
            end


            %--- Write creation date to header ---
            Time = celestial.time.get_atime([],0,0);
            HeaderInfo = FITS.cellhead_addkey(HeaderInfo,...
                                              'CRDATE',Time.ISO,'Creation date of FITS file',...
                                              'COMMENT','','File Created by MATLAB fitswrite.m written by Eran Ofek');

            %--- Convert default keywords to int32 ---
            KeysInt32 = {'BITPIX','NAXIS','NAXIS1','NAXIS2','BZERO','BSCALE'};
            Ni32     = length(KeysInt32);
            for Ii32=1:1:Ni32
               I = find(strcmp(HeaderInfo(:,1),KeysInt32{Ii32})==1);
               if (~isempty(I))
                  HeaderInfo{I,2} = int32(HeaderInfo{I,2});
               end
            end


            %--- Convert default keywords to logical ---
            KeysLogi = {'SIMPLE','EXTEND'};
            Nlogi     = length(KeysLogi);
            for Ilogi=1:1:Nlogi
               I = find(strcmp(HeaderInfo(:,1),KeysLogi{Ilogi})==1);
               if (~isempty(I))
                  if (~islogical(HeaderInfo{I,2}))
                     switch HeaderInfo{I,2}
                      case 'F'
                         HeaderInfo{I,2} = false(1,1); 
                      case 'T'
                         HeaderInfo{I,2} = true(1,1);
                      otherwise
                         error('Keyword type is not logical');
                     end
                  end
               end
            end


            % check if last keyword is END - if not add.
            if (~strcmp(HeaderInfo{end,1},'END'))
                HeaderInfo{end+1,1} = 'END';
                HeaderInfo{end,2} = '';
                HeaderInfo{end,3} = '';
            end

            %--- Prepare string of header information to write to header ---
            HeaderBlock = '';
            [Nline,Nr] = size(HeaderInfo);

            Counter = 0;
            for Iline=1:1:Nline
               if (~isempty(HeaderInfo{Iline,1}) && strcmpi(HeaderInfo{Iline,1},'END')==0)
                  %--- write keyword name ---
                  HeaderLine = sprintf('%-8s',upper(HeaderInfo{Iline,1}));
                  switch upper(HeaderInfo{Iline,1})
                   case {'COMMENT','HISTORY'}
                      % do not write "=" sign
                      HeaderLine = sprintf('%s',HeaderLine);
                   otherwise
                      % write "=" sign
                      HeaderLine = sprintf('%s= ',HeaderLine);
                  end
                  %--- write keyword value ---
                  switch upper(HeaderInfo{Iline,1})
                   case {'COMMENT','HISTORY'}
                      % do not write value
                   otherwise
                      if (ischar(HeaderInfo{Iline,2}))
                         HeaderLine = sprintf('%s''%s''',HeaderLine,HeaderInfo{Iline,2});
                         Nblanks    = 20-(length(HeaderInfo{Iline,2})+2);

                         if (Nblanks<0)
                            Nblanks = 0;
                         end
                         HeaderLine = sprintf('%s%s',HeaderLine,blanks(Nblanks));

                      elseif (islogical(HeaderInfo{Iline,2}))
                          switch HeaderInfo{Iline,2}
                          case 1
                         HeaderLine = sprintf('%s%20s',HeaderLine,'T');
                          case 0
                         HeaderLine = sprintf('%s%20s',HeaderLine,'F');                 
                          otherwise
                         error('DataType is not logical');
                         end
                      elseif (isinteger(HeaderInfo{Iline,2}))
                         HeaderLine = sprintf('%s%20i',HeaderLine,HeaderInfo{Iline,2});
                      elseif (isfloat(HeaderInfo{Iline,2}))

                         switch InPar.IdentifyInt
                          case 'y'
                             % attempt to identify integers automatically
                             if (fix(HeaderInfo{Iline,2})==HeaderInfo{Iline,2})
                                % integer identified - print as integer
                                HeaderLine = sprintf('%s%20i',HeaderLine,HeaderInfo{Iline,2});
                             else
                                % float or double
                                HeaderLine = sprintf('%s%20.8f',HeaderLine,HeaderInfo{Iline,2});
                             end
                          case 'n'
                               % float or double
                               HeaderLine = sprintf('%s%20.8f',HeaderLine,HeaderInfo{Iline,2});
                          otherwise
                                 error('Unknown IdentifyInt option');
                         end
                      else
                         error('unknown Format in header information');
                      end
                  end
                  %--- write comment to header ---
                  if (Nr>2)
                     switch upper(HeaderInfo{Iline,1})
                      case {'COMMENT','HISTORY'}
                         % do not write "/"
                         HeaderLine = sprintf('%s%-s',HeaderLine,HeaderInfo{Iline,3});
                      otherwise
                         if (isempty(HeaderInfo{Iline,3}))
                            % do nothing - do not add comment
                         else
                            HeaderLine = sprintf('%s /%-s',HeaderLine,HeaderInfo{Iline,3});
                         end
                      end
                   end
                   % cut line if longer than 80 characters or paded with spaces
                   if (length(HeaderLine)>80)
                      HeaderLine = HeaderLine(1:80);
                   end
                   HeaderLine = sprintf('%-80s',HeaderLine);

                   %%HeaderBlock = sprintf('%s%s',HeaderBlock,HeaderLine);
                   Counter = Counter + 1;
                   HeaderBlock(Counter,:) = sprintf('%s',HeaderLine);
                end   
            end
            %--- Add End keyword to end of header ---
            %%HeaderBlock = sprintf('%s%-80s',HeaderBlock,'END');
            Counter = Counter + 1;
            HeaderBlock(Counter,:) = sprintf('%-80s','END');

            else
              %--- HeaderInfo is already provided as char array ---
              % assume contains also the 'END' keyword
              HeaderBlock = HeaderInfo;
            end


            %--- pad by spaces up to the next 2880 character block ---
            PadBlock = sprintf('%s',blanks(2880 - mod(numel(HeaderBlock),2880)));


            %--- pad by spaces up to the next 2880 character block ---
            %HeaderBlock = sprintf('%s%s',HeaderBlock,blanks(2880 - 80 -
            %mod(length(HeaderBlock),2880)));
            %--- Add End keyword to end of header ---
            %HeaderBlock = sprintf('%s%-80s',HeaderBlock,'END');


            %--- Open file for writing in rigth byte order in all platforms.
            switch lower(InPar.OverWrite)
               case 'n'
                  if (exist(FileName,'file')==0)
                     % File does not exist - continue
                  else
                     error('Output file already exist');
                  end
               otherwise
                  % do nothing
            end

            Fid = fopen(FileName,'w','b');  % machineformat: ieee-be
            if (Fid==-1)
               fprintf('Error while attempting to open FITS file for writing\n');
               Flag = 0;
            else   
               %--- Write header to file ---

               %fwrite(Fid,HeaderBlock,'char');

               for Ic=1:1:size(HeaderBlock,1)
                 fprintf(Fid,'%-80s',HeaderBlock(Ic,:));
               end
               fprintf(Fid,'%s',PadBlock);

               %--- Write image data ---
               switch DataType
                case {'int8',8}
                   fwrite(Fid,Image.','uchar');
                case {'int16',16}
                   fwrite(Fid,Image.','int16');
                case {'int32',32}
                   fwrite(Fid,Image.','int32');
                case {'int64',64}
                   fwrite(Fid,Image.','int64');
                case {'single','float32',-32}
                   fwrite(Fid,Image.','single');
                case {'double','float64',-64}
                   fwrite(Fid,Image.','double');
                otherwise
                   fclose(Fid);
                   error('Unknown DataType option');
               end
               fclose(Fid);
            end
        end % FITS.write_old

        % SIM to FITS
        function OutName=sim2fits(Sim,varargin)
            % Write SIM images as FITS files
            % Package: @FITS
            % Description: Convert SIM object to FITS files
            % Input  : - A Sim object
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'OutName     - Output name, or cell array of output
            %                           names. If empty generate a
            %                           tmp file. Default is empty.
            %            'WriteHead'  - Write header.
            %                           Default is true.
            %            'DataType'   - data type. Default is 'single'.
            %            'OverWrite'- Overwrite an existing image. Default
            %                       is false.
            %            'WriteTime'- Add creation time to image header.
            %                       Default is false.
            %            'Prefix'     - File name prefix.
            %                           Default is ''.
            %            'Suffix'     - File name suffix.
            %                           Default is ''.
            %            'OutDir'     - Directory in which to write the
            %                           files. Default is ''.
            %            'SuffixFITS' - Add 'fits' suffix to file name.
            %                           Default is true.
            %            'ExecField'  - SIM fields to save.
            %                           Default is 'Im'.
            % Output : - Cell array of output names.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Feb 2017
            %    URL : http://weizmann.ac.il/home/eofek/matlab/

            DefV.OutName             = [];
            DefV.WriteHead           = true;
            DefV.DataType            = 'single';
            DefV.OverWrite           = false;
            DefV.WriteTime           = false;
            DefV.Prefix              = '';
            DefV.Suffix              = '';
            DefV.OutDir              = '';
            DefV.SuffixFITS          = true;
            DefV.ExecField           = {SIM.ImageField};
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            if ~isempty(InPar.OutName)
                if ~iscell(InPar.OutName)
                    InPar.OutName = {InPar.OutName};
                end
            end
            
            Nsim = numel(Sim);
            OutName = cell(1,Nsim);
            for Isim=1:1:Nsim
                % for each SIM object
                if isempty(InPar.OutName)
                    OutName{Isim} = tempname;
                else
                    OutName{Isim} = InPar.OutName{Isim};
                end
                
                if (InPar.SuffixFITS)
                    OutName{Isim} = sprintf('%s%s%s%s.fits',InPar.OutDir,InPar.Prefix,OutName{Isim},InPar.Suffix);
                else
                    OutName{Isim} = sprintf('%s%s%s%s',InPar.OutDir,InPar.Prefix,OutName{Isim},InPar.Suffix);
                end
                
                If = 1;
                if (InPar.WriteHead)
                    Header = Sim(Isim).Header;
                else
                    Header = [];
                end
                FITS.write(Sim(Isim).(InPar.ExecField{If}),OutName{Isim},...
                                                           'Header',Header,...
                                                           'DataType',InPar.DataType,...
                                                           'OverWrite',InPar.OverWrite,...
                                                           'WriteTime',InPar.WriteTime);
                
                
            end
        end
        
    end % Static
    
    
    % FITS tables
    methods (Static)

        function [Out,Col]=read_table(TableName,varargin)
            % Read binary or ascii FITS tables.
            % Package: @FITS
            % Description: Read binary or ascii FITS tables.
            % Input  : - List of FITS tables to read. Any input valid to Util.files.create_list.m.
            %            If multiple files are provided then all the files hould be of
            %            the same type (e.g., fits binary table).
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'TableType'- FITS table type {'auto'|'bintable'|'table'}.
            %                         Default is 'auto'.
            %                         'auto' will attempt to read the table type
            %                         from the 'XTENSION' header keyword.
            %            'HDUnum'   - HDU number in which the FITS table header is.
            %                         Default is 2.
            %            'ModColName' - If the program failed because some columns
            %                         have name starting with numbers or invalid signs
            %                         set this to true (default is false).
            %                         This will modify the column names.
            %            'OutTable' - Type of table output:
            %                         'astcat' - AstCat object.
            %                         'astcat_t' - AstCat object in which
            %                                      the data is stored as a table.
            %            'XTENkey'  - Header keyword from which to read the table type.
            %                         Default is 'XTENSION'.
            %            'StartRow' - First row to read from FITS table. Default is [].
            %                         If empty, then read the entire table.
            %            'NRows'    - Number of rows to read from FITS table.
            %                         Default is []. If empty, then read the entire
            %                         table.
            %            'OutClass' - A function to use in order to force all the
            %                         columns to be of the same class (e.g., @single).
            %                         Default is @double. If empty, then will keep
            %                         the original class. This option shoyld be used
            %                         if you want to read the data into a matrix.
            %            'NullVal'  - If the column is of numeric type will attempt
            %                         to replace the FITS null value with this value.
            %                         Default is NaN.
            %            'BreakRepCol'- {true|false}. If true and FITS table columns
            %                         are repeating then will change column information
            %                         according to the matrix column count.
            %                         Default is true.
            % Output : - AstCat object containing the FITS table.
            %          - A structure array of additional columns
            %            information, like format.
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Jan 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [Out,Col]=FITS.read_table('asuWDs.fit','ModColName',true);
            % Reliable: 2
            %--------------------------------------------------------------------------

            HeaderField     = HEAD.HeaderField;
            CatField        = AstCat.CatField;
            ColField        = AstCat.ColField;
            ColCellField    = AstCat.ColCellField;
            ColUnitsField   = AstCat.ColUnitsField';


            DefV.TableType      = 'auto';    % {'auto'|'bintable'|'table'}
            DefV.HDUnum         = 2;
            DefV.ModColName     = false;
            DefV.OutTable       = 'astcat';  % {'astcat'|'astcat_t'|...}
            DefV.XTENkey        = 'XTENSION';
            DefV.StartRow       = [];
            DefV.NRows          = [];
            DefV.OutClass       = @double;
            DefV.NullVal        = NaN;       % [] do nothing
            DefV.BreakRepCol    = true;   

            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            % prep list of fits table names
            [~,ListTableName] = Util.files.create_list(TableName,NaN);
            Nfile = numel(ListTableName);
            
            
            Out = AstCat(Nfile,1);
            Col = Util.struct.struct_def({'Col','Cell','Units','TypeChar','Repeat','Scale','Zero','Nulval','Tdisp','Data'},Nfile,1);
            
            import matlab.io.*
            for Ifile=1:1:Nfile
                
                % get header
                Head     = FITS.get_head(ListTableName{1},InPar.HDUnum);
                % populate header
                Out(Ifile).(HeaderField) = Head.(HeaderField);
                
                % identify table type
                switch lower(InPar.TableType)
                    case 'auto'
                        HeadCell = Head.(HeaderField);
                        Ixten    = find(strcmp(HeadCell(:,1),InPar.XTENkey));
                        if (isempty(Ixten))
                            error('Automatic table identification mode was unable to access FITS table type');
                        end
                        InPar.TableType = HeadCell{Ixten,2};
                    otherwise
                        % do nothing
                end

                % Table type
                switch lower(InPar.TableType)
                    case 'bintable'
                        Fun_getColParms = @fits.getBColParms;
                    case 'table'
                        Fun_getColParms = @fits.getAColParms;
                    otherwise
                        error('Unknown TableType option');
                end

                if (isempty(InPar.StartRow) || isempty(InPar.NRows))
                    CellRowPar = {};
                else
                    CellRowPar = {InPar.StartRow,InPar.NRows};
                end

                
                % for each fits table
                Fptr = fits.openFile(ListTableName{Ifile});
                fits.movAbsHDU(Fptr,InPar.HDUnum);
                Ncol = fits.getNumCols(Fptr);
                %[ttype,tunit,typechar,repeat,scale,zero,nulval,tdisp]= fits.getBColParms(Fptr,2);
                Col(Ifile).Cell     = cell(1,Ncol);
                Col(Ifile).Units    = cell(1,Ncol);
                Col(Ifile).TypeChar = cell(1,Ncol);
                Col(Ifile).Repeat   = cell(1,Ncol);
                Col(Ifile).Scale    = cell(1,Ncol);
                Col(Ifile).Zero     = cell(1,Ncol);
                Col(Ifile).Nulval   = cell(1,Ncol);
                Col(Ifile).Tdisp    = cell(1,Ncol);
                Col(Ifile).Data     = cell(1,Ncol);
                for Icol=1:1:Ncol
                    [Col(Ifile).Cell{Icol},Col(Ifile).Units{Icol},Col(Ifile).TypeChar{Icol},...
                                           Col(Ifile).Repeat{Icol},Col(Ifile).Scale{Icol},Col(Ifile).Zero{Icol},...
                                           Col(Ifile).Nulval{Icol},Col(Ifile).Tdisp{Icol}]= Fun_getColParms(Fptr,Icol);
                                       
                                       
                    [Col(Ifile).Data{Icol}] = fits.readCol(Fptr,Icol,CellRowPar{:});
                    if (~isempty(InPar.OutClass))
                        Col(Ifile).Data{Icol} = InPar.OutClass(Col(Ifile).Data{Icol});
                    end
                    if (~isempty(InPar.NullVal) && ~isempty(Col(Ifile).Nulval{Icol}) && isnumeric(Col(Ifile).Data{Icol}))
                        Col(Ifile).Data{Icol}(Col(Ifile).Data{Icol}==Col(Ifile).Nulval{Icol}) = InPar.NullVal;

                    end
                    % override ColRepeat using the actual data
                    Col(Ifile).Repeat{Icol} = size(Col(Ifile).Data{Icol},2);
                end
                fits.closeFile(Fptr);

                % deal with repeating columns
                if (InPar.BreakRepCol)
                    Nnc  = sum(cell2mat(Col(Ifile).Repeat));
                    NewCol.Cell = cell(Nnc,1);

                    Icol1 = 1;
                    for Icol=1:1:Ncol            
                        IcolN = Icol1 + Col(Ifile).Repeat{Icol} - 1;
                        %Icol1 = Icol1 + Col.Repeat{Icol}; % at the end of the loop
                        for Irep=1:1:Col(Ifile).Repeat{Icol}
                            if (Col(Ifile).Repeat{Icol}>1)
                                NewCol.Cell{Icol1+Irep-1} = sprintf('%s_%d_',Col(Ifile).Cell{Icol},Irep);
                            else
                                NewCol.Cell{Icol1+Irep-1} = Col(Ifile).Cell{Icol};
                            end
                        end
                        [NewCol.Units{Icol1:IcolN}] = deal(Col(Ifile).Units{Icol});
                        [NewCol.TypcChar{Icol1:IcolN}] = deal(Col(Ifile).TypeChar{Icol});
                        [NewCol.Repeat{Icol1:IcolN}] = deal(1);
                        [NewCol.Scale{Icol1:IcolN}] = deal(Col(Ifile).Scale{Icol});
                        [NewCol.Zero{Icol1:IcolN}] = deal(Col(Ifile).Zero{Icol});
                        [NewCol.Tdisp{Icol1:IcolN}] = deal(Col(Ifile).Tdisp{Icol});
                        Icol1 = Icol1 + Col(Ifile).Repeat{Icol}; % at the end of the loop

                    end
                    Col(Ifile).Cell     = NewCol.Cell;
                    Col(Ifile).Units    = NewCol.Units;
                    Col(Ifile).TypeChar = NewCol.TypcChar;
                    Col(Ifile).Repeat   = NewCol.Repeat;
                    Col(Ifile).Scale    = NewCol.Scale;
                    Col(Ifile).Zero     = NewCol.Zero;
                    Col(Ifile).Tdisp    = NewCol.Tdisp;   
                end
                if (InPar.ModColName)
                    % modify column names to legal variable names
                    Col(Ifile).Cell = regexprep(Col(Ifile).Cell,{'-','/','(',')','&','@','#','^','%','*','~'},'');
                    Col(Ifile).Cell = strcat('T',Col(Ifile).Cell);
                end
                Col(Ifile).Col      = cell2struct(num2cell(1:1:length(Col(Ifile).Cell)),Col(Ifile).Cell,2);

                % output
                switch lower(InPar.OutTable)
                    case 'astcat'
                        Out(Ifile).(CatField)     = [Col(Ifile).Data{:}];
                        Out(Ifile).(ColField)     = Col(Ifile).Col;
                        Out(Ifile).(ColCellField) = Col(Ifile).Cell;
                        Out(Ifile).(ColUnitsField)= Col(Ifile).Units;
                        
                    case 'astcat_t'
                        Out(Ifile).(CatField)     = table(Col(Ifile).Data{:});
                        Out(Ifile).(ColField)     = Col(Ifile).Col;
                        Out(Ifile).(ColCellField) = Col(Ifile).Cell;
                        Out(Ifile).(ColUnitsField)= Col(Ifile).Units;
                    
                    otherwise
                        error('Unknown OuTable option');
                end

            end
        
        end 
        
    
    end % methods
    
    % FITS spectrum
    methods (Static)
        function [Table,ColCell]=read_sdss_spec(File,varargin)
            % Read SDSS spectra in FITS file
            % Package: FITS
            % Description: Read SDSS spectra in FITS file into an AstCat or
            %              a matrix.
            % Input  : - A FITS file, a cell array of FITS files, or a
            %            string with wild cards. See Util.files.create_list
            %            for options.
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'OutType' - Options: 'AstCat' | 'cellmat' | 'AstSpec'
            %                        Default is 'AstCat'.
            % Output : - Output spectra with columns
            %            [Wavelength [A], Flux, ivar, and_mask, or_mask,
            %            wdisp, sky, model]
            %          - ColCell
            % Example: [Table,ColCell]=FITS.read_sdss_spec('*.fits');
            
            
            CatField = AstCat.CatField;
            Units    = 1e-17;
            
            DefV.OutType              = 'AstCat';
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            
            [~,List] = Util.files.create_list(File,NaN);
            Nlist    = numel(List);
            
            ColOrder = [2, 1, 3:1:8];
            for Ilist=1:1:Nlist
                Table(Ilist) = FITS.read_table(List{Ilist});
                Table(Ilist).(CatField) = Table(Ilist).(CatField)(:,ColOrder);
                Table(Ilist).(CatField)(:,1) = 10.^(Table(Ilist).(CatField)(:,1));
                Table(Ilist).(CatField)(:,2) = Table(Ilist).(CatField)(:,2).*Units; 
                Table(Ilist).ColCell = Table(Ilist).ColCell(ColOrder);
                Table(Ilist).ColCell{1} = 'wavelength';
                Table(Ilist) = colcell2col(Table(Ilist));
                
            end
            ColCell = Table(1).ColCell;
            
            switch lower(InPar.OutType)
                case 'cellmat'
                    Table = astcat2cell(Table);
                    
                case 'astspec'
                     error('AstSpec option not supported yet');
                otherwise
                    % do nothing
            end
            
        end
    end
    
end % end class
            
 