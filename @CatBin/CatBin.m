% A class for binary catalog file access
% Package: @CatBin
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


classdef CatBin < handle
    properties (SetAccess = public)
        File
        HeaderFile
        FID
        ColTypes
        ColCell
        ColUnits
        SortedByCol
        Nrow
        Ncol
    end
  
    % CatBin File description:
    
    %--- Table file ---
    % 1 Byte - File Type [0 | 1] - 0 index; 1 - tables
    % 2 Byte - Header size
    % N Byte - single header [describes catalog columns]
    % Content -
    % 1 Byte - Number of tables: N_t
    % 2xN_t bytes - 2 bytes indicating the size of each table, per table 
    % N bytes - Table
    % ...
    % ...
    
    %--- index file ---
    % 1 Byte - File Type [0 | 1] - 0 index; 1 - tables1 Byte - File Type [0 | 1] - 0 index; 1 - tables1 Byte - File Type [0 | 1] - 0 index; 1 - tables
    
    %-------------------
    %--- Constractor ---
    %-------------------
    methods
        function CB=CatBin(FileName,varargin)
            % Construct a single element CatBin object
            % Input  : - A File name
            %          * Arbitrary number of ...,keyword,value,...
            %            The following keywords are possible:
            %            'HeaderFile'- Header file name. Default is [].
            %            'ColTypes'  - Coloumn types. Default is 'single'.
            %            'ColCell'   - Cell array of column names.
            %                          Default is [].
            %            'ColUnits'  - Cell array of column units.
            %                          Default is [].
            %            'SortedByCol- Column index by which catalog is sorted.
            %                          Default is [].
            %            'Nrow'      - Number of rows. Default is [].
            %            'Ncol'      - Number of columns. Default is [].
            % Output : - A CatBin object
            
            CB.File = FileName;
            
            DefV.HeaderFile           = [];
            DefV.ColTypes             = 'single';
            DefV.ColCell              = [];
            DefV.ColUnits             = [];
            DefV.SortedByCol          = [];
            DefV.Nrow                 = [];
            DefV.Ncol                 = 0;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);
            
            FN = fieldnames(InPar);
            for Ifn=1:1:numel(FN)
                CB.(FN{Ifn}) = InPar.(FN{Ifn});
            end
            
        end
        
    end
    
    % static methods
    methods (Static)
        function create(Files)
            % Create a multiple element CatBin object
            % Input : - A list of file names. See Util.files.cteate_list
            %           for options.
            % Output : - A CatBin object
            [~,List] = Util.files.create_list(Files,NaN);
            Nl = numel(List);
            for Il=1:1:Nl
                CB(Il) = CatBin(List{Il});
            end
            
        end
    end
    
    methods
        function CB=access_open(CB,varargin)
            % open CatBin file for access
            % Input  : - A CatBin object
            %          * Additional parameters to pass to fopen
            %            Default is 'w'.
            % Output : - A CatBin object
            
            if (numel(varargin)==0)
                varargin{1} = 'w';
            end
            
            Ncb = numel(CB);
            for Icb=1:1:Ncb
                CB(Icb).FID = fopen(CB(Icb).File,varargin{:});
            end
                
            
        end
        
        function CB=access_close(CB)
            % close CatBin file for access
            % Input  : - A CatBin object
            % Output : - A CatBin object
            
            Ncb = numel(CB);
            for Icb=1:1:Ncb
                fclose(CB(Icb).FID);
                CB(Icb).FID = [];
            end
                
            
        end
        
        function CB=write_header(CB,ColCell,ColUnits,SortedByCol)
            % Write header file
            % Package: @CatBin
            % Input  : - A CatBin object
            %          - An optional cell array of column names.
            %            By default it is to read it from CatBin object.
            %          - An optional cell array of column units.
            %            Default is empty.
            %            By Default it  is to read it from CatBin object.
            %          - An integer indicating by which column the catalog
            %            is sorted. If 0 then not sorted. Default is 0.
            %            By Default it  is to read it from CatBin object.
            % Output : - A CatBin object with the ColCell, ColUnits and
            %            SortedByCol fields populated.
            % Example: CB=CatBin('try.11'); CB=write_header(CB,{'RA','Dec'});
            
            
            
            if (nargin>1)
                CB.ColCell = ColCell;
            end
            if (nargin>2)
                CB.ColUnits = ColUnits;
            end
            if (nargin>3)
                CB.SortedByCol = SortedByCol;
            end
            
            % set the default header file name if not exist
            if (isempty(CB.HeaderFile))
                CB.HeaderFile = sprintf('%s.head',CB.File);
            end
            
            Nchar_ColCell = cellfun(@numel,CB.ColCell);
           
            % FORMAT structure
            % Number of columns - uint8
            % Number of chars in each column X Ncol - unit8
            % Number of ColUnits - uint8
            % Number of chars in each col units X Ncol - uint8
            % Number of ColTypes - uint8
            % Number of chars in each col type X Ncol - uint8
            % ColCell
            % ColUnits
            % ColTypes
            % Sorted by column - 0 if not sorted - uint8
            
            FIDh = fopen(CB.HeaderFile,'w');
            % number of columns
            fwrite(FIDh,numel(CB.ColCell),'uint8');
            % Number of chars in each column X Ncol
            fwrite(FIDh,Nchar_ColCell,'uint8');
            % Number of ColUnits 
            fwrite(FIDh,numel(CB.ColUnits),'uint8');
            % Number of chars in each col units X Ncol 
            if (~isempty(CB.ColUnits))
                Nchar_ColUnits = cellfun(@numel,CB.ColUnits);
                fwrite(FIDh,Nchar_ColUnits,'uint8');
            end
            % Number of ColTypes - uint8
            if (ischar(CB.ColTypes))
                NcolType = 1;
                ColTypes = {CB.ColTypes};
            else
                NcolType = numel(CB.ColTypes);
                ColTypes = CB.ColTypes;
            end
            fwrite(FIDh,NcolType,'uint8');
            % Number of chars in each col type X Ncol - uint8
            Nchar_ColTypes = cellfun(@numel,ColTypes);
            fwrite(FIDh,Nchar_ColTypes,'uint8');
            
            % ColCell
            fwrite(FIDh,cell2mat(CB.ColCell),'char');
            % ColUnits
            if (~isempty(CB.ColUnits))
                fwrite(FIDh,cell2mat(CB.ColUnits),'char');
            end
            
            % ColTypes
            fwrite(FIDh,cell2mat(ColTypes),'char');
            
            % Sorted by column
            fwrite(FIDh,CB.SortedByCol,'uint8');
            
            fclose(FIDh);
            
        end
        
        function CB=read_header(CB)
            % Read CatBin header file
            % Package: @CatBin
            
            % set the default header file name if not exist
            if (isempty(CB.HeaderFile))
                CB.HeaderFile = sprintf('%s.head',CB.File);
            end
            
            % FORMAT structure
            % Number of columns - uint8
            % Number of chars in each column X Ncol - unit8
            % Number of ColUnits - uint8
            % Number of chars in each col units X Ncol - uint8
            % Number of ColTypes - uint8
            % Number of chars in each col type X Ncol - uint8
            % ColCell
            % ColUnits
            % ColTypes
            % Sorted by column - 0 if not sorted - uint8
            
            % read the header file
            FIDh = fopen(CB.HeaderFile,'r');
            NcolCell      = fread(FIDh,1,'uint8');
            Nchar_ColCell = fread(FIDh,NcolCell,'uint8');
            NcolUnits     = fread(FIDh,1,'uint8');
            Nchar_ColUnits= fread(FIDh,NcolUnits,'uint8');
            NcolTypes     = fread(FIDh,1,'uint8');
            Nchar_ColTypes= fread(FIDh,NcolTypes,'uint8');
            ColCellChar   = fread(FIDh,sum(Nchar_ColCell),'*char');
            ColUnitsChar  = fread(FIDh,sum(Nchar_ColUnits),'*char');
            ColTypesChar  = fread(FIDh,sum(Nchar_ColTypes),'*char');
            CB.SortedByCol= fread(FIDh,1,'uint8');
            fclose(FIDh);
            
            % re order the ColCell
            CB.ColCell = cell(1,NcolCell);
            I1 = 1;
            for Ic=1:1:NcolCell
                I2 = I1+Nchar_ColCell(Ic)-1;
                CB.ColCell{Ic} = ColCellChar(I1:I2).';
                I1 = I2 + 1;
            end
            
            % re order ColUnits
            CB.ColUnits = cell(1,NcolUnits);
            I1 = 1;
            for Ic=1:1:NcolUnits
                I2 = I1+Nchar_ColUnits(Ic)-1;
                CB.ColUnits{Ic} = ColUnitsChar(I1:I2).';
                I1 = I2 + 1;
            end
            
            % re order ColTypes
            CB.ColTypes = cell(1,NcolTypes);
            I1 = 1;
            for Ic=1:1:NcolTypes
                I2 = I1+Nchar_ColTypes(Ic)-1;
                CB.ColTypes{Ic} = ColTypesChar(I1:I2).';
                I1 = I2 + 1;
            end
            
            
        end
        
        function CB=write(CB,DataCat,Type)
            % Write data into CatBin object file
            
            
            Ncb = numel(CB);
            if (AstCat.isastcat(DataCat))
                Ncat = numel(DataCat);
                if (Ncb~=Ncat)
                    error('CatBin and AstCat object must have the same number of elements');
                end
                for Icat=1:1:Ncat
                    if (isempty(CB(Icat).FID))
                        % open CatBin
                        CB(Icat) = access_open(CB(Icat),'w');
                    end
                    % write the data in row order(!)
                    fwrite(CB(Icat).FID,DataCat.Cat,Type);
                    
                    CB(Icat) = access_vlose(CB(Icat));
                    CB(Icat).ColTypes = Type;
                end
            else
                if (Ncb>1)
                    error('CatBin object must have a single element');
                end
                Icat = 1;
                if (isempty(CB(Icat).FID))
                    % open CatBin
                    CB(Icat) = access_open(CB(Icat),'w');
                end
                % write the data in rwo order(!)
                fwrite(CB(Icat).FID,DataCat,Type)

                CB(Icat) = access_close(CB(Icat));
            end
        end
        
        function [Data,CB]=read(CB,Lines)
            
            if (ischar(CB.ColTypes))
                % all data is kept in the same format
                if (isempty(CB.FID))
                    CB = access_open(CB,'r');
                    CloseAfterRead = true;
                else
                    CloseAfterRead = false;
                end
                
                Size = (Lines(2) - Lines(1)+1).*CB.Ncol;
                Size
                fseek(CB.FID,(Lines(1)-1).*CB.Ncol+1,'bof');
                Data = fread(CB.FID,Size,CB.ColTypes);
                if (CloseAfterRead)
                    CB = access_close(CB);
                end
            end
        end
        
        
        
    end % method
    
    
end
