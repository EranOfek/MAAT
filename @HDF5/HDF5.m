% HDF5 class: High level access to hdf5 file
% Package: class/@HDF5
% Description: A class that allows to access the variables in HDF5 files
%              datasets by index of the elements in the array.
%              The class also include some basic tools to save and read
%              HDF5 files.
%              See also: the Util.IO.( package.
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef HDF5
    properties (SetAccess = public)
        File     % File Name
        Var      % Dataset path name
        FID
        DSetID
        Size
        %IsWrite = false
    end
         
    % Basic high level static methods: save, load
    methods (Static)
        % simple save array
        function save(Data,FileName,VarName,Attrib,WriteMode)
            % save a new array into an HDF5 file dataset
            % Package: @HDF5
            % Description: A simple and fast save HDF5 for arrays
            % Input  : - An array to save.
            %          - HDF5 file name to save.
            %          - HDF5 dataset name. Default is '/V'.
            %          - A two column cell array of attribute key,val.
           
            
            
            Def.VarName = '/V';
            Def.Attrib  = {};
           
           
            if (nargin<4)
                Attrib  = Def.Attrib;
                if (nargin<3)
                    VarName = Def.VarName;
                end
            end
            
            Plist  = 'H5P_DEFAULT';
            Dcpl   = 'H5P_DEFAULT';
            
            TypeID=HDF5.set_type(Data);
            
            Size   = size(Data);
            
            % create/open an HDF5 file
            if (exist(FileName,'file')>0)
                % file exist
                FID = H5F.open(FileName,'H5F_ACC_RDWR','H5P_DEFAULT');
            else
                FID = H5F.create(FileName);
            end
            
            % Dataset ID
            
            % create new dataset
            H5_Dims    = fliplr(Size);
            H5_MaxDims = H5_Dims;
            SpaceID = H5S.create_simple(ndims(Data),H5_Dims,H5_MaxDims);
            
            DID = H5D.create(FID,VarName,TypeID,SpaceID,Dcpl);
            %DID = H5D.open(FID,VarName);
            H5D.write(DID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',Plist,Data);
            
            
            % Wtite attributes
            if (~isempty(Attrib))
                Nline = size(Attrib);
                for Iline=1:1:Nline
                    Acpl = H5P.create('H5P_ATTRIBUTE_CREATE');
                    Value = Attrib{Iline,2};
                    AttTypeID=HDF5.set_type(Value);
                    AttSpaceID = H5S.create('H5S_SCALAR');
                    AttID = H5A.create(DID,Attrib{Iline,1},AttTypeID,AttSpaceID,Acpl);
                    H5A.write(AttID,'H5ML_DEFAULT',Attrib{Iline,2})
                    H5A.close(AttID);
                    H5T.close(AttTypeID);
                end
            end
            
            H5D.close(DID);
            H5F.close(FID);
            
        end
        
        function Data=h5read(FileName,VarName,varargin)
            % load HDF5 data using the built in h5read
            % Package: @HDF5
            % Description: Identical to h5read
            % Input  : - HDF5 file name
            %          - Dataset name
            %          - Optional start index
            %          - Optional count index
            % Output : - The data
            
            Data=h5read(FileName,VarName,varargin{:});
            
        end
        
        function Data=load(FileName,VarName,Offset,Block)
            % load an array from a HDF5 
            % Package: @HDF5
            % Input  : - HDF5 file name.
            %          - Dataset name
            %          - [I,J] offset in the array from which to start
            %            uploading. If not given than get the entire array.
            %          - [I,J] block size to upload from array.
            % Output : - An array
            
            FID = H5F.open(FileName); 
            DSetID = H5D.open(FID,VarName);
            if (nargin>2)
                % User requested for specific location in the array
                Plist = 'H5P_DEFAULT';
                %Space = H5D.get_space(DSetID);
                %[~,Dims] = H5S.get_simple_extent_dims(Space);
                
                Dims = double(fliplr(Block));
                MemSpaceID = H5S.create_simple(numel(Dims),Dims,[]);
                FileSpaceID = H5D.get_space(DSetID);
                Offset = double(fliplr(Offset-1));  % index start at 0
                Block = double(fliplr(Block));
                H5S.select_hyperslab(FileSpaceID,'H5S_SELECT_SET',Offset,[],[],Block);
                Data = H5D.read(DSetID,'H5ML_DEFAULT',MemSpaceID,FileSpaceID,Plist);
            else
                Data = H5D.read(DSetID);
            end
            
            % Read attributes
            % FFU
            
            H5D.close(DSetID);
            H5F.close(FID);
            
        end
        
        function writeatt(FileName,DatasetName,AttribCell)
            % Write attributes into HDF5 file
            % Input  : - FileName
            %          - Dataset name
            %          - Cell array of pairs of attributes: key,val,...
            fileattrib(FileName,'+w');
            
            N = numel(AttribCell);
            if ((N./2)~=floor(N./2))
                error('Number of key,val attributes must be even');
            end
            for I=1:2:N-1
                h5writeatt(FileName,DatasetName,AttribCell{I},AttribCell{I+1});
            end
        end
        
        function Data=load_muti_datasets(FileName,VarNames)
            % Load full multiple datasets from a single HDF5 file
            % Package: @HDF5
            % Description: Load the full content of multiple datasets
            %              in a single HDF5 file into the same variable.
            %              Assuming all datasets has the same number of
            %              columns. This is a little bit faster then using
            %              h5read multiple times as the file is opend once.
            % Input  : - HDF5 file name.
            %          - Cell array of dataset names (starts with '/').
            % Output : - Data
            % Example: Data=HDF5.load_muti_datasets('UCAC4_htm_019100.hdf5',{'/htm_019100'});
            % Reliable: 2
            
            
            FID = H5F.open(FileName); 
            Nds = numel(VarNames);
            for Ids=1:1:Nds
                DSetID = H5D.open(FID,VarNames{Ids});
                % User requested for specific location in the array
                if (Ids==1)
                    Data = H5D.read(DSetID);
                else
                    Data = [Data; H5D.read(DSetID)];
                end
                H5D.close(DSetID);
            end
            
            H5F.close(FID);
            
        end
        
        function [Data,Att]=loadh(FileName,VarName,GetAtt)
            % Load HDF5 file and attributes
            % Package: @HDF5
            % Description: load a matrix from HDF5 file. 
            %              If dataset name is not provided than will read all
            %              datasets into a structure. This function doesn't support
            %              groups.
            %              This is becoming faster than matlab (2014a) for matices with
            %              more than ~10^4 elements.
            % Input  : - File name.
            %          - variable name (dataset). If empty or not provided than will
            %            attempt to read all datasets.
            %          - Get attribute. If empty then do not get attributes.
            %            'h' - store attributes in an HEAD object. 
            %            's' - store attributes in a structure array.
            %            Default is empty.
            % Output : - Datasets.
            % License: GNU general public license version 3
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    May 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Data=loadh('R.hd5','R');
            % Reliable: 2
            %--------------------------------------------------------------------------

            Def.VarName = [];
            Def.GetAtt  = [];
            if (nargin==1)
                VarName = Def.VarName;
                GetAtt  = Def.GetAtt;
            elseif (nargin==2)
                GetAtt  = Def.GetAtt;
            else
                % do nothing
            end

            if (isempty(VarName))
                % check for available datasets in HD5 file
                Info = h5info(FileName);
                if (~isempty(Info.Datasets))
                    Nds = numel(Info.Datasets);
                    for Ids=1:1:Nds
                        Data.(Info.Datasets(Ids).Name) = HDF5.loadh(FileName,Info.Datasets(Ids).Name,GetAtt);
                    end   
                    Ind = (1:1:Nds);  % indices of all datasets
                end
            else
                DataSet = sprintf('/%s',VarName);
                Data    = h5read(FileName,DataSet);

                if (GetAtt)
                    Info = h5info(FileName);
                    Ind  = strcmp({Info.Datasets.Name},VarName);
                end
            end

            if (~isempty(GetAtt))
                % Get all attributes for datasets
                Nds = numel(Ind);
                switch lower(GetAtt)
                    case 'h'
                        % attribute will be stored in an HEAD object
                        Att = HEAD(Nds,1);
                        for Ids=1:1:Nds
                            Att(Ids).Header = [{Info.Datasets(Ids).Attributes.Name}.', {Info.Datasets(Ids).Attributes.Value}.'];
                        end          
                    case 's'
                        % attributes will be stored in a structure
                        for Ids=1:1:Nds
                            Att(Ids) = cell2struct({Info.Datasets(Ids).Attributes.Value},{Info.Datasets(Ids).Attributes.Name},2);
                        end
                    otherwise
                        error('Unknown attributes output type');
                end
            end

    
        end
        
        function Var=load_check(FileName,VarName,StoreWS,WS)
            % Load a variable from an HDF5 file or session (if exist)
            % Package: @HDF5
            % Description: Load a variable from an HDF5 file, but first
            %              check if the variable name exist in the matlab
            %              session, and if it exist load it from there.
            %              This is faster for loading large files.
            % Input  : - HDF5 file name.
            %          - Variable name in the HDF5 file.
            %          - Store the variable in the main work space, if doesn't exist,
            %            {true|false}. Default is true.
            %          - Workspace {'base'|'caller'}, Default is 'base'.
            % Output : - Loaded variable.
            % Example: Var=load_check('try.hdf5','HTM');
            % Reliable: 2

            Def.StoreWS = true;
            Def.WS = 'base';
            if (nargin==2)
               StoreWS = Def.StoreWS;
               WS      = Def.WS;
            elseif (nargin==3)
               WS      = Def.WS;
            elseif (nargin==4)
               % do nothing
            end
            
            try
                Var = evalin(WS,sprintf('%s;',VarName));
            catch
                % variable doesn't exist in workspace
                Var = HDF5.load(FileName,VarName);
               if (StoreWS)
                    assignin(WS,VarName,Var);
               end
            end
            
%             if (evalin(WS,sprintf('exist(''%s'')',VarName))==1)
%                % variable exist in workspace
%                Var = evalin(WS,sprintf('%s;',VarName));
%             else
%                Var = HDF5.load(FileName,VarName);
%                if (StoreWS)
%                     assignin(WS,VarName,Var);
%                end
%             end

            
            
        end
        
        function AttTypeID=set_type(Value)
            % Set HDF5 data type using the H5T.copy command
            % Package: @HDF5
            % Input  : - Data value
            % Output : - TypeID
            
            
            switch class(Value)
                case 'double'
                    AttTypeID = H5T.copy('H5T_NATIVE_DOUBLE');
                case 'single'
                    AttTypeID = H5T.copy('H5T_NATIVE_FLOAT');
                case 'int64'
                    AttTypeID = H5T.copy('H5T_NATIVE_LLONG');
                case 'uint64'
                    AttTypeID = H5T.copy('H5T_NATIVE_ULLONG');
                case 'int32'
                    AttTypeID = H5T.copy('H5T_NATIVE_INT');
                case 'uint32'
                    AttTypeID = H5T.copy('H5T_NATIVE_UINT');
                case 'int16'
                    AttTypeID = H5T.copy('H5T_NATIVE_SHORT');
                case 'uint16'
                    AttTypeID = H5T.copy('H5T_NATIVE_USHORT');
                case 'int8'
                    AttTypeID = H5T.copy('H5T_NATIVE_SCHAR');
                case 'uint8'
                    AttTypeID = H5T.copy('H5T_NATIVE_UCHAR');
                case 'char'
                    AttTypeID = H5T.copy('H5T_C_S1');
                    if ~isempty(Value)
                        H5T.set_size(AttTypeID,numel(Value));
                    end
                    H5T.set_strpad(AttTypeID,'H5T_STR_NULLTERM');
                otherwise
                    error('Unknown H5T type');
            end
        end
         
        
    end % end methods
    
    % static methods for writing/reading HDF5 containing catalog data
    methods (Static)
        
        function save_cat(FileName,VarName,Data,SortCol,StepRows)
            % save catalog data in HDF5 file
            % Package: @HDF5
            % Description: save catalog data in HDF5 file
            %              Given a matrix containing a catalog, save the
            %              data in an HDF5 file. THe data will be saved
            %              under two variable names in the HDF5 file.
            %              /<base>_Cat will contain the catalog, while
            %              /<base>_Ind will contain an index data.
            %              The index data contains two columns [Ind Val],
            %              where Val is the values of the sorted column
            %              at the line index specified by Ind. Ind are in
            %              steps given by the StepRows parameter.
            % Input  : - File name
            %          - Base variable name.
            %            The actual name will be <base>_Cat and <base>_Ind.
            %          - Matrix containing the data to save
            %          - Column index by which to sort the catalog.
            %          - Number of rows step by for which to save the index
            %            data. Default is 30.
            % Outout : null
            % Example: HDF5.save_cat('try_cat.hdf5','V',Data,2,1000);
            % Reliable: 2
            
            if (nargin<5)
                StepRows = 30;
            end
            
            % prep Data
            Data  = sortrows(Data,SortCol);
            Nrows = size(Data,1);
            VecInd       = [(1:StepRows:Nrows), Nrows]';
            VecSortedCol = Data(VecInd,SortCol);
            IndexData    = [VecInd, VecSortedCol];
            
            % save index data
            VarNameInd = sprintf('/%s_Ind',VarName);
            HDF5.save(IndexData,FileName,VarNameInd);
            
            % save catalog
            SizeData = size(Data);
            Attrib = {'NCOL',SizeData(2); 'NROW',SizeData(1)};
            
            HDF5.save(Data,FileName,sprintf('/%s',VarName),Attrib);
            
        end
        
        function Cat=load_cat(FileName,VarName,SearchParValue,Ncol)
            % Load catalog stored in an HDF5 file
            % Package: @HDF5
            % Description: Load catalog stored in an HDF5 file. Given a
            %              and catalog in HDF5 file created by
            %              HDF5.save_cat, load the catalog. The catalog is
            %              sorted by one of the columns and it is possible
            %              to retrieve only line in some range. The search
            %              is done using the index data.
            % Input  : - HDF5 File name.
            %          - Variable name from which to load the catalog.
            %          - A two element vector of lower and upper value.
            %            Only lines in which the sorted parameter is
            %            between the low and high value will be retrieved.
            %          - Number of columns in the catalog.
            %            Default is empty (will attempt to find it).
            % Output : - A matrix containing the catalog
            % Example: Cat=HDF5.load_cat('try_cat6.hdf5','V',[0.1 0.101],20);
            % Reliable: 2
            
            Def.SearchParValue = [];
            Def.Ncol           = [];
            if (nargin<3)
                SearchParValue = Def.SearchParValue;
                Ncol           = Def.Ncol;
            elseif (nargin<4)
                Ncol           = Def.Ncol;
            else
                % do nothing
            end
            
            VarNameStr = sprintf('/%s',VarName);
            if (isempty(SearchParValue))
                % read entire catalog
                Cat = HDF5.load(FileName,VarNameStr);
            else
                % read index data first
                VarIndStr = sprintf('/%s_Ind',VarName);
                DataInd   = HDF5.load(FileName,VarIndStr);
                
                Ndi = size(DataInd,1);
                
                % search the index
                I1 = Util.find.bin_sear(DataInd(:,2),SearchParValue(1));
                I2 = Util.find.bin_sear(DataInd(:,2),SearchParValue(2));
                
                if (isempty(Ncol))
                    % get number of columns from HDF5 file attributes
                    error('Get Ncol from attributes not implemented yet');
                end
                
                % read data
                Offset = [DataInd(I1,1), 1];
                if (I1==I2)
                    I2 = I2 + 1;
                end
                I2 = min(I2,Ndi);
                
                Block  = [DataInd(I2,1)-DataInd(I1,1), Ncol];
                Cat = HDF5.load(FileName,VarNameStr,Offset,Block);
                
            end
            
        end

        function Nsrc=get_nsrc(CatName)
            % Count number of sources over all HTM in HDF5 files
            % Package: @HDF5
            % Input  : - Catalog name (e.g., 'APASSS')
            % Output : - A matrix of [HTM_index, Nsrc]
            % Example: Nsrc=HDF5.get_nsrc(CatName);
            % Reliable: 2
            
            Dir = dir(sprintf('%s_htm_*.hdf5',CatName));
            Ndir = numel(Dir);
            Nsrc = zeros(100.*Ndir,2);
            K = 0;
            for Idir=1:1:Ndir
                
                Info = h5info(Dir(Idir).name);
                IndH = find(cellfun(@numel,strfind({Info.Datasets.Name},'_'))==1);
                Nih  = numel(IndH);
                for Iih=1:1:Nih
                    K = K + 1;
                    IndHTM = str2double(Info.Datasets(IndH(Iih)).Name(5:end));
                    Nsrc(K,:) = [IndHTM size(h5read(Dir(Idir).name,sprintf('/%s',Info.Datasets(IndH(Iih)).Name)),1)];
                end
            end
            Nsrc = Nsrc(1:K,:);
            
        end
        
        function save_htm_ind(HTM,FileName,VarName,Attrib,Nsrc)
            % Save HTM indinces of the celestial sphere in an HDF5 file
            % Package: @HDF5
            % Description: Generate HDF5 file with HTM indices.
            %              The HTM indices contains the HTM tree and the 3
            %              poles of the 3 great circles that defines each
            %              HTM.
            % Input  : - Either a structure of HTM to save (created by
            %            celestial.htm.htm_build), or the HTM level.
            %          - HDF5 File name in which to store the HTM indices.
            %          - Variable name in which to store the data.
            %            Default is '<CatName>_HTM'
            %          - Cell array of attribute {Key,Val} to save
            %            in a 'ColCell' variable name.
            %          - Nsrc matrix [IndHTM Nsrc]
            % Output : null
            % Example: HDF5.save_htm_ind(7,'try_htm.hdf5',[],{},Nsrc)
            % Reliable: 2
            
            Tmp = regexp(FileName,'_','split');
            Def.HTM = sprintf('%s_HTM',Tmp{1});
            
            if (nargin<3)
                VarName = Def.HTM;
                Attrib  = {};
                Nsrc    = [];
            end
            
            if (isempty(VarName))
                VarName = Def.HTM;
            end
            
            if (isnumeric(HTM))
                % generate HTM index
                [HTM]=celestial.htm.htm_build(HTM);
            end
            
            Nhtm = numel(HTM);
            
            Data = nan(Nhtm,13);
            for Ihtm=1:1:Nhtm
                Nlev = numel(HTM(Ihtm).id);
                ID   = sum(logspace(1,Nlev,Nlev).*HTM(Ihtm).id)./10;
                % Level, Father, Son1, Son2, Son3, Son4, Poles 1 long, poles 1 lat, ...
                if isempty(HTM(Ihtm).father)
                    Father = NaN;
                else
                    Father = HTM(Ihtm).father;
                end
                if (isempty(HTM(Ihtm).son))
                    Son = [NaN NaN NaN NaN];
                else
                    Son = HTM(Ihtm).son;
                end
                
                if (isempty(Nsrc))
                    Ns = NaN;
                else
                    Ns = Nsrc(Nsrc(:,1)==Ihtm,2);
                    if (isempty(Ns))
                        Ns = NaN;
                    end
                end
                Data(Ihtm,:) = [HTM(Ihtm).level, Father, Son, HTM(Ihtm).PolesCoo(1,:), HTM(Ihtm).PolesCoo(2,:), HTM(Ihtm).PolesCoo(3,:), Ns];
            end
            
            % save HTM
            AttribHTM = {'Table.Col.1','Level'; ...
                      'Table.Col.2','Father'; ...
                      'Table.Col.3','Son1'; ...
                      'Table.Col.4','Son2'; ...
                      'Table.Col.5','Son3'; ...
                      'Table.Col.6','Son4'; ...
                      'Table.Col.7', 'Poles1Lon';...
                      'Table.Col.8', 'Poles1Lat';...
                      'Table.Col.9', 'Poles2Lon';...
                      'Table.Col.10','Poles2Lat';...
                      'Table.Col.11','Poles3Lon';...
                      'Table.Col.12','Poles3Lat';...
                      'Table.Col.13','Nsrc'};
            HDF5.save(single(Data),FileName,VarName,AttribHTM);
            % save column names
            HDF5.save([],FileName,'ColNames',Attrib);
                
        end
        
        function save_cat_colcell(CatName,ColCell,ColUnits)
            % Save ColCell cell array of an HTM catalog into a MAT file
            % Package: @HDF5
            % Input  : - Catalog name (e.g., 'APASS')
            %          - ColCell cell array
            %          - ColUnits cell array
            % Reliable : 2
            
            if (nargin<3)
                ColUnits = {};
            end
            FileName = sprintf('%s_htmColCell.mat',CatName);
            save(FileName,'ColCell','ColUnits')
            
            
        end
        
        function [ColCell,Col]=read_colnames(FileName,VarName)
            % read HDF5 catalog column names from index file
            % Package: @HDF5
            % Input  : - HDF5 file name.
            %          - Variable name. Default is '/ColNames'.
            % Output : - Cell array of column names.
            %          - Structure array of column indices.
            % Example: [ColCell,Col]=HDF5.read_colnames('GAIADR1_htm.hdf5');
            
            if (nargin<2)
                VarName = '/ColNames';
            end
            
            Ncol = h5readatt('GAIADR1_htm.hdf5','/ColNames','Table.Ncol');
            ColCell = cell(1,Ncol);
            for Icol=1:1:Ncol
                ColCell{Icol} = h5readatt('GAIADR1_htm.hdf5','/ColNames',sprintf('Table.Col.%d',Icol));
                Col.(ColCell{Icol}) = Icol;
            end
        end
       
        function HTM=load_htm_ind(FileName,VarName)
            % load HTM data into structure from an HDF5 file
            % Package: @HDF5
            % Description: load HTM data into structure from an HDF5 file
            % Input  : - HDF5 file name containing the HTM data.
            %          - Variable name. Default is '<CatName>_HTM'.
            % Output : - A structure array containing the HTM structure.
            % Example: HTM=HDF5.load_htm_ind('try_htm.hdf5','HTM');
            % Reliable :2 
            
            if (nargin<2)
                Tmp = regexp(FileName,'_','split');
                VarName = sprintf('%s_HTM',Tmp{1});
            
            end
            
            % read data from HDF5 file
            Data = HDF5.load(FileName,VarName);
            
            % load into HTM structure
            Nhtm = size(Data,1);
            HTM  = Util.struct.struct_def({'level','father','son','PolesCoo'},1,Nhtm);
            for Ihtm=1:1:Nhtm
                HTM(Ihtm).level = Data(Ihtm,1);
                %HTM(Ihtm).id    = [];
                %HTM(Ihtm).coo   = [];
                %HTM(Ihtm).cosd  = [];
                %HTM(Ihtm).center_coo = [];
                %HTM(Ihtm).center_cosd = [];
                if (isnan(Data(Ihtm,2)))
                    HTM(Ihtm).father  = [];
                else
                    HTM(Ihtm).father  = Data(Ihtm,2);
                end
                if (isnan(Data(Ihtm,3)))
                    HTM(Ihtm).son  = [];
                else
                    HTM(Ihtm).son  = Data(Ihtm,3:6);
                end
                HTM(Ihtm).PolesCoo = [Data(Ihtm,7:8); Data(Ihtm,9:10); Data(Ihtm,11:12)];
                
            end
            
        end
        
        function ID=search_htm_ind(FileName,VarName,Long,Lat,Radius)
            % A coordinate cone search in an HTM stored in HDF5 file.
            % Package: @HDF5
            % Description: A coordinate cone search in an HTM stored in HDF5 file.
            %              See also: celestial.htm.htm_search_cone
            % Input  : - An HDF5 file name or an open HDF5 object, in which
            %            the HTM indices are stored.
            %          - Variable name. If empty, default is <CatName>_HTM.
            %          - Search longitude [radians].
            %          - Search latitude [radians].
            %          - Search radius [radians].
            % Example: ID=HDF5.search_htm_ind('UCAC4_htm.hdf5',[],1,1,0.001)
            % Reliable: 2
            
            
            Check = true;
            if (isempty(VarName))
                Tmp = regexp(FileName,'_','split');
                VarName = sprintf('%s_HTM',Tmp{1});
            end
                     
            if (Check)
                DataHTM = HDF5.load_check(FileName,VarName);
            else
                DataHTM = HDF5.load(FileName,VarName);
            end
            
            ID=HDF5.htm_search_cone(DataHTM,Long,Lat,Radius);
            
            % check that HTM contains sources
            ID = ID(DataHTM(ID,13)>0);

        end
       
        
        function ID=htm_search_cone(DataHTM,Long,Lat,Radius,Ind)
            % Search for all HTM leafs interscting a small circle (cone search)
            % Package: celestial.htm
            % Description: Search for all HTM leafs interscting a small circle
            %              (i.e., cone search).
            % Input  : - Either a table of HTM data or an open HDF5 object
            %            in which the HTM data is stored.
            %          - Longitude [radians] to search.
            %          - Latitude [radians] to search.
            %          - Search radius [radians].
            % Example:  [HTM,LevList]=celestial.htm.htm_build(4);
            %           ID=HDF5.htm_search_cone(HTM,1,1,0.0001)
            % Reliable : 2
            
            Col.Father = 2;
            Col.Son    = [3 4 5 6];
            Col.PolesLong  = [7 9  11];
            Col.PolesLat   = [8 10 12];

            if (nargin<5)
                Ind = [];
            end


            if isempty(Ind)
                % first iteration
                Sons  = (1:1:8);
                %Nsons = 8;
            else
                Sons  = Ind;
                %Nsons = 4;
            end

            ID = [];
            Nsons = numel(Sons);
            PolesLong = zeros(3,Nsons);
            PolesLat  = zeros(3,Nsons);

            % DataHTM is the full HTM table
            for Isons=1:1:Nsons
                %CSon = Sons(Isons);
                PolesLong(:,Isons) = DataHTM(Sons(Isons),Col.PolesLong); %   HTM(Sons(Isons)).PolesCoo(:,1);
                PolesLat(:,Isons)  = DataHTM(Sons(Isons),Col.PolesLat);  % HTM(Sons(Isons)).PolesCoo(:,2);
            end
            Flag = celestial.htm.cone_in_polysphere(PolesLong,PolesLat,Long,Lat,Radius);


            for Isons=1:1:Nsons
                if (Flag(Isons))
                    % cone overlap HTM
                    CSon = Sons(Isons);
                    if isnan(DataHTM(CSon,Col.Son))
                        % arrived at last leaf
                        % return ID
                        ID = [ID, CSon];
                    else
                        Ind = DataHTM(CSon,Col.Son); % HTM(CSon).son;
                        ID  = [ID, HDF5.htm_search_cone(DataHTM,Long,Lat,Radius,Ind)];
                        %ID = cat(2,ID,celestial.htm.htm_search_cone(HTM,Long,Lat,Radius,Ind));
                    end
                end
            end  
            
        end
        
        function [FileName,DataName]=get_file_var_from_htmid(FileBase,ID,NfilesInHDF)
            % Construct file and var name for HTM file stored in HDF5
            % Package: @HDF5
            % Description: Given a file base (e.g., 'UCAC4') and HTM ID
            %              and number of files in HDF5 file, construct the
            %              HDF5 file name (e.g., UCAC4_htm_032400.hdf5),
            %              and the data variable name (e.g., htm_032412).
            % Input  : - Catalog base name (e.g., 'UCAC4').
            %          - HTM index.
            %          - Number of variables in file (default is 100).
            % Output : - File name.
            %          - Variable name.
            % Example: [FileName,DataName]=HDF5.get_file_var_from_htmid('UCAC4',45661,100)
            % Reliable: 2
           
            
            if (nargin<3)
                NfilesInHDF = 100;
            end
            
            FileID    = floor(ID./NfilesInHDF).*NfilesInHDF;
            FileName  = sprintf('%s_htm_%06d.hdf5',FileBase,FileID);
            DataName  = sprintf('htm_%06d',ID);
        end
        
    end % end methods
    
    % online variable
    methods (Static)
        function H=access_open(FileName,VarName)
            % Open an HDF5 file/dataset for reading
            % Input  : - HDF5 file name.
            %          - HDF5 dataset name.
            % Output : - An HDF5 object containing the HDF5 and dataset ID.
            
            
            H=HDF5;
            Info = h5info(FileName,VarName);
            H.Size   = Info.Dataspace.Size;
            H.File   = FileName;
            H.Var    = VarName;
            H.FID    = H5F.open(FileName); 
            H.DSetID = H5D.open(H.FID,VarName);
            
        end
    end
        
    % non static
    methods
            
        function access_close(H)
            % close and HDF5/dataset 
            % Input  : - An HDF5 object
            % Output : null
            
            H5D.close(H.DSetID);
            H5F.close(H.FID);
        end
            
        function Data=get_data(H,Offset,Block)
            % Get data from an HDF5 object
            % Description: After creation of an HDF5 object, you can access
            %              the data in the corresponding HDF5/dataset using
            %              this function.
            % Input  : - An HDF5 object.
            %          - Data offset.
            %          - data block.
            % Output : - Data
            % Example: get_data(H,[1000 1],[10 20])
            
            
            % User requested for specific location in the array
            Plist = 'H5P_DEFAULT';
            %Space = H5D.get_space(DSetID);
            %[~,Dims] = H5S.get_simple_extent_dims(Space);
                
            Dims = fliplr(Block);
            MemSpaceID = H5S.create_simple(numel(Dims),Dims,[]);
            FileSpaceID = H5D.get_space(H.DSetID);
            Offset = fliplr(Offset-1);  % index start at 0
            Block = fliplr(Block);
            H5S.select_hyperslab(FileSpaceID,'H5S_SELECT_SET',Offset,[],[],Block);
            Data = H5D.read(H.DSetID,'H5ML_DEFAULT',MemSpaceID,FileSpaceID,Plist);
        end
        
    end % end methods
    
    methods
        
    end % end methods
end

            
