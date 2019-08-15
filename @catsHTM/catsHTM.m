% catsHTM static class. Read / write HDF5/HTM catalogs
% Package: @catsHTM
% Description: A static class for catsHTM related functions.
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef catsHTM 
             
    
    % file and variable names
    methods (Static)
        function CatName=filename2base(FileName)
            % Convert HDF5/HTM file name to catalog name (file base name)
            % Package: @catsHTM
            % Description: Convert HDF5/HTM file name (e.g., 'APASS_htm.hdf5')
            %              to catalog name (file base name; e.g., 'APASS').
            % Input  : - HDF5 file name that contains the catalog name.
            %            The file name is composed of strings seperated by
            %            "_", where the first string is the catalog name.
            % Output : - Catalog name.
            % Example: CatName=catsHTM.filename2base('SDSSDR10_htm.hdf5')
            % Reliable: 
            
             Tmp = regexp(FileName,'_','split');
             CatName = Tmp{1};
             
        end
        
        function [FileName,DataName]=get_file_var_from_htmid(FileBase,ID,NfilesInHDF)
            % Construct file and var name for HTM file stored in HDF5
            % Package: @catsHTM
            % Description: Given a file base (e.g., 'UCAC4') and HTM ID
            %              and number of files in HDF5 file, construct the
            %              HDF5 file name (e.g., UCAC4_htm_032400.hdf5),
            %              and the data variable name (e.g., htm_032412).
            % Input  : - Catalog base name (e.g., 'UCAC4').
            %          - HTM index.
            %          - Number of variables in file (default is 100).
            % Output : - File name.
            %            If ID is a vector then this is a cell array.
            %          - Variable name.
            %            If ID is a vector then this is a cell array.
            % Example: [FileName,DataName]=catsHTM.get_file_var_from_htmid('UCAC4',45661,100)
            % Reliable: 2


            if (nargin<3)
                NfilesInHDF = 100;
            end

            FileID    = floor(ID./NfilesInHDF).*NfilesInHDF;
            Nid       = numel(FileID);
            if (Nid==1)
                FileName  = sprintf('%s_htm_%06d.hdf5',FileBase,FileID);
                DataName  = sprintf('htm_%06d',ID);
            else
                FileName = cell(Nid,1);
                DataName = cell(Nid,1);
                for Iid=1:1:Nid
                    FileName{Iid}  = sprintf('%s_htm_%06d.hdf5',FileBase,FileID(Iid));
                    DataName{Iid}  = sprintf('htm_%06d',ID(Iid));
                end
            end
        end
       
        function [IndexFileName,IndexVarName]=get_index_filename(CatName)
            % Get HDF5/HTM index file name and variable name from CatName
            % Package: @catsHTM
            % Description: Get HDF5/HTM index file name and variable name
            %              from CatName.
            % Input  : - Catalog name (e.g., 'APASS').
            % Output : - File name (e.g., 'APASS_htm.hdf5')
            %          - Variable name (e.g., 'APASS_HTM')
            % Example: [IndexFileName,IndexVarName]=catsHTM.get_index_filename('PS1')
            % Reliable: 2
            
            IndexFileName = sprintf('%s_htm.hdf5',CatName);
            IndexVarName = sprintf('%s_HTM',CatName);
            
        end
        
        
        
    end
    
    % Create/save data into HDF5/HTM related files
    methods (Static)
           
        function save_cat(FileName,VarName,Data,SortCol,StepRows)
            % save catalog data in HDF5 file
            % Package: @catsHTM
            % Description: save catalog data in HDF5 file
            %              Given a matrix containing a catalog, save the
            %              data in an HDF5 file. The data will be saved
            %              under two variable names in the HDF5 file.
            %              /<base>_Cat will contain the catalog, while
            %              /<base>_Ind will contain an index data.
            %              The index data contains two columns [Ind Val],
            %              where Val is the values of the sorted column
            %              at the line index specified by Ind. Ind are in
            %              steps given by the StepRows parameter.
            % Input  : - File name
            %          - Base variable name.
            %            The actual name will be <base> and <base>_Ind.
            %          - Matrix containing the data to save
            %          - Column index by which to sort the catalog.
            %          - Number of rows step by for which to save the index
            %            data. Default is 30.
            % Outout : null
            % Example: catsHTM.save_cat('try_cat.hdf5','V',Data,2,1000);
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
        
        function save_htm_ind(HTM,FileName,VarName,Attrib,Nsrc)
            % Save HTM indinces of the celestial sphere in an HDF5 file
            % Package: @catsHTM
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
            % Example: catsHTM.save_htm_ind(7,'try_htm.hdf5',[],{},Nsrc)
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
                % Level, Father, Son1, Son2, Son3, Son4, Poles 1 long,
                % poles 1 lat, ..., Nsrc
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
            % Save ColCell cell array of an HTM catalog
            % Package: @catsHTM
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
        
        function count_edge_in_cat(CatName,SearchRadius,NfilesInHDF)
            %
            % Example: catsHTM.count_edge_in_cat('APASS');
            RAD = 180./pi;
            
            if (nargin<3)
                NfilesInHDF = 100;
                if (nargin<2)
                    SearchRadius = 3./3600./RAD;
                end
            end
            
            % load HTM index
            [IndexFileName,IndexVarName]=catsHTM.get_index_filename(CatName);
            [~,DataHTM] = catsHTM.load_htm_ind(IndexFileName,IndexVarName);
            Level=celestial.htm.nhtm2level(size(DataHTM,1));
            HTM = celestial.htm.htm_build(Level);
            
            % for each HTM that contain sources
            Ihtm = find(DataHTM(:,13)>0);
            Nh = numel(Ihtm);
            for Ih=1:1:Nh
                tic;
                I = Ihtm(Ih);
                
                % load catalog of HTM
                Cat = catsHTM.load_cat(CatName,I);
                Nsrc = size(Cat,1);
                % Search for all sources in HTM tile that are near the
                % edges.
                FlagEdge = false(Nsrc,1);
                for Isrc=1:1:Nsrc
                    FlagEdge(Isrc) = numel(celestial.htm.htm_search_cone(HTM,Cat(Isrc,1),Cat(Isrc,2),SearchRadius))>1;
                end
                toc
                sum(FlagEdge)
                
            end
            
        end
        
        function generate_edge_cat(CatName,SearchRadius,NfilesInHDF)
            % OBSOLOTE
            
            RAD = 180./pi;
            if (nargin<3)
                NfilesInHDF = 100;
                if (nargin<2)
                    SearchRadius = 5./3600./RAD;
                end
            end
            
            % load HTM index
            [IndexFileName,IndexVarName]=catsHTM.get_index_filename(CatName);
            [~,DataHTM] = catsHTM.load_htm_ind(IndexFileName,IndexVarName);
            Level=celestial.htm.nhtm2level(size(DataHTM,1));
            HTM = celestial.htm.htm_build(Level);
            
            % for each HTM that contain sources
            Ihtm = find(DataHTM(:,13)>0);
            Nh = numel(Ihtm);
            for Ih=1:1:Nh
                I = Ihtm(Ih);
                
                % search for all HTMs in Cat2 that may opverlap with
                % Cat1 current HTM
                MeanRA  = mean(HTM(I).coo(:,1));
                MeanDec = mean(HTM(I).coo(:,2));

                D = celestial.coo.sphere_dist_fast(MeanRA,MeanDec,HTM(I).coo(:,1),HTM(I).coo(:,2));
                CircRadius = max(D) + SearchRadius; % [rad]

                ID2 = celestial.htm.htm_search_cone(HTM,MeanRA,MeanDec,CircRadius);

                % load all ID2 from HTM2
                Nid2 = numel(ID2);
                for Iid2=1:1:Nid2
                    if (Iid2==1)
                        Cat   = catsHTM.load_cat(CatName,ID2(Iid2));
                        N2     = size(Cat,1);
                    else
                        Cat   = [Cat; catsHTM.load_cat(CatName,ID2(Iid2))];
                        N2     = size(Cat,1);
                    end
                end

                % search for sources in edge of HTM
                FlagInHTM   = celestial.htm.in_polysphere(Cat(:,1:2),HTM(I).coo);
                FlagNearHTM = celestial.htm.cone_in_polysphere(HTM(I).PolesCoo(:,1),HTM(I).PolesCoo(:,2),Cat(:,1),Cat(:,2),SearchRadius);
                FlagEdge    = ~FlagInHTM(:) & FlagNearHTM(:);
                %sum(FlagEdge)
                % store the sources
                [FileName,DataName]=catsHTM.get_file_var_from_htmid(CatName,I,NfilesInHDF);
                HDF5.save(Cat(FlagEdge,:),FileName,sprintf('/htm_%06d_Edge',I));
                
            end
        end
        
        function create_catalog_lists4wget(Dir,WriteDir)
            % Create list of catalogs foe wget including checsums
            % Input  : - Directory in which the catsHTM catalog resides
            %            (e.g., '/raid/eran/catsHTM').
            %          - Directory in which to write wget lists.
            %            Default is '' - i.e., current dir.
            % Example:
            % catsHTM.create_catalog_lists4wget('/data/euler/catsHTM','/home/eran/');
           
            if (nargin<2)
                WriteDir = '';
            end
            
            URL  = 'https://astro.weizmann.ac.il/catsHTM/';
            Pars = '-U Mozilla/5.0 --no-check-certificate';
            Nc   = numel(Dir);
            
            PWD = pwd;
            cd(Dir);
            
            F1 = Util.IO.rdir('*.hdf5');
            F2 = Util.IO.rdir('*.mat');
            F  = [F1;F2];
            
            Nf = numel(F);
            FIDw = fopen(sprintf('%s%s%s',WriteDir,filesep,'list.euler.wget'),'w');
            FIDc = fopen(sprintf('%s%s%s',WriteDir,filesep,'list.euler.checksum'),'w');
            tic;
            for If=1:1:Nf
                Pars1 = sprintf('%s -P .%s',Pars,F(If).folder(Nc+1:end));
                fprintf(FIDw,'wget %s %s%s/%s\n',Pars1,URL,F(If).folder(Nc+1:end),F(If).name);
                [~,Str] = system(sprintf('md5sum %s%s%s',F(If).folder,filesep,F(If).name));
                fprintf(FIDc,'%s',Str);
            end
            fclose(FIDw);
            fclose(FIDc);
            toc
            
            cd(PWD);
            
        end
        
    end
    
    % Load and search HDF5/HTM files
    methods (Static)
        
        function [Cat,Ind]=load_cat(FileName,VarName,SearchParValue,Ncol,NfilesInHDF)
            % Load catalog stored in an HDF5 file
            % Package: @catsHTM
            % Description: Load catalog stored in an HDF5 file. Given a
            %              and catalog in HDF5 file created by
            %              HDF5.save_cat, load the catalog. The catalog is
            %              sorted by one of the columns and it is possible
            %              to retrieve only line in some range. The search
            %              is done using the index data.
            % Input  : - HDF5 File name, or catalog name.
            %            If catalog name, then second argument must be
            %            numeric index.
            %          - Variable name from which to load the catalog,
            %            or a numeric HTM index.
            %          - A two element vector of lower and upper value.
            %            Only lines in which the sorted parameter is
            %            between the low and high value will be retrieved.
            %            If empty, retrieve all lines. Default is empty.
            %          - Number of columns in the catalog.
            %            Default is empty (will attempt to find it).
            %          - Number of HTM data matrix in each hdf5 file.
            %            Default is 100.
            % Output : - A matrix containing the catalog
            % Example: Cat=catsHTM.load_cat('APASS_htm_010000.hdf5','htm_010001',[0.1 0.101],20);
            %          Cat=catsHTM.load_cat('APASS_htm_043600.hdf5','htm_043601');
            %          Cat=catsHTM.load_cat('APASS',43601);
            % Reliable: 2
            
            Def.SearchParValue = [];
            Def.Ncol           = [];
            Def.NfilesInHDF    = 100;
            if (nargin<3)
                SearchParValue = Def.SearchParValue;
                Ncol           = Def.Ncol;
                NfilesInHDF    = Def.NfilesInHDF;
            elseif (nargin<4)
                Ncol           = Def.Ncol;
                NfilesInHDF    = Def.NfilesInHDF;
            elseif (nargin<5)
                NfilesInHDF    = Def.NfilesInHDF;
            else
                % do nothing
            end
            
            if (isnumeric(VarName))
                % assume FileName is CatName and VarName in HTM index
                [FileName,VarName] = catsHTM.get_file_var_from_htmid(FileName,VarName,NfilesInHDF);
            end
            
            VarNameStr = sprintf('/%s',VarName);
            if (isempty(SearchParValue))
                % read entire catalog [only if exist]
                Ind = 1;
                try
                    Cat = HDF5.load(FileName,VarNameStr);
                catch
                    Cat = [];
                end
            else
                % read index data first
                try
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
                    Ind    = DataInd(I1,1);
                    Offset = [DataInd(I1,1), 1];
                    if (I1==I2)
                        I2 = I2 + 1;
                    end
                    I2 = min(I2,Ndi);

                    Block  = [1+DataInd(I2,1)-DataInd(I1,1), Ncol];
                    Cat = HDF5.load(FileName,VarNameStr,Offset,Block);
                catch
                    Ind = 1;
                    Cat = [];
                end
            end
            
        end

        function Cat=load_multiple_cats(CatName,ID,NfilesInHDF)
            % Load HDF5/HTM catalog from multiple files/datasets
            % Package: @catsHTM
            % Description: Load HDF5/HTM catalog from multiple files/datasets
            %              Not as fast as expected.
            % Input  : - CatName
            %          - Vector of htm indices.
            %          - Number of datasets in HDF5 file. Default is 100.
            % Output : - Joint catalog.
            % Example: Data=catsHTM.load_multiple_cats('UCAC4',[19100:1:191002]')
            % Reliable: 2
            
            if (nargin<3)
                NfilesInHDF = 100;
            end
            
            %Nid = numel(ID);
            % get file/dataset name for all IDs
            [FileName,DataName] = catsHTM.get_file_var_from_htmid(CatName,ID,NfilesInHDF);
            
            % select unique files 
            FileID = floor(ID./NfilesInHDF).*NfilesInHDF;
            UniqueFID = unique(FileID);
            Nufid = numel(UniqueFID);
            for Iufid=1:1:Nufid
                % select all IDs in file
                Ifile = find(UniqueFID(Iufid)==FileID);
                if (Iufid==1)
                    Cat = HDF5.load_muti_datasets(FileName{Ifile(1)}, DataName(Ifile));
                else
                    Cat = [Cat; HDF5.load_muti_datasets(FileName{Ifile(1)}, DataName(Ifile))];
                end
            end
            
        end
        
        function [Cat,EdgeOk]=load_cat_edge(CatName,IndHTM,NfilesInHDF)
            % Load and concat HDF5/HTM catalog and its edge catalog 
            % Package: @catsHTM
            % Description: Load and concat HDF5/HTM catalog and its edge catalog 
            % Input  : - Catalog base name.
            %          - HTM index.
            %          - Number of HTM datasets in file. Default is 100.
            % Output : - THe catalog.
            %          - A logical flag indicating if the Edge catalog was
            %            sucessfully uploaded.
            % Example: Cat=catsHTM.load_cat_edge('APASS',19000);
            % Reliable: 2
            
            ColDec = 2;
            if (nargin<3)
                NfilesInHDF = 100;
            end
            
            [FileName,DataName]=catsHTM.get_file_var_from_htmid(CatName,IndHTM,NfilesInHDF);
            EdgeOK = true;
            try
                Cat  = HDF5.load(FileName,sprintf('/%s',DataName));
            catch
                Cat = [];
            end
            if (~isempty(Cat))
                try
                    CatE = HDF5.load(FileName,sprintf('/%s_Edge',DataName));
                catch
                    CatE = [];
                    EdgeOK = false;
                end
            else
                CatE = [];
            end
            Cat = [Cat; CatE];
            
            if (~isempty(Cat))
                Cat = sortrows(Cat,ColDec);
            end
            
        end
        
        function [Cat,ColCell]=load_1htm(CatName,IndexHTM,NfilesInHDF)
            % Load a single tile of HDF5/HTM catalog
            % Package: @catsHTM
            % Description: Load a single HTM tile of HDF5/HTM catalog based
            %              on its HTM index.
            %              This is slower relative to catsHTM.load_cat,
            %              since it also loads the index file.
            % Input  : - Catalog name (e.g., 'APASS').
            %          - HTM index.
            %          - Number of data varaible in HDF5 file.
            %            Default is 100.
            % Output : - Catalog matrix.
            %          - Cell array of column names.
            % Example: [Cat,ColCell]=catsHTM.load_1htm('APASS',25000)
            % Reliable: 2
            
            if (nargin<3)
                NfilesInHDF = 100;
            end
            
            FileName = sprintf('%s_htm.hdf5',CatName);
            DataName = sprintf('%s_HTM',CatName);
            Data     = HDF5.load(FileName,DataName);
            
            if ~(IndexHTM>0 && IndexHTM<=size(Data,1))
                error('IndexHTM was not found in index file');
            end
            if (Data(IndexHTM,13)>0)
                [FileName,DataName]=catsHTM.get_file_var_from_htmid(CatName,IndexHTM,NfilesInHDF);
                Cat = catsHTM.load_cat(FileName,DataName);
            else
                Cat = [];
            end
            
            if (nargout>1)
                File = sprintf('%s_htmColCell.mat',CatName);
                load(File);
            end
            
        end
        
        function [ColCell,ColUnits,Col] = load_colcell(CatName)
            % Load ColCell and ColUnits for an HDF5/HTM catalog
            % Package: @catsHTM
            % Input  : - Catalog base name (e.g., 'DECaLS').
            % Output : - Cell array of column names.
            %          - Cell array of column units
            %          - Structure with column names and indices
            % Example: [ColCell,ColUnits]=catsHTM.load_colcell('APASS')
            % Reliable: 2
            
            File = sprintf('%s_htmColCell.mat',CatName);
            load(File);
            
            if (nargout>2)
                Col = cell2struct(num2cell(1:1:numel(ColCell)),ColCell,2)
            end
        end
        
        function [ColCell,Col]=read_colnames(FileName,VarName)
            % read HDF5 catalog column names from index file
            % Package: @catsHTM
            % Input  : - HDF5 file name.
            %          - Variable name. Default is '/ColNames'.
            % Output : - Cell array of column names.
            %          - Structure array of column indices.
            % Example: [ColCell,Col]=catsHTM.read_colnames('GAIADR1_htm.hdf5');
            
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
       
        function [HTM,Data]=load_htm_ind(FileName,VarName)
            % load HTM data into structure from an HDF5 file
            % Package: @catsHTM
            % Description: load HTM data into structure from an HDF5 file
            % Input  : - HDF5 file name containing the HTM data.
            %          - Variable name. Default is '<CatName>_HTM'.
            % Output : - A structure array containing the HTM structure.
            %          - Thr matrix containing the HTM data.
            % Example: HTM=catsHTM.load_htm_ind('try_htm.hdf5','HTM');
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
            % Package: @catsHTM
            % Description: A coordinate cone search in an HTM stored in HDF5 file.
            %              See also: celestial.htm.htm_search_cone
            % Input  : - An HDF5 file name or an open HDF5 object, in which
            %            the HTM indices are stored.
            %          - Variable name. If empty, default is <CatName>_HTM.
            %          - Search longitude [radians].
            %          - Search latitude [radians].
            %          - Search radius [radians].
            % Example: ID=catsHTM.search_htm_ind('UCAC4_htm.hdf5',[],1,1,0.001)
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
            
            ID=catsHTM.htm_search_cone(DataHTM,Long,Lat,Radius);
            
            % check that HTM contains sources
            ID = ID(DataHTM(ID,13)>0);

        end
       
        function ID=htm_search_cone(DataHTM,Long,Lat,Radius,Ind)
            % Search for all HTM leafs interscting a small circle (cone search)
            % Package: @catsHTM
            % Description: Search for all HTM leafs interscting a small circle
            %              (i.e., cone search).
            % Input  : - Either a table of HTM data or an open HDF5 object
            %            in which the HTM data is stored.
            %          - Longitude [radians] to search.
            %          - Latitude [radians] to search.
            %          - Search radius [radians].
            % Example:  [HTM,LevList]=celestial.htm.htm_build(4);
            %           ID=catsHTM.htm_search_cone(HTM,1,1,0.0001)
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
                        ID  = [ID, catsHTM.htm_search_cone(DataHTM,Long,Lat,Radius,Ind)];
                        %ID = cat(2,ID,celestial.htm.htm_search_cone(HTM,Long,Lat,Radius,Ind));
                    end
                end
            end  
            
        end
        
        
      
    end % Static
    
    % utilities
    methods (Static)
        function Nsrc=get_nsrc(CatName)
            % Count number of sources over all HTM in HDF5 files
            % Package: @catsHTM
            % Input  : - Catalog name (e.g., 'APASSS')
            % Output : - A matrix of [HTM_index, Nsrc]
            % Example: Nsrc=catsHTM.get_nsrc(CatName);
            % Reliable: 2
            
            Dir = dir(sprintf('%s_htm_*.hdf5',CatName));
            Ndir = numel(Dir);
            Nsrc = zeros(100.*Ndir,3);
            K = 0;
            
            for Idir=1:1:Ndir
                Info = h5info(Dir(Idir).name);
                IndH = find(cellfun(@numel,strfind({Info.Datasets.Name},'_'))==1);
                Nih  = numel(IndH);
                for Iih=1:1:Nih
                    K = K + 1;
                    IndHTM = str2double(Info.Datasets(IndH(Iih)).Name(5:end));
                    Nsrc(K,:) = [IndHTM size(h5read(Dir(Idir).name,sprintf('/%s',Info.Datasets(IndH(Iih)).Name)),1),Idir];
                end
            end
            Nsrc = Nsrc(1:K,:);
            
        end
        
        function [Nsrc,SumN]=nsrc(CatName)
            % Count sources in the HDF5/HTM index file
            % Package: @catsHTM
            % Description: Count sources in the HDF5/HTM index file
            % Input  : - Catalog name (e.g., 'SDSSDR10').
            % Output : - Matrix of [HTMindex, Nsrc].
            %          - Total number of sources in catalog.
            % Example: [Nsrc,SumN]=catsHTM.nsrc('SDSSDR10');
            % Reliable: 2
            
            FileName = sprintf('%s_htm.hdf5',CatName);
            DataName = sprintf('%s_HTM',CatName);
            %HTM = catsHTM.load_htm_ind(FileName);
            Data = HDF5.load(FileName,DataName);
            
            Nsrc = Data(:,[2 13]);
            SumN = nansum(Nsrc(:,2));
            
        end
        
        function Ref=reference(CatName)
            % Get references for an HDF5/HTM catalog
            % Package: @catsHTM
            % Description: Get references for an HDF5/HTM catalog
            % Input  : - Catalog base name (e.g., 'GAIADR1').
            % Output : - Structure containing reference and acknowledgment
            %            for the catalog.
            % Example: catsHTM.reference('SDSSDR10')
            %
            
            switch lower(CatName)
                case 'SDSSDR10'
                    Ref.CatName = 'SDSSDR10';
                    Ref.Name    = 'SDSS-DR10 sources';
                    Ref.Ref{1}  = 'Ahn et al. (2014)';
                    Ref.Link{1} = 'http://adsabs.harvard.edu/abs/2014ApJS..211...17A';
                    Ref.Ack     = 'http://www.sdss3.org/collaboration/boiler-plate.php';
                case 'GAIADR1'
                case 'GALEX'
                case 'DECaLS'
                case 'TMASS'
                case 'WISE'
                case 'FIRST'
                case 'NVSS'
                case 'APASS'
                case 'UCAC4'
                case 'XMM'
                    
                otherwise
                    error('Unknown CatName option');
            end
            Ref.CatAck  = 'HDF5/HTM catalog from Ofek (2018)';
            
        end
        
    end
    
    % search
    methods (Static)
        function [Cat,ColCell,ColUnits]=cone_search(CatName,RA,Dec,Radius,varargin)
            % Cone earch on local HDF5/HTM catalog
            % Package: @catsHTM
            % Description: Perform a cone search around RA/Dec on a local catalog in
            %              HDF5 format sorted into HTM.
            % Input  : - Catalog name (e.g., 'GAIADR1').
            %            see VO.search.htmcat_names for options.
            %          - J2000.0 R.A. [radians, [H M S], or sexagesimal string].
            %          - J2000.0 Dec. [radians, [sign D M S], or sexagesimal string].
            %          - Search radius [arcsec].
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'RadiusUnits' - Radius units. Default is 'arcsec'.
            %            'IndexFileTemplate' - Index Catalog name template.
            %                            Default is '%s_htm.hdf5'.
            %            'CatFileTemplate' - Catalog name template.
            %                            Default is '%s_htm_%06d.hdf5'.
            %            'htmTemplate' - HTM dataset template name.
            %                            Default is 'htm_%06d'.
            %            'NcatInFile'  - Number of Datasets in file.
            %                            Default is 100.
            %            'IndexVarName' - Default is [].
            %            'ColRA'       - Default is 1.
            %            'ColDec'      - Default is2.
            %            'OnlyCone'    - Return only sources within cone.
            %                            If false will return also some
            %                            objects outside cone.
            %                            Default is true.
            %            'ColCellFile' - Default is '%s_htmColCell.mat'.
            %            'OutType'     - Output type {'mat'|'astcat'}.
            %                            Default is 'mat'.
            % Output : - Catalog of source within cone.
            %          - Cell array of column names.
            % License: GNU general public license version 3
            %     By : Eran O. Ofek                    Dec 2017
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Cat=catsHTM.cone_search('UCAC4',1,1,10);
            %          Cat=catsHTM.cone_search('GAIADR1',1,1,10);
            %          Cat=catsHTM.cone_search('GALEX',1,1,10);
            % Reliable: 2
            %--------------------------------------------------------------------------

            RAD = 180./pi;

            DefV.RadiusUnits          = 'arcsec';  % do not change this default!
            DefV.IndexFileTemplate    = '%s_htm.hdf5';
            DefV.CatFileTemplate      = '%s_htm_%06d.hdf5';
            DefV.htmTemplate          = 'htm_%06d';
            DefV.NcatInFile           = 100;
            DefV.IndexVarName         = [];
            DefV.ColRA                = 1;
            DefV.ColDec               = 2;
            DefV.OnlyCone             = true;
            DefV.ColCellFile          = '%s_htmColCell.mat';
            DefV.OutType              = 'mat';

            if (isempty(varargin))
                InPar  = DefV;
                Radius = Radius./(RAD.*3600);  % arcsec to [radians]
            else
                InPar  = InArg.populate_keyval(DefV,varargin,mfilename);
                Radius = convert.angular(InPar.RadiusUnits,'rad',Radius);  % [radians]
            end

            if (ischar(RA))
                RA = celestial.coo.convertdms(RA,'SH','r');
            end
            if (ischar(Dec))
                Dec = celestial.coo.convertdms(Dec,'SD','R');
            end
            
            InPar.ColCellFile = sprintf(InPar.ColCellFile,CatName);

            %Ncol = 45;
            %Ncol = 8;
            load(InPar.ColCellFile);
            Ncol  = numel(ColCell);


            MinDec = Dec - Radius;
            MaxDec = Dec + Radius;

            IndexFileName = sprintf(InPar.IndexFileTemplate,CatName);
            ID     = catsHTM.search_htm_ind(IndexFileName,InPar.IndexVarName,RA,Dec,Radius);
            FileID = floor(ID./InPar.NcatInFile).*InPar.NcatInFile;
            Nid = numel(ID);
            Cat = zeros(0,Ncol);
            %C = Util.struct.struct_def({'Cat'},Nid,1);
            for Iid=1:1:Nid

                %FileID    = floor(ID(Iid)./InPar.NcatInFile).*InPar.NcatInFile;
                FileName  = sprintf(InPar.CatFileTemplate,CatName,FileID(Iid));
                DataName  = sprintf(InPar.htmTemplate,ID(Iid));

                %Cat = [Cat; HDF5.load_cat(FileName,DataName,[MinDec, MaxDec],Ncol)];
                if (Iid==1)
                    Cat = HDF5.load(FileName,DataName);
                    %Cat = HDF5.load_cat(FileName,DataName,[MinDec, MaxDec],Ncol);
                    
                    %Ncol = size(Cat,2);
                else
                    Cat = [Cat; HDF5.load(FileName,DataName)];
                    %Cat = [Cat; HDF5.load_cat(FileName,DataName,[MinDec, MaxDec],Ncol)];
                end

                %C(Iid).Cat = HDF5.load_cat(FileName,DataName,[MinDec, MaxDec],Ncol).';

            end

            %Cat = [C.Cat]';

            % select only sources in Cone
            if (InPar.OnlyCone)
                D = celestial.coo.sphere_dist_fast(RA,Dec,Cat(:,InPar.ColRA),Cat(:,InPar.ColDec));
                Cat = Cat(D<Radius,:);
            end



            switch lower(InPar.OutType)
                case 'mat'
                    % do nothing
                case 'astcat'
                    AstC = AstCat;
                    AstC.Cat = Cat;
                    AstC.ColCell = ColCell;
                    AstC = colcell2col(AstC);
                    Cat  = AstC;
                otherwise
                    error('Unknown OutType option');
            end
            
        end

        function CatM=sources_match(CatName,Cat,varargin)
            % Match sources in an input catalog with catsHTM catalog
            % Package: @catsHTM
            % Description: Given a catalog of sources with their RA/Dec,
            %              match each one of them to a source in an
            %              catsHTM catalog.
            % Input  : - catsHTM catalog name (e.g., 'UCAC4').
            %          - An AstCat object with sources.
            %          * Arbitrary number of key,val pairs:
            %            'OutType' - Output catalog type {'mat'|'astcat'}.
            %                       Default is 'AstCat'.
            %            'SearchRadius' - Search radius. Default is 2.
            %            'SearchRadiusUnits' - Search radius units.
            %                       Default is 'arcsec'.
            %            'ColCell' - Default is {}.
            %            'ColRA' - Default is {'RA','ALPHAWIN_J2000'}.
            %            'ColDec' - Default is {'Dec','DELTAWIN_J2000'}.
            %            'CooUnits' - Input catalog coordinates units.
            %                       Default is 'rad'.
            %            'ColDecHTM' - Default is 2.
            %            'ColRAHTM'  - Default is 1.
            % Output : -
            % Example: CatM=catsHTM.sources_match('GAIADR2',CoaddSim);
            
            CatField     = AstCat.CatField;
            ColCellField = AstCat.ColCellField;
            
            DefV.OutType               = 'AstCat';
            DefV.SearchRadius          = 2;
            DefV.SearchRadiusUnits     = 'arcsec';
            DefV.ColCell               = {};
            DefV.ColRA                 = {'RA','ALPHAWIN_J2000'};
            DefV.ColDec                = {'Dec','DELTAWIN_J2000'};
            DefV.CooUnits              = 'rad';  % in the AstCat object
            DefV.ColDecHTM             = 2;
            DefV.ColRAHTM              = 1;
            
            InPar  = InArg.populate_keyval(DefV,varargin,mfilename);

            InPar.SearchRadius = convert.angular(InPar.SearchRadiusUnits,'rad',InPar.SearchRadius);  % [rad]
            
            % Convert input catalog to an AstCat object
            if (~AstCat.isastcat(Cat))
                Tmp = Cat;
                Cat = AstCat;
                Cat.(CatField)     = Tmp;
                Cat.(ColCellField) = InPar.ColCell;
                Cat                = colcell2col(Cat);
            end
            % RA/Dec columns
            [~,Col.RA,~]     = select_exist_colnames(Cat,InPar.ColRA(:));
            [~,Col.Dec,~]    = select_exist_colnames(Cat,InPar.ColDec(:));
            
            RA  = Cat.(CatField)(:,Col.RA);
            Dec = Cat.(CatField)(:,Col.Dec);
            % convert to radians;
            ConvCoo = convert.angular(InPar.CooUnits,'rad');
            RA      = RA.*ConvCoo;
            Dec     = Dec.*ConvCoo;
            
            MedRA   = nanmedian(RA);
            MedDec  = nanmedian(Dec);
            D       = celestial.coo.sphere_dist_fast(MedRA,MedDec,RA,Dec);
            Radius  = max(D).*(1+10.*eps);  % [rad]
            Radius  = convert.angular('rad','arcsec',Radius); % [arcsec]
            
            [CatH,ColCellH] = catsHTM.cone_search(CatName,MedRA,MedDec,Radius);
            
            CatH = sortrows(CatH,InPar.ColDecHTM);
            
            Nsrc = size(Cat.(CatField),1);
            CatM.Match  = nan(Nsrc,numel(ColCellH));
            CatM.Dist   = nan(Nsrc,1);
            CatM.Nmatch = zeros(Nsrc,1);
            for Isrc=1:1:Nsrc
                % search match for Cat.Cat(Isrc,:)
                Ind = VO.search.search_sortedlat(CatH,RA(Isrc),Dec(Isrc),InPar.SearchRadius);
            
                if (~isempty(Ind))
                    Dist = celestial.coo.sphere_dist_fast(RA(Isrc),Dec(Isrc),CatH(Ind,InPar.ColRAHTM),CatH(Ind,InPar.ColDecHTM));
                    Nmatch = numel(Ind);
                    if (Nmatch>1)
                        [Dist,MinInd] = min(Dist);
                        Ind = Ind(MinInd);
                    end
                        
                    CatM.Match(Isrc,:) = CatH(Ind,:);
                    CatM.Dist(Isrc)    = Dist;
                    CatM.Nmatch(Isrc)  = Nmatch;
                end                
            end
            CatM.ColCell = ColCellH;
            
            switch lower(InPar.OutType)
                case 'astcat'
                    Cat = AstCat;
                    Cat.(CatField) = CatM.Match;
                    Cat.(ColCellField) = CatM.ColCell;
                    Cat = col_insert(Cat,[CatM.Nmatch],numel(CatM.ColCell)+1,'Nmatch');
                    Cat = col_insert(Cat,[CatM.Dist],  numel(CatM.ColCell)+2,'Dist');
                    CatM = Cat;
                case 'mat'
                    % do nothing
                otherwise
                    error('Unknown OutType option');
            end
            
        end
        
        function [ColCell,ConcatRes]=serial_search(CatName,Fun,varargin)
            % Execute a function on entire HDF5/HTM catalog
            % Package: @catsHTM
            % Description: Execute a function on entire HDF5/HTM catalog.
            %              This can be used for selection of sources based
            %              on any parameters.
            % Input  : - Catalog name (e.g., 'GAIADR1').
            %          - Function name to execute:
            %            Fun(Cat,FunPar{:})
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'Concat' - Concat results to previous results.
            %                       Default is true.
            %                       Concat result will be outputed as
            %                       second output argument.
            %            'FunPar' - Cell array of additional parameters to
            %                       pass to the function.
            %            'NparPool' - Number of parallel processes to run.
            %                       Default is 24.
            %            'Xmatch'   - A logical flag indicating if to
            %                       prepare a list of all sources in the
            %                       current HTM and neighbooring HTMs.
            %                       This will be used by FunX.
            %                       Default is false.
            %            'FunX'     - A function to call if Xmatch is true.
            %                       Default is [].
            %                       FunX(Cat,CatNeigh,varargin)
            %            'FunXPar'  - A cell array of key,val atruments to
            %                       pass to FunX as additional parameters.
            %            'SearchRadius' - Search radius for FunX. Default
            %                       is 100 arcsec.
            %            'SearchRadiusUnits' - Default is 'arcsec'.
            %            'Verbose' - Default is true.
            % Output : - Cell array of column names in catalog.
            %          - Optional concat results.
            % Example: [ColCell]=catsHTM.serial_search('APASS',@sin)
            % Reliable: 2
            
            DefV.Concat                = true;
            DefV.FunPar                = {};
            DefV.NparPool              = 24;
            DefV.Xmatch                = false;
            DefV.FunX                  = [];
            DefV.FunXPar               = {};
            DefV.SearchRadius         = 100;  % [arcsec]
            DefV.SearchRadiusUnits    = 'arcsec';
            DefV.ColDec               = 2;          
            DefV.Verbose              = true;
            InPar  = InArg.populate_keyval(DefV,varargin,mfilename);
            
            SearchRadius = convert.angular(InPar.SearchRadiusUnits,'rad',InPar.SearchRadius);  % [rad]

            % load HTM data for Cat1
            [IndexFileName,IndexVarName] = catsHTM.get_index_filename(CatName);
            % HTM1 is the HTM index file
            [HTM,DataHTM] = catsHTM.load_htm_ind(IndexFileName,IndexVarName);
            % Level, Father, Son1, Son2, Son3, Son4, Poles 1 long,
            % poles 1 lat, ..., Nsrc

            %
            Nhtm = numel(HTM);
            
            L = celestial.htm.nhtm2level(Nhtm);
            
            [HTM,Level] = celestial.htm.htm_build(L);
            Level = Level(L);
            Nh    = numel(Level.ptr);
            
            [ColCell] = catsHTM.load_colcell(CatName);
            Ncol      = numel(ColCell);
            
            % number of parallel processes
            %parpool(InPar.NparPool);
            
            %parfor Ih=1:1:Nh
            %Sum{1}=0;
            %Sum{2}=0;
            %tic;
            First = true;
             for Ih=1:1:Nh
                 %Ih
%                  if (Ih./1000)==floor(Ih./1000)
%                     Ih
%                     toc
%                     tic;
%                  end
                % for each HTM in Cat1
                %Cat1    = [];
                Ihtm   = Level.ptr(Ih);
                
                % if HTM in Cat1 contain sources
                if (DataHTM(Ihtm,13)>0)
                    % load Cat
      
                    Cat = catsHTM.load_cat(CatName,Ihtm);
                    
                   
                    if (~isempty(Fun))
                        if (InPar.Concat)
                            CR = Fun(Cat,InPar.FunPar{:});
                            if (InPar.Verbose)
                                fprintf('HTM index: %d    Number of objects: %d\n',Ih,size(CR,1));
                            end
                            if ~isempty(CR)
                                if (First)
                                    First = false;
                                    ConcatRes(1).Cat = CR.';
                                else
                                    ConcatRes(end+1).Cat = CR.';
                                end
                            end
                        else
                            Fun(Cat,InPar.FunPar{:});
                        end
                    end
                    
                end
            end
            
        end
        
        
        function [ColCell]=serial_search_x(CatName,Fun,varargin)
            % Execute a function on entire HDF5/HTM catalog
            % Package: @catsHTM
            % Description: Execute a function on entire HDF5/HTM catalog.
            %              This can be used for selection of sources based
            %              on any parameters.
            % Input  : - Catalog name (e.g., 'GAIADR1').
            %          - Function name to execute:
            %            Fun(Cat,FunPar{:})
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'Istart' - Default is 1.
            %            'Iend'   - Default is inf.
            %            'FunPar' - Cell array of additional parameters to
            %                       pass to the function.
            %            'NparPool' - Number of parallel processes to run.
            %                       Default is 24.
            %            'Xmatch'   - A logical flag indicating if to
            %                       prepare a list of all sources in the
            %                       current HTM and neighbooring HTMs.
            %                       This will be used by FunX.
            %                       Default is false.
            %            'FunX'     - A function to call if Xmatch is true.
            %                       Default is [].
            %                       FunX(Cat,CatNeigh,varargin)
            %            'FunXPar'  - A cell array of key,val atruments to
            %                       pass to FunX as additional parameters.
            %            'SearchRadius' - Search radius for FunX. Default
            %                       is 100 arcsec.
            %            'SearchRadiusUnits' - Default is 'arcsec'.
            % Output : - Cell array of column names in catalog.
            % Example: catsHTM.serial_search_x('GAIADR2',[],'FunX',@search_allml,'Xmatch',true)
            % Reliable: 2
            
            DefV.Istart                = 1;
            DefV.Iend                  = Inf;
            DefV.FunPar                = {};
            DefV.NparPool              = 24;
            DefV.Xmatch                = false;
            DefV.FunX                  = [];
            DefV.FunXPar               = {};
            DefV.SearchRadius         = 100;  % [arcsec]
            DefV.SearchRadiusUnits    = 'arcsec';
            DefV.ColDec               = 2;          
            InPar  = InArg.populate_keyval(DefV,varargin,mfilename);
            
            SearchRadius = convert.angular(InPar.SearchRadiusUnits,'rad',InPar.SearchRadius);  % [rad]

            % load HTM data for Cat1
            [IndexFileName,IndexVarName] = catsHTM.get_index_filename(CatName);
            % HTM1 is the HTM index file
            [HTM,DataHTM] = catsHTM.load_htm_ind(IndexFileName,IndexVarName);
            % Level, Father, Son1, Son2, Son3, Son4, Poles 1 long,
            % poles 1 lat, ..., Nsrc

            %
            Nhtm = numel(HTM);
            
            L = celestial.htm.nhtm2level(Nhtm);
            
            [HTM,Level] = celestial.htm.htm_build(L);
            Level = Level(L);
            Nh    = numel(Level.ptr);
            
            [ColCell] = catsHTM.load_colcell(CatName);
            Ncol      = numel(ColCell);
            
            % number of parallel processes
            %parpool(InPar.NparPool);
            
            %parfor Ih=1:1:Nh
            %Sum{1}=0;
            %Sum{2}=0;
            %ResML = [];
            
            if (isinf(InPar.Iend))
                % do nothing - use Nh
            else
                Nh = InPar.Iend;
            end
            
            
            tic;
             for Ih=InPar.Istart:1:Nh
                 %Ih
                 if (Ih./1000)==floor(Ih./1000)
                    Ih
                    toc
                    tic;
                 end
                % for each HTM in Cat1
                %Cat1    = [];
                Ihtm   = Level.ptr(Ih);
                
                % if HTM in Cat1 contain sources
                if (DataHTM(Ihtm,13)>0)
                    % load Cat
      
                    Cat = catsHTM.load_cat(CatName,Ihtm);
                    
                    if (InPar.Xmatch)
                        
                        % search for all HTMs in Cat2 that may opverlap with
                        % Cat1 current HTM
                        MeanRA  = mean(HTM(Ihtm).coo(:,1));
                        MeanDec = mean(HTM(Ihtm).coo(:,2));
                        MinDec  = min(HTM(Ihtm).coo(:,2))-SearchRadius;
                        MaxDec  = max(HTM(Ihtm).coo(:,2))+SearchRadius;

                        D = celestial.coo.sphere_dist_fast(MeanRA,MeanDec,HTM(Ihtm).coo(:,1),HTM(Ihtm).coo(:,2));
                        CircRadius = max(D) + SearchRadius; % [rad]

                        ID2 = celestial.htm.htm_search_cone(HTM,MeanRA,MeanDec,CircRadius);

                        % load all ID2 from HTM2
                        Nid2 = numel(ID2);
                        for Iid2=1:1:Nid2
                            if (Iid2==1)
                                [Cat2,Ind2]   = catsHTM.load_cat(CatName,ID2(Iid2),[MinDec MaxDec],Ncol);
                                N2     = size(Cat2,1);
                                Cat2ID = [ID2(Iid2).*ones(N2,1), Ind2-1+(1:1:N2)'];
                            else
                                [Cat2tmp, Ind2] = catsHTM.load_cat(CatName,ID2(Iid2),[MinDec MaxDec],Ncol);
                                Cat2   = [Cat2; Cat2tmp];
                                N2     = size(Cat2,1);
                                Cat2ID = [Cat2ID; [ID2(Iid2).*ones(N2,1), Ind2-1+(1:1:N2)']];
                            end
                        end
   
                        if (isempty(Cat2))
                            % if Cat2 is empty - skip
                        else

                            % sort Cat2 and Cat2ID
                            [Cat2,SI] = sortrows(Cat2,InPar.ColDec);
                            %Cat2ID = Cat2ID(SI,:);

                            % cross match Cat1 and Cat2
                            %[Match,Ind,IndCatMinDist] = VO.search.match_cats(Cat2,Cat1,'Radius',SearchRadius,'RadiusUnits','rad');

                            InPar.FunX(Cat,Cat2,ColCell,InPar.FunXPar{:});
                            
                        end
                        
                    end
                    
                    %
                    if (~isempty(Fun))
                        Fun(Cat,ColCell,InPar.FunPar{:});
                    end
                    
                end
             end   % end of for=Ih loop
             %save Sum.mat Sum
             %save ResML.mat ResML
             
        end
    end
    
    % cross matching
    methods (Static)
        function xmatch_2cats(CatName1,CatName2,varargin)
            % Cross match two HDF5/HTM catalogs
            % Package: @catsHTM
            % Description: Cross match two HDF5/HTM catalogs. For each
            %              source in the first catalog the index of the
            %              nearest source, within some distance, in the
            %              second catalog is saved.
            % Input  : - Catalog base name.
            %          - Catalog base name.
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'SearchRadius' - Search radius. Default is 2''.
            %            'SearchRadiusUnits' - Search radius units.
            %                       Default is 'arcsec'.
            %            'QueryFun'    - Optional function handle to
            %                       execute on the matched catalog.
            %                       Syntax:
            %                       Flag = QueryFun(Cat1,Cat2matched,QueryFunPar{:});
            %                       where Flag is a vector logicals
            %                       indicating the selected rows to be
            %                       saved.
            %            'QueryFunPar' - Cell array of additional arguments
            %                       to pass to QueryFun.
            %                       Default is {}.
            %            'SaveFun' - Optional function handle to execute on
            %                       the queried matched catalog.
            %                       SaveFun(Cat1,Cat2matched,SaveFunPar{:});
            %                       E.g., to save the data.
            %            'SaveFunPar' - Cell array of additional arguments
            %                       to pass to SaveFun.
            %                       Default is {}.
            %            'Cat2_ColDec' - Declination column in second
            %                       catalog. Default is 2.
            %            'NparPool' - Number of parallel processes to run.
            %                       Default is 24.
            %            'DeleteParPool' - Delete existing parpool.
            %                       Default is false.
            % Output : null [output is written as an HDF5/HTM catalog].
            % Example: catsHTM.xmatch_2cats('APASS','APASS')
            % Reliable: 2
            
            DefV.SearchRadius         = 2;  % [arcsec]
            DefV.SearchRadiusUnits    = 'arcsec';
            DefV.QueryAllFun          = [];
            DefV.QueryAllFunPar       = {};
            DefV.QueryFun             = [];
            DefV.QueryFunPar          = {};
            DefV.SaveFun              = [];
            DefV.SaveFunPar           = {};
            DefV.Cat2_ColDec          = 2;
            DefV.NparPool             = 24;
            DefV.DeleteParPool        = false;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);
            
            SearchRadius = convert.angular(InPar.SearchRadiusUnits,'rad',InPar.SearchRadius);  % [rad]
            
            % load HTM data for Cat1
            [IndexFileName1,IndexVarName1] = catsHTM.get_index_filename(CatName1);
            % HTM1 is the HTM index file
            [HTM1,DataHTM1] = catsHTM.load_htm_ind(IndexFileName1,IndexVarName1);
            % Level, Father, Son1, Son2, Son3, Son4, Poles 1 long,
            % poles 1 lat, ..., Nsrc

            % load HTM data for Cat2
            [IndexFileName2,IndexVarName2] = catsHTM.get_index_filename(CatName2);
            % HTM2 is the HTM index file
            [HTM2,DataHTM2] = catsHTM.load_htm_ind(IndexFileName2,IndexVarName2);
            % Level, Father, Son1, Son2, Son3, Son4, Poles 1 long,
            % poles 1 lat, ..., Nsrc

            %
            Nhtm1 = numel(HTM1);
            Nhtm2 = numel(HTM2);
            
            L1 = celestial.htm.nhtm2level(Nhtm1);
            L2 = celestial.htm.nhtm2level(Nhtm2);
            LMax = max(L1,L2);
            
            [HTM,Level] = celestial.htm.htm_build(LMax);
            Level1 = Level(L1);
            Level2 = Level(L2);
            Nh1    = numel(Level1.ptr);
            Nh2    = numel(Level2.ptr);
            
            [ColCell1] = catsHTM.load_colcell(CatName1);
            [ColCell2] = catsHTM.load_colcell(CatName2);
            Ncol2      = numel(ColCell2);
            
            % number of parallel processes
            if (InPar.DeleteParPool)
                delete(gcp('nocreate'));
            end
            
            % comment out if needed
            %parpool(InPar.NparPool);
            
            % replace parfor with for if needed
            Nh1
            
            Istart=1;
            for Ih1=Istart:1:Nh1
                Ih1
                % for each HTM in Cat1
                Cat1    = [];
                Cat2    = [];
                Ihtm1   = Level1.ptr(Ih1);
                
                % if HTM in Cat1 contain sources
                if (DataHTM1(Ihtm1,13)>0)
                    % load Cat1
      
                    Cat1 = catsHTM.load_cat(CatName1,Ihtm1);

                    %[Cat2,EdgeOK] = catsHTM.load_cat_edge(CatName2,Ihtm1);
                    %if (~EdgeOk)
                    
                    % search for all HTMs in Cat2 that may opverlap with
                    % Cat1 current HTM
                    MeanRA  = mean(HTM(Ihtm1).coo(:,1));
                    MeanDec = mean(HTM(Ihtm1).coo(:,2));
                    MinDec  = min(HTM(Ihtm1).coo(:,2))-SearchRadius;
                    MaxDec  = max(HTM(Ihtm1).coo(:,2))+SearchRadius;

                    %%
                    %if ((MeanDec.*180./pi)>-30)
                    
                    D = celestial.coo.sphere_dist_fast(MeanRA,MeanDec,HTM(Ihtm1).coo(:,1),HTM(Ihtm1).coo(:,2));
                    CircRadius = max(D) + SearchRadius; % [rad]

                    ID2 = celestial.htm.htm_search_cone(HTM2,MeanRA,MeanDec,CircRadius);
      
                    % load all ID2 from HTM2
                    Nid2 = numel(ID2);
                    for Iid2=1:1:Nid2
                        if (Iid2==1)
                            [Cat2,Ind2]   = catsHTM.load_cat(CatName2,ID2(Iid2),[MinDec MaxDec],Ncol2);
                            N2     = size(Cat2,1);
                            Cat2ID = [ID2(Iid2).*ones(N2,1), Ind2-1+(1:1:N2)'];
                        else
                            [Cat2tmp, Ind2] = catsHTM.load_cat(CatName2,ID2(Iid2),[MinDec MaxDec],Ncol2);
                            Cat2   = [Cat2; Cat2tmp];
                            N2     = size(Cat2,1);
                            Cat2ID = [Cat2ID; [ID2(Iid2).*ones(N2,1), Ind2-1+(1:1:N2)']];
                        end
                    end
   
                    if (isempty(Cat2))
                        % if Cat2 is empty - skip
                    else
                        
                        % sort Cat2 and Cat2ID
                        [Cat2,SI] = sortrows(Cat2,InPar.Cat2_ColDec);
                        %Cat2ID = Cat2ID(SI,:);

                        % cross match Cat1 and Cat2
                        % return list of size of Cat1
                        [Match,Ind,IndCatMinDist] = VO.search.match_cats(Cat2,Cat1,'Radius',SearchRadius,'RadiusUnits','rad');

                        if (~isempty(InPar.QueryAllFun))
                            % execute InPar.QueryAllFun
                            %  QueryAllFun(Cat1,Ind,Cat2,varargin)
                            if (Ih1==Istart)
                                Data = [];
                            end
                            
                            
                            Data = InPar.QueryAllFun(Cat1,Ind,Cat2,IndCatMinDist,InPar.QueryAllFunPar{:},'Data',Data,'Ih1',Ih1,'Nh1',Nh1,'SearchRadius',InPar.SearchRadius);
                        end
                        
                        %Cat2(IndCatMinDist,:)
                        IsN = isnan(IndCatMinDist);
                        IndCatMinDist(IsN) = 1;

                        %DataInd = Cat2ID(IndCatMinDist,:);
                        DataInd = Cat2ID(SI(IndCatMinDist),:);
                        %DataInd(1:2,:)
                        DataInd(IsN,:) = NaN;

                        Cat2matched        = Cat2(IndCatMinDist,:);
                        Cat2matched(IsN,:) = NaN;
                        if (~isempty(InPar.QueryFun))
                            % execute InPar.QueryFun
                            % QueryFun can select specific sources (by some
                            % attributes) from the matched Cat1 and Cat2

                            FlagSelected       = InPar.QueryFun(Cat1,Cat2matched,InPar.QueryFunPar{:});
                            % what to do with FlagSelected?
                            Cat1        = Cat1(FlagSelected,:);
                            Cat2matched = Cat2matched(FlagSelected,:);

                        end

                        if (~isempty(InPar.SaveFun))
                            % execute InPar.SaveFun
                            % Fun(Cat1,Cat2matched)
                            InPar.SaveFun(Cat1,Cat2matched,InPar.SaveFunPar{:});
                        end
                    end
                    %%
                    %end
                end
            end
            
            
            save(sprintf('Data_%d.mat',Ih1),'Data');
     
            
        end
        
        
       
%         function xmatch_save_index(Cat1,Cat2matched,Cat1ID,Cat2matchedID,CatBaseName)
%             % 
%           
%         end
    end
    
    % plots
    methods (Static)
        function H=plot_density(CatName,varargin)
            % Plot a catsHTM catalog surface density
            % Package: @catsHTM
            % Description: Plot a catsHTM catalog surface density in
            %              sources per deg^2 or sources per HTM on a
            %              celestial sphere map.
            % Input  : - Catalog name (e.g., 'NVSS');
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'PerDeg2'  - plot density per deg^2.
            %                         Default is true.
            %                         Otherwise plot per HTM.
            %            'Step'     - Interpolation step in deg.
            %                         Default is 0.3 deg.
            %            'PlotType' - Options are: 'trisurf'|'scatterm'
            %                         Default is 'scatterm'
            %            'MarkerSize'- Marker size for scatterm.
            %                         Default is 5.
            %            'Projection'- Map projection.
            %                         Default is 'aitoff'.
            %            'LogN'     - Plot log10 number of sources.
            %                         Default is false.
            % Output : - Plot handle.
            % Example: H=catsHTM.plot_density('SDSSDR10')
            
            
            RAD = 180./pi;
            
            DefV.PerDeg2              = true;  % otherwise per HTM
            DefV.Step                 = 0.3;  % [deg]
            DefV.PlotType             = 'scatterm';
            DefV.MarkerSize           = 5;
            DefV.Projection           = 'aitoff';
            DefV.LogN                 = false;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            Col.Level    = 1;
            Col.PolesLon = [7 9 11];
            Col.PolesLat = [8 10 12];
            Col.Nsrc     = 13;
            [IndexFileName,IndexVarName]=catsHTM.get_index_filename(CatName);
            [HTM,Data]=catsHTM.load_htm_ind(IndexFileName,IndexVarName);
            F = Data(:,Col.Level) == max(Data(:,Col.Level));
            Data1 = Data(F,:);
            
            Level = celestial.htm.nhtm2level(size(Data,1));
            [HTM,LevelList] = celestial.htm.htm_build(Level);
            Nhtm = numel(LevelList(end).ptr);
            for Ihtm=1:1:Nhtm
                IndHTM  = LevelList(end).ptr(Ihtm);
                MeanRA  = mean(HTM(IndHTM).coo(:,1));
                MeanDec = mean(HTM(IndHTM).coo(:,2));
                Table(Ihtm,:) = [MeanRA, MeanDec, Data1(Ihtm,Col.Nsrc)];
            end
            if (InPar.PerDeg2)
                % Area of HTM triangle [deg^2]
                Area = 4.*pi.*RAD.^2./Nhtm;
                Table(:,3) = Table(:,3)./Area;  % convert to sources per deg^2
            end
            Table(:,1:2) = Table(:,1:2).*RAD;
            
            if (InPar.LogN)
                Table(:,3) = log10(Table(:,3));
            end
            
            switch InPar.PlotType
                case 'trisurf'
                    Tri = delaunay(Table(:,1),Table(:,2));
                    H=trisurf(Tri,Table(:,1), Table(:,2), Table(:,3));
                    view(0,90);
                    shading interp
                    colorbar
                case 'scatterm'
                    axesm(InPar.Projection); 
                    framem
                    H=scatterm(Table(:,2),Table(:,1),InPar.MarkerSize,Table(:,3),'filled');
                    colorbar
                otherwise
                    error('Unknown PlotType option');
            end
            
%             F = scatteredInterpolant(Table(:,2).*RAD,Table(:,1).*RAD,Table(:,3));
%             Lon = (-180:InPar.Step:180);
%             Lat = (-90:InPar.Step:90);
%             [MLon,MLat] = meshgrid(Lon,Lat);
%              
%             F.Method = 'nearest';
%             Ninterp = F(MLat,MLon);
%             surface(Lat,Lon,Ninterp');
%             shading interp
%             colorbar
%             
   

        end
    end
end % end class
            
