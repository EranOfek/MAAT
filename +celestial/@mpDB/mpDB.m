classdef mpDB < handle
% mpDB static class                       
% Package: +celestial/@mpDB
% Description: A static class for local minor planets ephemerides files.
% Instellation: The OS environment variable mpDB_dir contains the directory
%               in which the AsteroidsEphemDB HDF5 files resides
%               (e.g., '/euler/eran/work/AsteroidsEphemDB').
%               We recomend to set the value of this environment variable
%               in the starup.m file (i.e., using setenv.m).
% Input  : null
% Output : null
% Tested : Matlab R2018b
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2


    properties
       PolyYear = []; 
    
        %mpDB_dir = ''; 
        DB       = [];  % Table [ID, Name, AstType, Y1, Y2, Obj1, Obj2]
        PolyDB   = [];
        
        
    end
% 
%     % Constructor
%     methods
%         function H=mpDB(H)
%             % mpDB class constructor
%             % Package: +celestial/@mpDB
%             % Input  : 
%             % Output : - mpDB object
%             
%             H.mpDB_dir = getenv('mpDB_dir');
%             
%         end
%         
%         
%     end % end methods
%     

    % getter/setter
    methods
        function H=get.PolyDB(H)
            %
            
            if isempty(H.PolyYear)
                error('In order to populate the DB first set the PolyYear');
            end
            
            if isempty(H.PolyDB)
                H.PolyDB = celestial.mpDB.load_all_poly(H.PolyYear);
            end
            
        end
        
    end

    % mpDB / aux function
    methods (Static)
        function Dir=getenvDir
            % get AsteroidsEphemDB files directory from env. var.
            % Package: +celestial/@mpDB
            % Output : - String of directory
            % Example: celestial.mpDB.getenvDir
            
            Dir = getenv('mpDB_dir');
        end % end getenvDir function
        
        function setenvDir(Dir)
            % set AsteroidsEphemDB files directory in env. var. 'mpDB_dir'
            % Package: +celestial/@mpDB
            % Input : - String of directory
            %           Default is '/euler/eran/work/AsteroidsEphemDB'
            % Example: celestial.mpDB.setenvDir
            
            if nargin<1
                Dir = '/euler/eran/work/AsteroidsEphemDB';
            end
            
            setenv('mpDB_dir',Dir);
        end % end setenvDir function
        
        function Out=astType2index(In)
            % AstType file 'n','u','c','p', to index and vise versa.
            % Package: +celestial/@mpDB
            % Input  : - A character, or an index.
            % Output : - Index or character.
            % Example: celestial.mpDB.astType2index('c')
            %          celestial.mpDB.astType2index(3)
            
            Str = 'nucp';
            
            if ischar(In)
                Out = strfind(Str,In);
            else
                Out = Str(In);
            end
            
        end
        
        function [FileName,DataSet,Data]=loadEphemFile(ID,Type,Year)
            % Construct Ephem file name and load data
            % Package: +celestial/@mpDB
            % Description:
            % Input  : - ID [AstTypeIndex, StartYear, EndYear, StartObj, EndObj]
            %            where AstTypeIndex is the index corresponding to
            %            'nucp' (see astType2index.m).
            %          - Dataset type:
            %            'poly' - polynomial approximation of RA/Dec
            %            'ephem' - Full epehmerides. (Default).
            %            'ind' - ephemerides index file.
            %          - The year of the polynomial file.
            %            Mandatory if If 'poly' option is selected.
            % Output : - File name.
            %          - Dataset name.
            %          - Data.
            % Example:
            % [FileName,DataSet,Data]=celestial.mpDB.loadEphemFile([1 2017 2020 363001 364000],'ephem');
            % Reliable: 2
            
            if (nargin<2)
                Type = 'ephem';
            end
            
            Dir      = celestial.mpDB.getenvDir;
            % File name example
            % AstEphem_2017_2020_363001_364000.hdf5
            
            AstType = celestial.mpDB.astType2index(ID(1));
            
            FileName = sprintf('%s%sAstEphem_%s_%04d_%04d_%06d_%06d.hdf5',Dir,filesep,AstType,ID(2:end));
            switch lower(Type)
                case 'poly'
                    if nargin<3
                        error('when Type=poly is selected - Year must be provided');
                    end
                    DataSet = sprintf('/Poly%04d',Year);
                case 'ind'
                    DataSet = '/Ind';
                case 'ephem'
                    DataSet = '/Ephem';
                otherwise
                    error('Unknown Type option');
            end
            if (nargout>2)
                Data = h5read(FileName,DataSet);
            end
            
        end % end function loadEphemFile
        
        function EphemCatalog=createEphemCatalog(varargin)
            % Generate a catalog of AstEphem files in the mpDB_dir
            % Package: +celestial/@mpDB
            % Description: Generate a catalog of AstEphem files in the
            %              mpDB_dir directory.
            %              The ID vector contains:
            %              [AstType, YearStart, YearEnd, Nstart, Nend]
            % Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'BaseFileName' - Default is 'AstEphem_[nucp]_\w+.hdf5'.
            %            'Save'         - File name in which to save the
            %                             catalog. If empty, then do  not
            %                             save. Default is
            %                             'EphemCatalog.mat'.
            % Output : - Structure array of all identified files, with the
            %            following fields:
            %            'ID' - The ID vector contains:
            %                   [AstType, YearStart, YearEnd, Nstart, Nend]
            %            'FileName' - File name.
            % Example: Cat=celestial.mpDB.createEphemCatalog
            % Reliable: 2
            
            
            DefV.BaseFileName         = 'AstEphem_[nucp]_\w+.hdf5';
            DefV.Save                 = 'EphemCatalog.mat';
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            
            Dir      = celestial.mpDB.getenvDir;
            Files    = dir(sprintf('%s%s*.hdf5',Dir,filesep));
            Tmp      = regexp({Files.name},InPar.BaseFileName,'match');
            Flag     = ~Util.cell.isempty_cell(Tmp);
            Files    = Files(Flag);
            
            Nf = numel(Files);
            EphemCatalog = Util.struct.struct_def({'ID'},Nf,1);
            for If=1:1:Nf
                EphemCatalog(If).FileName = Files(If).name;
                
                [~,File,Ext]=fileparts(Files(If).name);
                
                switch lower(Ext)
                    case '.hdf5'
                        Tmp = regexp(File,'_','split');
                        AstType = Tmp{2};
                        Year1 = str2double(Tmp{3});
                        Year2 = str2double(Tmp{4});
                        Obj1  = str2double(Tmp{5});
                        Obj2  = str2double(Tmp{6});
                        
                        ATIndex = celestial.mpDB.astType2index(AstType);
                        
                        EphemCatalog(If).ID   = [ATIndex; Year1; Year2; Obj1; Obj2];
                    otherwise
                        error('Was suppose to select only hdf5 files');
                end
            end
            
            if ~isempty(InPar.Save)
                PWD = pwd;
                cd(Dir);
                save(InPar.Save,'EphemCatalog','-v7.3');
                cd(PWD);
            end
            
        end % end createEphemCatalog function
        
        
    end % methods
    
    % mpDB / read functions
    methods (Static)
        function PolyS=load_all_poly(Year,varargin)
            % Load the 10th deg polynomial of MP coordinates for one tear. 
            % Package: +celestial/@mpDB
            % Description: Load the 10th deg polynomial of MP coordinates
            %              for one tear.
            % Input  : - Year for which the polynomials are relevant.
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'EphemCatalog' - EphemCatalog as input, in this
            %                      case the PolyS will not be loaded.
            % Output : - A structure containing the following fields:
            %            'CD1' - A mtaritx of polynomials representing the
            %                    1st cosine direction. Line per asteroid,
            %                    10 columns.
            %            'CD2' - A mtaritx of polynomials representing the
            %                    2nd cosine direction.
            %            'CD3' - A mtaritx of polynomials representing the
            %                    3rd cosine direction.
            %            'MaxDist' - Maximal deviation between the
            %                    polynomial representation and the true
            %                    position of the asteroid. [Radians].
            % Example: PolyS=celestial.mpDB.load_all_poly(2017);
            % Reliable: 2
            
            
            DefV.EphemCatalog         = [];
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            if isempty(InPar.EphemCatalog)
                % Generate EphemCatalog
                EphemCatalog = celestial.mpDB.createEphemCatalog('Save',[]);
            else
                EphemCatalog = InPar.EphemCatalog;
            end
            
            Nf = numel(EphemCatalog);
            for If=1:1:Nf
                [~,~,Data]=celestial.mpDB.loadEphemFile(EphemCatalog(If).ID.','poly',Year);
                if (If==1)
                    AllData = Data;
                else
                    AllData = [AllData; Data];
                end
            end
            
            PolyS.CD1     = AllData(:,1:10);
            PolyS.CD2     = AllData(:,11:20);
            PolyS.CD3     = AllData(:,21:30);
            PolyS.MaxDist = AllData(:,31);
            
        end
        
    end % methods
    
    
end % end class
            
 