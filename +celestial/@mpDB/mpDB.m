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
%--------------------------------------------------------------------------

%     properties
%         mpDB_dir = ''; 
%         
%         
%     end
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


    % mpDB
    methods (Static)
        function Dir=getenvDir
            % get AsteroidsEphemDB files directory from env. var.
            % Package: +celestial/@mpDB
            % Output : - String of directory
            
            Dir = getenv('mpDB_dir');
        end % end getenvDir function
        
        function setenvDir(Dir)
            % set AsteroidsEphemDB files directory in env. var. 'mpDB_dir'
            % Package: +celestial/@mpDB
            % Input : - String of directory
            
            setenv('mpDB_dir',Dir);
        end % end setenvDir function
        
        function [FileName,DataSet,Data]=loadEphemFile(ID,AstType,Type,Year)
            % Construct Ephem file name and load data
            % Package: +celestial/@mpDB
            % Description:
            % Input  : - ID [StartYear, EndYear, StartObj, EndObj]
            %          - AstType: 'n', 'u', 'c'. Default is 'n'.
            %            n - numbered asteroids, u- unnumbered, c-comets
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
            % [FileName,DataSet,Data]=celestial.mpDB.loadEphemFile([2017 2020 363001 364000],'ephem')
            % Reliable: 2
            
            if (nargin<3)
                Type = 'ephem';
                if nargin<2
                    AstType = 'n';
                end
            end
            
            Dir      = celestial.mpDB.getenvDir;
            % File name example
            % AstEphem_2017_2020_363001_364000.hdf5
            
            FileName = sprintf('%s%sAstEphem_%s_%04d_%04d_%06d_%06d.hdf5',Dir,filesep,AstType,ID);
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
        
        function Cat=createEphemCatalog
            %
            % Example: Cat=celestial.mpDB.createEphemCatalog
            
            Dir      = celestial.mpDB.getenvDir;
            Files    = dir(sprintf('%s%s*.hdf5',Dir,filesep));
            Nf = numel(Files);
            Cat = Util.struct.struct_def({'ID'},Nf,1);
            for If=1:1:Nf
                Cat(If).FileName = Files(If).name;
                
                [~,File,Ext]=fileparts(Files(If).name);
                
                switch lower(Ext)
                    case '.hdf5'
                        Tmp = regexp(File,'_','split');
                        Year1 = str2double(Tmp{2});
                        Year2 = str2double(Tmp{3});
                        Obj1  = str2double(Tmp{4});
                        Obj2  = str2double(Tmp{5});
                        
                        Cat(If).ID   = [Year1; Year2; Obj1; Obj2];
                    otherwise
                        error('Was suppose to select only hdf5 files');
                end
            end
            
        end % end createEphemCatalog function
        
    end % methods
    
    
end % end class
            
 