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
        Observatory = '675';  % Palomar
        
        
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
    end % methods s=(Static)
    
    methods (Static)
        function generateAllAstEphemFiles(varargin)
            % generate all AstEphem files
            % Package: +celestial/@mpDB
            % Description: generate all AstEphem files, including
            %              the ObjReferenceFile.mat file.
            % Example: celestial.mpDB.generateAllAstEphemFiles('ObjType',2)
            
            DefV.ObjType              = [1,2,3];
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            DataObj = celestial.SolarSys.get_orbit_files('load');
            
            F1 = DataObj.ObjType==1;
            F2 = DataObj.ObjType==2;
            F3 = DataObj.ObjType==3;
            
            AstIndex = zeros(size(DataObj.ObjType));
            AstIndex(F1) = (1:1:sum(F1))';
            AstIndex(F2) = (1:1:sum(F2))';
            AstIndex(F3) = (1:1:sum(F3))';
            
            % save objects reference file
            ObjReferenceFile = table(DataObj.Number,DataObj.Name,DataObj.ObjType,AstIndex,true(size(AstIndex)));
            ObjReferenceFile.Properties.VariableNames={'Number','Name','ObjType','ObjIndex','Status'};
            save -v7.3 ObjReferenceFile.mat ObjReferenceFile
            
            Ntype = numel(InPar.ObjType);
            for Itype=1:1:Ntype
                
                Flag     = DataObj.ObjType == InPar.ObjType(Itype);
                
                FileType = celestial.mpDB.astType2index(InPar.ObjType(Itype));
                
                switch lower(FileType)
                    case 'n'
                        % numbered asteroids
                        ObjName = DataObj.Number(Flag);
                    case 'u'
                        % unnumbered asteroids
                        ObjName = DataObj.Name(Flag);

                    case 'c'
                        % comets
                        ObjName = DataObj.Name(Flag);

                    case 'p'
                        % planets
                        error('p option not implemented yet');
                    otherwise
                        error('Unknown FileType option');
                end
               
                celestial.mpDB.generateAstEphemFiles('ObjectNames',ObjName,'AstIndex',AstIndex(Flag),'FileType',FileType);

                
                
            end
            
            
            
        end
        
        function generateAstEphemFiles(varargin)
            % generate AstEphem files with ephmerides and polynomial fits
            % Package: +celestial/@mpDB
            % Description: generate an AstEphem files in HDF5 format.
            %              Each file contains the ephemerides for 'NinFile'
            %              objects (/Cat dataset) and polynomial fits of
            %              positions in each year (/Poly%04d dataset).
            % Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'ObjName' - Object name. Name, or number.
            %                        For numbered asteroids use e.g., '1;'.
            %            'AstIndex'- This is an index number (e.g., object ID)
            %                        that will be added as column
            %                        'AstIndex' to the catalog.
            %                        Default is 1.
            %            'FileType'  - Default is 'n'.
            %            'YearStart' - Default is 2017.
            %            'YearEnd'   - Default is 2021.
            %            'StepSize'  - Step size. Default is 8.
            %            'StepSizeUnits' - Step size units. Default is 'h'.
            %            'CENTER'    - Observer position.
            %                          '500' - Geocentric.
            %                          Default is '675' (i.e., Palomar).
            %            'TimeBuffer' - Additional days before/after
            %                         start/end year.
            %            'PauseAfterError' - Default is 300 s.
            %            'NinFile' - Default is 1000.
            % Example:
            % celestial.mpDB.generateAstEphemFiles('ObjectNames',(1:1:10),'AstIndex',(1:1:10),'NinFile',5);
            % celestial.mpDB.generateAstEphemFiles('ObjectNames',(541001:1:541132),'AstIndex',(541001:1:541132),'NinFile',1000);

            
            
            DefV.ObjectNames          = (1:1:1000);
            DefV.AstIndex             = [];  % must be provided
            DefV.FileType             = 'n';
            
            DefV.YearStart            = 2017;
            DefV.YearEnd              = 2021;
            DefV.StepSize             = 8;
            DefV.StepSizeUnits        = 'h';
            DefV.CENTER               = '675';  % Palomar
            DefV.TimeBuffer           = 3;   % days
            DefV.PauseAfterError      = 300;  % s
            
            
            DefV.NinFile              = 1000;
            
            DefV.BaseFileName         = 'AstEphem_[nucp]_\w+.hdf5';
            DefV.Save                 = 'EphemCatalog.mat';

            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            
            if isempty(InPar.AstIndex)
                error('AstIndex must be provided');
                if (numel(InPar.AstIndex)~=numel(ObjectNames))
                    error('AstIndex and ObjectNames must have the same number of elements');
                end
            end
            
            Nobj = numel(InPar.ObjectNames);
            %Nobj = 541132
            
            K = 0;
            for Iobj=1:1:Nobj
                if K==0
                    Istart = Iobj;
                    AstIndex1 = InPar.AstIndex(Istart);
                    
                    Cat = AstCat;
                    
                end
                AstIndex2 = InPar.AstIndex(Iobj);
                

                K = K + 1;
                if (iscell(InPar.ObjectNames))
                    Object = InPar.ObjectNames{Iobj};
                elseif (isnumeric(InPar.ObjectNames))
                    Object = sprintf('%d;',InPar.ObjectNames(Iobj));
                else
                    error('Unknown ObjectNames format');
                end

                Iobj
                Object

                [Cat(K)]=celestial.mpDB.jpl_generate_1obj_ephem('ObjName',Object,...
                            'AstIndex',InPar.AstIndex(Iobj),...
                            'YearStart',InPar.YearStart,...
                            'YearEnd',InPar.YearEnd',...
                            'StepSize',InPar.StepSize,...
                            'StepSizeUnits',InPar.StepSizeUnits,...
                            'CENTER',InPar.CENTER,...
                            'TimeBuffer',InPar.TimeBuffer,...
                            'PauseAfterError',InPar.PauseAfterError);

                 
                 
                 PolyDB(K,:)=celestial.mpDB.fit_poly2cat(Cat(K),'YearStart',InPar.YearStart,...
                            'YearEnd',InPar.YearEnd',...
                            'TimeBuffer',InPar.TimeBuffer);


                        
                if (Iobj./InPar.NinFile)==floor(Iobj./InPar.NinFile) || Iobj==Nobj

                    SizeCat = sizecat(Cat);
                    Ind = [1, 1+cumsum(SizeCat)];
                    IndS = Ind(1:end-1).';
                    IndE = Ind(2:end).'-1;

                    CatM = merge(Cat);
                    CatM.Cat = real(CatM.Cat);

                    FileName = sprintf('AstEphem_%s_%04d_%04d_%06d_%06d.hdf5',InPar.FileType,InPar.YearStart,InPar.YearEnd, AstIndex1  , AstIndex2);
                    delete(FileName);
                    VarName  = sprintf('/Ephem');
                    Attrib = [Cat(1).ColCell', num2cell(1:1:numel(Cat(1).ColCell))'];

                    try
                        HDF5.save(CatM.Cat,FileName,VarName,Attrib);

                        VarName  = sprintf('/Ind');
                        Data = [(Is:1:I).', [IndS, IndE]];
                        HDF5.save(Data,FileName,VarName,{});

                        for Iyears=1:1:Nyears
                            VarName  = sprintf('/Poly%04d',Years(Iyears));
                            DataX = [PolyDB(:,Iyears).ParCX].';
                            DataY = [PolyDB(:,Iyears).ParCY].';
                            DataZ = [PolyDB(:,Iyears).ParCZ].';
                            DataD = [PolyDB(:,Iyears).MaxD].';
                            Data  = [DataX, DataY, DataZ, DataD];

                            HDF5.save(Data,FileName,VarName,{});
                        end
                    catch
                        % do nothing
                    end

                    clear Cat;
                    K = 0;
                end

                  
            end
                
            
            
        end
        
        function [Cat]=jpl_generate_1obj_ephem(varargin)
            % Generate JPL ephemerides for one object
            % Package: +celestial/@mpDB
            % Description: Generate JPL ephemerides for one object
            % Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'ObjName' - Object name. Name, or number.
            %                        For numbered asteroids use e.g., '1;'.
            %            'AstIndex'- This is an index number (e.g., object ID)
            %                        that will be added as column
            %                        'AstIndex' to the catalog.
            %                        Default is 1.
            %            'YearStart' - Default is 2017.
            %            'YearEnd'   - Default is 2021.
            %            'StepSize'  - Step size. Default is 8.
            %            'StepSizeUnits' - Step size units. Default is 'h'.
            %            'CENTER'    - Observer position.
            %                          '500' - Geocentric.
            %                          Default is '675' (i.e., Palomar).
            %            'TimeBuffer' - Additional days before/after
            %                         start/end year.
            %            'PauseAfterError' - Default is 300 s.
            % Output : - A catalog with object epehmerides.
            % Example: [Cat]=celestial.mpDB.jpl_generate_1obj_ephem
            
            % Defaults
            DefV.ObjName              = '1;';
            DefV.AstIndex             = 1;
            %DefV.ObjType              = 'n';
            DefV.YearStart            = 2017;
            DefV.YearEnd              = 2021;
            DefV.StepSize             = 8;
            DefV.StepSizeUnits        = 'h';
            DefV.CENTER               = '675';  % Palomar
            DefV.TimeBuffer           = 3;   % days
            DefV.PauseAfterError      = 300;  % s
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);


            YearStart = InPar.YearStart;
            YearEnd   = InPar.YearEnd;
            
            Years     = (YearStart:YearEnd);
            Nyears    = numel(Years);

            TimeBuffer = InPar.TimeBuffer;
            StartJD    = celestial.time.julday([1   1 YearStart])-TimeBuffer;
            StopJD     = celestial.time.julday([31 12 YearEnd])+TimeBuffer;

            Object     = InPar.ObjName;


            K = 1;
            
            try
                tic
                [Cat(K)]=celestial.SolarSys.jpl_horizons('StartJD',StartJD,'StopJD',StopJD,'ObjectInd',Object,'CENTER',InPar.CENTER,'StepSizeUnits',InPar.StepSizeUnits,'StepSize',InPar.StepSize);
                toc
            catch
                pause(InPar.PauseAfterError);
                
                tic
                [Cat(K)]=celestial.SolarSys.jpl_horizons('StartJD',StartJD,'StopJD',StopJD,'ObjectInd',Object,'CENTER',InPar.CENTER,'StepSizeUnits',InPar.StepSizeUnits,'StepSize',InPar.StepSize);
                toc
               
            end

            pause(1);
            
            if isempty(Cat(K).Cat)
               Object
               warning('empty Cat');
               
            else

            

                % add column
                Nep = size(Cat(K).Cat,1);
                Cat(K).Cat = [ones(Nep,1).*InPar.AstIndex, Cat(K).Cat];
                Cat(K).ColCell = ['AstIndex', Cat(K).ColCell];
                Cat(K).ColUnits = ["", Cat(K).ColUnits];
                Cat(K) = colcell2col(Cat(K));

                % select column
                Cat(K) = col_select(Cat(K),{'AstIndex','JD','RA','Dec','APmag','r','rdot','Delta','Deltadot','SubObsTargetAng','TrailLeadSun','SubTargetObsAng'});

                % add angular rate column
                Col = Cat(K).Col;
                D = celestial.coo.sphere_dist_fast(Cat(K).Cat(1:end-1,Col.RA),...
                                               Cat(K).Cat(1:end-1,Col.Dec),...
                                               Cat(K).Cat(2:end,Col.RA),...
                                               Cat(K).Cat(2:end,Col.Dec));
                D = [D; D(end)];
                Cat(K).Cat = [Cat(K).Cat, D];
                Cat(K).ColCell = [Cat(K).ColCell, "RateDay"];
                Cat(K).ColUnits = [Cat(K).ColUnits, "rad/day"];
                Cat(K) = colcell2col(Cat(K));
                Col = Cat(K).Col;
            end

        end
            
        function PolyDB=fit_poly2cat(Cat,varargin)
            % Fit a yearly position polynomials to JPL ephemerides
            % Package: +celestial/@mpDB
            % Description: Given an AstCat object generated by 
            %              celestial.mpDB.jpl_generate_1obj_ephem and
            %              contains ephemerides of single object, fit
            %              yearly polynomials to the cosine directions of
            %              the object.
            % Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'YearStart' - Default is 2017.
            %            'YearEnd'   - Default is 2021.
            %            'TimeBuffer' - Additional days before/after
            %                         start/end year.
            % Output : - A structure array in which each element
            %            corresponds to one year. The following fields are
            %            available:
            %            'MaxD' - maximum deviation between the polynomial
            %                     and ephmerides [radians].
            %            'ParCX' - Polynomials for the cosine-direction in
            %                     X coordinate.
            %                     Time is measured in Julian years relative
            %                     to Jan 1st of the year.
            % Example: [Cat]=celestial.mpDB.jpl_generate_1obj_ephem;
            %          PolyDB=celestial.mpDB.fit_poly2cat(Cat);
            
            

            DefV.YearStart            = 2017;
            DefV.YearEnd              = 2021;
            DefV.TimeBuffer           = 3;   % days
            DefV.PolyOrder            = 9;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);




            YearStart = InPar.YearStart;
            YearEnd   = InPar.YearEnd;

            Years     = (YearStart:YearEnd);
            Nyears    = numel(Years);
            TimeBuffer = InPar.TimeBuffer;
            PolyOrder  = InPar.PolyOrder;
            Col = Cat.Col;

         
            PolyDB = Util.struct.struct_def({'MaxD','ParCX','ParCY','ParCZ'},1,Nyears);

            K = 1;
            if ~isempty(Cat.Cat)
                
                % fit polynomial to asteroid motion per year
                for Iyears=1:1:Nyears
                    RefTime = celestial.time.julday([1 1 Years(Iyears)]);
                    RefYear = 365.25;

                    Ftime = Cat(K).Cat(:,Col.JD)>=(celestial.time.julday([1 1 Years(Iyears)])-TimeBuffer) & ...
                            Cat(K).Cat(:,Col.JD)<=(celestial.time.julday([31 12 Years(Iyears)])+TimeBuffer);
                    [CX,CY,CZ] = celestial.coo.coo2cosined(Cat(K).Cat(Ftime,Col.RA), Cat(K).Cat(Ftime,Col.Dec));

                    Time = (Cat(K).Cat(Ftime,Col.JD) - RefTime)./RefYear;

                    ParCX  = polyfit(Time,CX,PolyOrder);
                    ParCY  = polyfit(Time,CY,PolyOrder);
                    ParCZ  = polyfit(Time,CZ,PolyOrder);

                    PredCX = polyval(ParCX,Time);
                    PredCY = polyval(ParCY,Time);
                    PredCZ = polyval(ParCZ,Time);
                    [PredRA,PredDec] = celestial.coo.cosined2coo(PredCX,PredCY,PredCZ);

                    D = celestial.coo.sphere_dist_fast(Cat(K).Cat(Ftime,Col.RA), Cat(K).Cat(Ftime,Col.Dec),...
                                                       PredRA, PredDec);
                    MaxD = max(D); 
                    PolyDB(K,Iyears).MaxD  = MaxD;  %.*RAD
                    PolyDB(K,Iyears).ParCX = ParCX.';
                    PolyDB(K,Iyears).ParCY = ParCY.';
                    PolyDB(K,Iyears).ParCZ = ParCZ.';

                end
            else
                % empty Cat
                for Iyears=1:1:Nyears
                    PolyDB(K,Iyears).MaxD  = NaN;  %.*RAD
                    PolyDB(K,Iyears).ParCX = nan(InPar.PolyOrder+1,1);
                    PolyDB(K,Iyears).ParCY = nan(InPar.PolyOrder+1,1);
                    PolyDB(K,Iyears).ParCZ = nan(InPar.PolyOrder+1,1);
                end
                
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
        
        function EphemCatalog=generateCatOfAstEphemFiles(varargin)
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
            % Example: Cat=celestial.mpDB.generateCatOfAstEphemFiles
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
                
                %EphemCatalog = celestial.mpDB.createEphemCatalog('Save',[]);
                % load EphemCatalog
                load EphemCatalog.mat;
                
            else
                EphemCatalog = InPar.EphemCatalog;
            end
            
            Nf = numel(EphemCatalog);
            for If=1:1:Nf
                If
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
            
 