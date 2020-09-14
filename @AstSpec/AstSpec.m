% A class for astronomical spectra
% Package: @AstSpec
% Description: A class of structure array of astronomical spectra.
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef AstSpec < HEAD
    properties (SetAccess = public)
        Wave
        Int
        Err
        Back
        Mask
        WaveUnits
        IntUnits
        AddCol
        ObjName
        comments
        source
        FileName
        z
    end
    
    
    % Constructors
    methods

         function AstS=AstSpec(N,M)
             % AstSpec constructor
             % Package: @AstSpec
             % Description: AstSpec constructor method
            
            WaveField = 'Wave';
            
            if (nargin==0)
                N = 1;
                M = 1;
            elseif (nargin==1)
                if (numel(N)>1)
                    M = N(2);
                else
                    M = 1;
                end
            else
                % do nothing
            end

            for I=1:1:N
                for J=1:1:M
                    AstS(I,J).(WaveField) = [];
                end
            end
         end
        
        

    end
    
    %----------------------
    %--- Static classes ---
    %----------------------
    
    % convert classes to AstSpec
    methods (Static)
        
        function Ans=isastspec(Obj)
            %--------------------------------------------------------------------------
            % isastspec function                                             AstroSpec
            % Description: Check if an object is of AstSpec class.
            % Input  : - An object.
            % Output : - true if the object is AstSpec class, otherwise false.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Feb 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: AstSpec.isastspec(Obj)
            % Reliable: 2
            %--------------------------------------------------------------------------

            Ans = isa(Obj,'AstSpec');
         end
        
        function AstS=array2astspec(Array)
            % Convert an array to an AstSpec object
            % Package: @AstSpec
            % Description: Convert a matrix into an AstSpec object.
            % Input  : - An array with 2 or 3 columns [Wave, Int, Err].
            % Output : - An AstSpec object.
            % Example: A=AstSpec.array2astspec(rand(100,3));
            % Reliable: 2
            
            AstS = AstSpec;
            AstS.Wave = Array(:,1);
            AstS.Int  = Array(:,2);
            if (size(Array,2)>2)
                AstS.Err  = Array(:,3);
            end
            
        end
        
        function AstS=mat2spec(Mat,Columns,Units,Header)
            % Convert an array to AstSpec object, with more advance options
            % Description: Convert an array to AstSpec object, with more advance options
            % Input  : - A matrix
            %          - Cell array of column names in the matrix.
            %            Default is {'Wave','Int','Err'}.
            %          - Cell array of units associated with the columns.
            %            Default is
            %            {'Ang','erg*cm^-2*s^-1*Ang^-1','erg*cm^-2*s^-1*Ang^-1'}.
            %          - Header. Default is {}.
            
            AstS = AstSpec;
            
            Ncol = size(Mat,2);
           
            Def.Columns = {'Wave','Int','Err'};
            Def.Units   = {'Ang','erg*cm^-2*s^-1*Ang^-1','erg*cm^-2*s^-1*Ang^-1'};
            Def.Header  = {};
            if (nargin==2)
                Columns   = Def.Columns;
                Units     = Def.Units;
                Header    = Def.Header;
            elseif (nargin==3)
                Units     = Def.Units;
                Header    = Def.Header;
            elseif (nargin==4)
                Header    = Def.Header;
            elseif (nargin==5)
                % do nothing
            else
                error('Illegal number of input arguments');
            end
            % make sure that Units and columns have the correct num of
            % elements
            Units   = Units(1:Ncol);
            Columns = Columns(1:Ncol);
            
            Iwave    = find(strcmp(Columns,'Wave'));
            Iint     = find(strcmp(Columns,'Int'));
            Ierr     = find(strcmp(Columns,'Err'));
            Iback    = find(strcmp(Columns,'Back'));
            Imask    = find(strcmp(Columns,'Mask'));
            
            Nunits = numel(Units);
            if (isempty(Iwave))
                error('Must supply wavelength');
            else
                AstS(1).Wave      = Mat(:,Iwave);
                if (Iwave>Nunits)
                    AstS(1).WaveUnits = '';
                else
                    AstS(1).WaveUnits = Units{Iwave};
                end
            end
            if (isempty(Iint))
                error('Must supply intensity');
            else
                AstS(1).Int    = Mat(:,Iint);
                if (Iint>Nunits)
                    AstS(1).IntUnits = '';
                else
                    AstS(1).IntUnits = Units{Iint};
                end
            end
            if (~isempty(Ierr))
                AstS(1).IntErr = Mat(:,Ierr);
            end
            if (~isempty(Iback))
                AstS(1).Back   = Mat(:,Iback);
            end
            if (~isempty(Imask))
                AstS(1).Mask   = Mat(:,Imask);
            end
            
            AstS(1).WaveUnits  = Units{Iwave};
            AstS(1).IntUnits   = Units{Iint};
            AstS(1).Header = Header;
           
         end
        
        function AstS=mat2astspec(Mat,Columns,Units,Header)
            % Convert a matrix into an AstSpec class object
            % Package: @AstSpec
            % Input  : - Matrix.
            %          - Cell array containing columns description.
            %            Default is {'Wave','Int','Err'}.
            %            The Wave column will populate the spectrum column,
            %            etc.
            %          - Cell array of column units.
            %            Default is
            %            {'Ang','erg*cm^-2*s^-1*Ang^-1','erg*cm^-2*s^-1*Ang^-1'}.
            %          - Optional header. Default is empty.
            % Output : - An AstSpec class object.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: AstS = AstSpec.mat2astspec(rand(100,3));
            % Reliable: 2


            Ncol = size(Mat,2);

            Def.Columns = {'Wave','Int','Err'};
            Def.Units   = {'Ang','erg*cm^-2*s^-1*Ang^-1','erg*cm^-2*s^-1*Ang^-1'};
            Def.Header  = {};
            if (nargin==1)
                Columns   = Def.Columns;
                Units     = Def.Units;
                Header    = Def.Header;
            elseif (nargin==2)
                Units     = Def.Units;
                Header    = Def.Header;
            elseif (nargin==3)
                Header    = Def.Header;
            elseif (nargin==4)
                % do nothing
            else
                error('Illegal number of input arguments');
            end

            AstS = AstSpec;  % define AstSpec

            % make sure that Units and columns have the correct num of
            % elements
            Units   = Units(1:Ncol);
            Columns = Columns(1:Ncol);

            Iwave    = find(strcmp(Columns,'Wave'));
            Iint     = find(strcmp(Columns,'Int'));
            Ierr     = find(strcmp(Columns,'Err'));
            Iback    = find(strcmp(Columns,'Back'));
            Imask    = find(strcmp(Columns,'Mask'));

            Nunits = numel(Units);
            if (isempty(Iwave))
                error('Must supply wavelength');
            else
                AstS(1).Wave      = Mat(:,Iwave);
                if (Iwave>Nunits)
                    AstS(1).WaveUnits = '';
                else
                    AstS(1).WaveUnits = Units{Iwave};
                end
            end
            if (isempty(Iint))
                error('Must supply intensity');
            else
                AstS(1).Int    = Mat(:,Iint);
                if (Iint>Nunits)
                    AstS(1).IntUnits = '';
                else
                    AstS(1).IntUnits = Units{Iint};
                end
            end
            if (~isempty(Ierr))
                AstS(1).Err = Mat(:,Ierr);
            end
            if (~isempty(Iback))
                AstS(1).Back   = Mat(:,Iback);
            end
            if (~isempty(Imask))
                AstS(1).Mask   = Mat(:,Imask);
            end

            AstS(1).WaveUnits  = Units{1};
            AstS(1).IntUnits   = Units{2};
            AstS(1).Header = Header;

        end    
        
        function AS=spec_read_mat(File,varargin)
            %--------------------------------------------------------------------------
            % spec_read_mat function                                         AstroSpec
            % Description: Read a spectrum from a mat file.
            % Input  : - File name, file name with wild cards, or a cell array of
            %            file names. Each file is a mat file in the current directory
            %            containing spectra in mtarix, cell or AstSpec formats.
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            --- Additional parameters
            %            Any additional key,val, that are recognized by one of the
            %            following programs:
            % Output : - File names to read.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Feb 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example:  AS=spec_read_mat('rg5iii.mat')
            % Reliable: 
            %--------------------------------------------------------------------------

            DefV.SpecOutType      = 'AstSpec';   % 'AstSpec','cell','mat'
            DefV.MaxCol           = 3;
            DefV.ColNames         = {'Wave','Int','Err','Back','Mask'};
            DefV.ColCorr          = [1     ,1    ,1    ,1     ,1     ];
            DefV.WaveUnits        = 'Ang';
            DefV.IntUnits         = 'erg*cm^{-2}*s^{-1}*Ang^{-1}';
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            [~,List] = Util.files.create_list(File,NaN);
            Nl       = numel(List);

            ColNames = InPar.ColNames(1:1:InPar.MaxCol);
            Ncol     = numel(ColNames);

            % initialization
            switch lower(InPar.SpecOutType)
                case 'astspec'
                    AS = AstSpec(Nl,1);   % define AstSpec class
                case 'cell'
                    AS = cell(Nl,1);
                case 'mat'
                    if (Nl>1)
                        error('mat output is possible only when reading a single file');
                    end
                otherwise
                    error('Unknown SpecOutType option');
            end


            for Il=1:1:Nl
                % for each file

                % read mat file
                Mat = Util.IO.load2(List{Il});


                % check type
                switch lower(InPar.SpecOutType)
                    case 'astspec'
                        if (isnumeric(Mat))
                            for Imc=1:1:min(InPar.MaxCol,size(Mat,2))
                                AS(Il).(ColNames{Imc}) = Mat(:,Imc);
                            end
                            Ind = Il;
                        elseif (iscell(Mat))
                            for Imc=1:1:InPar.MaxCol
                                AS(Il).(ColNames{Imc}) = Mat{Il}(:,Imc);
                            end
                            Ind = Il;
                        elseif (AstSpec.isastspec(Mat))
                            Ind = (numel(AS)+1:numel(AS)+numel(Mat))';
                            AS  = [AS, Mat];
                        else
                            error('Uknwon object in Mat file');
                        end

                        % populate fields
                        AS(Ind).WaveUnits = InPar.WaveUnits;
                        AS(Ind).IntUnits  = InPar.IntUnits;
                        AS(Ind).FileName  = List{Il};

                    case 'mat'
                        if (isnumeric(Mat))
                            AS = Mat;
                        elseif (iscell(Mat))
                            error('Mat file contain cell array and requested output is mat');
                        elseif (AstSpec.isastspec(Mat))
                            AS = astspec2mat(Mat);
                        else
                            error('Uknwon object in Mat file');
                        end
                    case 'cell'
                        if (isnumeric(Mat))
                            AS{Il} = Mat;
                        elseif (iscell(Mat))
                            AS{end+1:end+numel(Mat)} = Mat;
                        elseif (AstSpec.isastspec(Mat))
                            AS{Il} = astspec2mat(Mat);
                        else
                            error('Uknwon object in Mat file');
                        end
                    otherwise
                        error('Unknwon SpecOytType option');
                end
            end
        end
    end
    
    % Read spectra from databases
    methods (Static)

        function Spec=get_gaia_synspec(Temp,Grav,Metal,Rot,OutType)
            % Get GAIA synthetic spectra from local DB
            % Package: @AstSpec
            % Description: get a synthetic stellar spectrum from the local GAIA
            %              spectral library. Spectra are in the range 2500-10500A
            %              and 1A resolution.
            %              Assuming alpha enhanement 0, and micro-turbulence 2km/s.
            % Input  : - Effective temperature [K] vector (in range 3500-47500 K).
            %          - Gravity [log g] (in range 0 to 5) vector.
            %          - Metallicity [log solar] (in rahe -2.5 to 2.5)
            %            vector.
            %          - Rotation velocity [km/s] (in the range 0 to
            %            500km/s) vector.
            %          - Output type: 'astspec' | 'mat'. Default is 'astspec'.
            %            'mat' option can read only one file.
            % Output : - Spectrum in flux units [wavelength[Ang], Flux].
            %            Flux units [erg cm^-2 s^-1 A^-1 on star]
            %            Return NaN if spectrum doesn't exist or web site is down.
            % Reference: http://gaia.esa.int/spectralib/spectralib1A/SpectraLib1a.cfm
            % Tested : Matlab 7.3
            %     By : Eran O. Ofek                    Nov 2008
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % See also: wget_gaia_synspec.m
            % Example: Spec=AstSpec.get_gaia_synspec(5000,0,0,0);
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (nargin<5)
                OutType = 'astspec';
            end

            Nt   = numel(Temp);
            Ng   = numel(Grav);
            Nm   = numel(Metal);
            Nr   = numel(Rot);
            Nmax = max([Nt,Ng,Nm,Nr]);
            Temp = Temp.*ones(Nmax,1);
            Grav = Grav.*ones(Nmax,1);
            Metal= Metal.*ones(Nmax,1);
            Rot  = Rot.*ones(Nmax,1);
            
            switch lower(OutType)
                case 'astspec'
                    AS   = AstSpec(Nmax,1);
            end
            
            MSDir      = Util.files.which_dir(mfilename);
            %DirLocation = sprintf('%s%s..%s..%s%s%s%s%s',MSDir,filesep,filesep,filesep,'data',filesep,'GAIA_SpecTemplates',filesep);
            %DirLocation = sprintf('%s%s..%s..%s..%s%s%s%s%s',MSDir,filesep,filesep,filesep,filesep,'data',filesep,'GAIA_SpecTemplates',filesep);
            DirLocation = sprintf('%s%s..%s..%s%s%s%s%s%s%s',MSDir,filesep,filesep,filesep,'data',filesep,'spec',filesep,'GAIA_SpecTemplate',filesep);
            
            DefaultPars = 'K2SNWNVD01F.mat';
            WaveFileName = 'GAIA_Wave1A.mat';
            for Is=1:1:Nmax
                if (Metal(Is)<0)
                    MetalSign = 'M';
                else
                    MetalSign = 'P';
                end

                SpecName = sprintf('%sT%05dG%02d%s%02dV%03d%s',DirLocation,round(Temp(Is)),round(Grav(Is).*10),MetalSign,round(abs(Metal(Is)).*10),round(Rot(Is)),DefaultPars);
                W        = Util.IO.load2(sprintf('%s%s',DirLocation,WaveFileName));
                
                SpecMat  = [W, Util.IO.load2(SpecName)];

                switch lower(OutType)
                    case 'astspec'
                        Spec(Is)=AstSpec.mat2spec(SpecMat,{'Wave','Int'},{'Ang','erg*cm^-2 *s^-1*Ang^-1'});
                        Spec(Is).z = 0;
                        Spec(Is).source = 'GAIA local DB';
                        Spec(Is).ObjName = sprintf('GAIA synspec T=%f, g=%f, M=%f, R=%f',Temp(Is),Grav(Is),Metal(Is),Rot(Is));
                    case 'mat'
                        % do nothing
                    otherwise
                        error('Unknown OutType option');
                end
            end
        end
        
        function Spec=get_all_gaia_synspec
            % Get a;; GAIA synthetic spectra from local DB
            % Package: @AstSpec
            % Description: get all synthetic stellar spectrum from the local GAIA
            %              spectral library. Spectra are in the range 2500-10500A
            %              and 1A resolution.
            %              Assuming alpha enhanement 0, and micro-turbulence 2km/s.
            % Input  : *
            % Output : - AstSpec object with the spectra in flux units [wavelength[Ang], Flux].
            %            Flux units [erg cm^-2 s^-1 A^-1 on star]
            %            Return NaN if spectrum doesn't exist or web site is down.
            % Reference: http://gaia.esa.int/spectralib/spectralib1A/SpectraLib1a.cfm
            % Tested : Matlab 7.3
            %     By : Eran O. Ofek                    Nov 2008
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % See also: wget_gaia_synspec.m
            % Example: Spec=AstSpec.get_gaia_synspec(5000,0,0,0);
            % Reliable: 2
            %--------------------------------------------------------------------------

            MSDir      = Util.files.which_dir(mfilename);
            %DirLocation = sprintf('%s%s..%s..%s%s%s%s%s',MSDir,filesep,filesep,filesep,'data',filesep,'GAIA_SpecTemplates',filesep);
            %DirLocation = sprintf('%s%s..%s..%s..%s%s%s%s%s',MSDir,filesep,filesep,filesep,filesep,'data',filesep,'GAIA_SpecTemplates',filesep);
            DirLocation = sprintf('%s%s..%s..%s%s%s%s%s%s%s',MSDir,filesep,filesep,filesep,'data',filesep,'spec',filesep,'GAIA_SpecTemplate',filesep);
            
            PWD = pwd;
            cd(DirLocation);
            
            WaveFileName = 'GAIA_Wave1A.mat';
            W        = Util.IO.load2(sprintf('%s%s',DirLocation,WaveFileName));
            
            Files = dir('T*.mat');
            Nf    = numel(Files);
            Spec  = AstSpec(Nf,1);
            for If=1:1:Nf
                SpecMat  = [W, Util.IO.load2(Files(If).name)];
            
                Spec(If)=AstSpec.mat2spec(SpecMat,{'Wave','Int'},{'Ang','erg*cm^-2 *s^-1*Ang^-1'});
                Spec(If).z = 0;
                Spec(If).source = 'GAIA local DB';
                Spec(If).ObjName = sprintf('GAIA synspec %s',Files(If).name);

            end
        end
        
        
        
        function Spec=wget_gaia_synspec(Temp,Grav,Metal,Rot)
            % wget GAIA synthetic spectra from web
            % Package: @AstSpec
            % Description: wget a synthetic stellar spectrum from the GAIA spectral
            %              library, in the range 2500-10500A and 1A resolution.
            %              Assuming alpha enhanement 0, and micro-turbulence 2km/s.
            % Input  : - Effective temperature [K] vector (in range 3500-47500 K).
            %          - Gravity [log g] (in range 0 to 5) vector.
            %          - Metallicity [log solar] (in rahe -2.5 to 2.5)
            %            vector.
            %          - Rotation velocity [km/s] (in the range 0 to 500km/s)
            %            vector.
            %          - Output type: 'astspec' | 'mat'. Default is 'astspec'.
            %            'mat; option can read a single spectrum.
            % Output : - Spectrum in flux units [wavelength[Ang], Flux].
            %            Flux units [erg cm^-2 s^-1 A^-1 on star]
            %            Return NaN if spectrum doesn't exist or web site is down.
            % Reference: http://gaia.esa.int/spectralib/spectralib1A/SpectraLib1a.cfm
            % Tested : Matlab 7.3
            %     By : Eran O. Ofek                    Nov 2008
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Spec=AstSpec.wget_gaia_synspec(5000,0,0,0);
            % Reliable: 2
            %--------------------------------------------------------------------------

            error('New location of the GAIA spectra is unknown - use local library');
            
            if (nargin<5)
                OutType = 'astspec';
            end

            URL = 'http://gaia.esa.int/spectralib/spectralib1A/spectrafitsfiles/D01/';
            TempDir = sprintf('T_%05d/',round(Temp));
            
            
            Nt   = numel(Temp);
            Ng   = numel(Grav);
            Nm   = numel(Metal);
            Nr   = numel(Rot);
            Nmax = max([Nt,Ng,Nm,Nr]);
            Temp = Temp.*ones(Nmax,1);
            Grav = Grav.*ones(Nmax,1);
            Metal= Metal.*ones(Nmax,1);
            Rot  = Rot.*ones(Nmax,1);
            
            switch lower(OutType)
                case 'astspec'
                    AS   = AstSpec(Nmax,1);
            end
            
            
            DefaultPars = 'K2SNWNVD01F.fits';

            for Is=1:1:Nmax
                if (Metal(Is)<0)
                   MetalSign = 'M';
                else
                   MetalSign = 'P';
                end

                SpecName = sprintf('T%05dG%02d%s%02dV%03d%s',round(Temp(Is)),round(Grav(Is).*10),MetalSign,round(abs(Metal(Is)).*10),round(Rot(Is)),DefaultPars);

                FullURL  = sprintf('%s%s%s',URL,TempDir,SpecName);
                [Status] = system(sprintf('wget %s',FullURL));

                switch Status
                 case 0
                    % ok
                    SpecT = fitsread(SpecName,'BinTable');
                    SpecMat  = [SpecT{1}, SpecT{2}];
                    delete(SpecName);
                 otherwise
                    % not ok - spec not exist or web is down
                    SpecMat = NaN;
                end


                switch lower(OutType)
                    case 'astspec'
                        Spec(Is)=AstSpec.mat2spec(SpecMat,{'Wave','Int'},{'Ang','erg*cm^-2 *s^-1*Ang^-1'});
                        Spec(Is).z = 0;
                        Spec(Is).source = 'GAIA local DB';
                        Spec(Is).ObjName = sprintf('GAIA synspec T=%f, g=%f, M=%f, R=%f',Temp(Is),Grav(Is),Metal(Is),Rot(Is));
                    case 'mat'
                        % do nothing
                    otherwise
                        error('Unknown OutType option');
                end
            end
        end
           
        function Spec=zodiac_spectrum(Wave,OutType)
            % Get the Zodiac light spectrum
            % Package: @AstSpec
            % Description: Return the zodiac spectrum as adopted from the HST STIS
            %              handbook. The high zodiacal ligh is defined where V=22.1
            %              mag/arcsec^-2.
            % Input  : - If empty than return the zodiacal spectrum.
            %            Otherwise, if a list of wavelength [A] are provided,
            %            interpolate the zodiacal flux at the requested wavelength [A].
            %          - Output type: 'mat'|'astspec'. Default is 'astspec'.
            % Output : - Zodiacal ligh spectrum
            %            [wavelength(Ang), Flux(erg/cm^2/s/A/arcsec^2)]
            % Reference: http://www.stsci.edu/hst/stis/documents/handbooks/currentIHB/c06_exptime7.html#695784
            %            but there is a discrepency with:
            %            http://www.stsci.edu/hst/wfc3/design/documents/handbooks/currentIHB/c09_exposuretime08.html#389841
            %            According to the HST help desk the STIS table should be used.
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Nov 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Spec=AstSpec.zodiac_spectrum;
            %          % to verify normalization: synphot(Spec,'Johnson','V','Vega')
            %          Flux=AstSpec.zodiac_spectrum([5000;5500])
            % Reliable: 1
            %--------------------------------------------------------------------------
            InterpMethod = 'linear';

            
            if (nargin==0)
                Wave = [];
                OutType = 'astspec';
            elseif (nargin==1)
                OutType = 'astspec';
            else
                % do nothing
            end

            % Wave, High Earthshine, High zodi, Total back
            Spec=[
            1000 2.41e-23 9.69e-29 2.41e-23
            1100 4.38e-22 1.04e-26 4.38e-22
            1200 4.01e-23 1.08e-25 4.03e-23
            1300 7.41e-25 6.59e-25 1.40e-24
            1400 4.29e-25 2.55e-24 2.98e-24
            1500 4.16e-25 9.73e-24 1.01e-23
            1600 2.55e-25 2.35e-22 2.35e-22
            1700 7.89e-25 7.21e-21 7.21e-21
            1800 9.33e-23 1.53e-20 1.54e-20
            1900 4.39e-22 2.25e-20 2.29e-20
            2000 1.01e-21 3.58e-20 3.68e-20
            2100 1.60e-21 1.23e-19 1.25e-19
            2200 7.49e-22 2.21e-19 2.22e-19
            2300 3.32e-22 1.81e-19 1.81e-19
            2400 2.50e-22 1.83e-19 1.83e-19
            2500 2.39e-22 2.53e-19 2.53e-19
            2600 5.62e-22 3.06e-19 3.06e-19
            2700 6.77e-21 1.01e-18 1.02e-18
            2800 2.03e-21 2.88e-19 2.90e-19
            2900 4.32e-20 2.08e-18 2.12e-18
            3000 9.34e-20 1.25e-18 1.35e-18
            3100 2.07e-19 1.50e-18 1.70e-18
            3200 3.60e-19 2.30e-18 2.66e-18
            3300 4.27e-19 2.95e-18 3.38e-18
            3400 6.40e-19 2.86e-18 3.50e-18
            3500 8.20e-19 2.79e-18 3.61e-18
            3600 1.06e-18 2.74e-18 3.80e-18
            3700 1.22e-18 3.32e-18 4.54e-18
            3800 1.23e-18 3.12e-18 4.35e-18
            3900 1.52e-18 3.34e-18 4.86e-18
            4000 2.38e-18 4.64e-18 7.01e-18
            4250 2.38e-18 4.65e-18 7.03e-18
            4500 2.86e-18 5.58e-18 8.44e-18
            4750 2.79e-18 5.46e-18 8.25e-18
            5000 2.63e-18 5.15e-18 7.77e-18
            5250 2.67e-18 5.37e-18 8.04e-18
            5500 2.58e-18 5.34e-18 7.92e-18
            5750 2.54e-18 5.40e-18 7.94e-18
            6000 2.42e-18 5.25e-18 7.67e-18
            6250 2.26e-18 5.02e-18 7.28e-18
            6500 2.17e-18 4.92e-18 7.09e-18
            6750 2.07e-18 4.79e-18 6.87e-18
            7000 1.93e-18 4.55e-18 6.48e-18
            7250 1.85e-18 4.43e-18 6.29e-18
            7500 1.74e-18 4.23e-18 5.97e-18
            7750 1.63e-18 4.04e-18 5.67e-18
            8000 1.56e-18 3.92e-18 5.49e-18
            8250 1.48e-18 3.76e-18 5.23e-18
            8500 1.35e-18 3.50e-18 4.85e-18
            8750 1.31e-18 3.43e-18 4.74e-18
            9000 1.22e-18 3.23e-18 4.44e-18
            9250 1.15e-18 3.07e-18 4.21e-18
            9500 1.10e-18 2.98e-18 4.08e-18
            9750 1.04e-18 2.86e-18 3.91e-18
            10000 1.00e-18 2.78e-18 3.78e-18
            10250 9.54e-19 2.67e-18 3.63e-18
            10500 9.04e-19 2.56e-18 3.46e-18
            10750 8.41e-19 2.41e-18 3.25e-18
            11000 8.03e-19 2.31e-18 3.11e-18];

            Spec = Spec(:,[1 3]);

            if (~isempty(Wave))
                Spec = interp1(Spec(:,1),Spec(:,2),Wave,InterpMethod);
            end
            
            switch lower(OutType)
                case 'astspec'
                    AS = AstSpec;
                    Is = 1;
                    AS(Is).Wave = Spec(:,1);
                    AS(Is).Int  = Spec(:,2);
                    AS(Is).WaveUnits = 'Ang';
                    AS(Is).IntUnits  = 'erg*cm^-2*s^-1*Ang^-1';
                    AS(Is).ObjName   = 'Zodiac spectrum';
                    AS(Is).source    = 'HST STIS handbook';
                    AS(Is).z         = 0;
                    Spec = AS;
                    
                otherwise
                    % do nothing
            end
        end
        
        function [ZodiMag,ZodiVmag,Flux,Counts]=zodiac_bck(Long,Lat,Date,Filter_family,Filter_name,FilterSys)
            % Get the Zodiac light surface brightness as a function of coordinates
            % Package: @AstSpec
            % Description: Calculate the zodiac magnitude and flux in a given filter
            %              and sky position. The zodiac spectrum and position
            %              dependent flux are adopted from the HST WFC3 handbook.
            % Input  : - Either Ecliptic longitude or Helio-ecliptic longitude
            %            (i.e., L - L_sun) in radians
            %            (see convertdm.m for additional options).
            %            By default this should be the Helio-ecliptic longitude.
            %          - Ecliptic latitude [radians].
            %            (see convertdm.m for additional options).
            %          - Date [D M Y] or JD at which to calculate the Helio-ecliptic
            %            longitude. Defaut is empty. If empty than treat the first
            %            input argument as Helio-ecliptic longitude.
            %          - Filter family (e.g., 'GALEX'). See get_filter.m for
            %            options. Default is 'SDSS'.
            %          - Filter band name (e.g., 'NUV'). See get_filter.m for
            %            options. Default is 'r'.
            %          - Mag system {'AB'|'Vega'}. Default is 'AB'.
            % Output : - Zodiacal light magnitude in filter and position.
            %          - Zodiacal light V-band magnitude in position.
            %          - Zodiacal light flux in filter [erg cm^-2 s^-1 arcsec^2]
            %          - Zodiacal light photon rate [photons cm^-2 s^-1 arcsec^2]
            % Tested : Matlab R2014a
            %     By : Ilan Sagiv                      Sep 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Reference: https://hst-docs.stsci.edu/display/WFC3IHB/9.7+Sky+Background#id-9.7SkyBackground-9.7.19.7.1ZodiacalLight,EarthShine,andLOW-SKY
            % Example: [F,V,Flux,Counts]=AstSpec.zodiac_bck(45./RAD,80./RAD,[],'LIM','NUV')
            % Reliable: NEED TO VERIFY zodi spec
            %--------------------------------------------------------------------------
            ReplaceNaNVal  = 21.5;

            InterpMethod   = 'linear';
            RAD            = 180./pi;
            h              = constant.h; %('h','cgs');
            c              = constant.c; %('c','cgs');
            hc             = h*c*1e8; %[erg]*[Ang]

            Def.Date          = [];
            Def.Filter_family = 'SDSS';
            Def.Filter_name   = 'r';
            Def.FilterSys     = 'AB';
            if (nargin==2)
                Date          = Def.Date;
                Filter_family = Def.Filter_family;
                Filter_name   = Def.Filter_name;
                FilterSys     = Def.FilterSys;
            elseif (nargin==3)
                Filter_family = Def.Filter_family;
                Filter_name   = Def.Filter_name;  
                FilterSys     = Def.FilterSys;
            elseif (nargin==4)
                Filter_name   = Def.Filter_name;   
                FilterSys     = Def.FilterSys;
            elseif (nargin==5)
                FilterSys     = Def.FilterSys;    
            else
                % do nothing
            end


            if (~isempty(Date))
                % assume ecliptic longitude as input
                % Transform Ecliptic coordinates to Helio-ecliptic coordinates
                HelioEcLong = celestial.coo.ecliptic2helioecliptic(Long,Date);
            else
                HelioEcLong = Long;
            end

            % convert to 0-pi range
            HelioEcLong(HelioEcLong>pi) = 2.*pi - HelioEcLong(HelioEcLong>pi);


            % get sodi spectrum
            Spec=AstSpec.zodiac_spectrum;


            %wavelength=data(:,1);
            %Zodiac    =data(:,3);
            %Spectrum=[wavelength,Zodiac];

            %if(0)% plot Spectrum
            %	figure()
            %	semilogy(wavelength,Zodiac,'r');
            %	xlabel('Wavelength [A]')
            %	ylabel('Flux [erg / sec\cdotcm^2 \cdot A]')
            %    axis([1000 12e3 1e-25 1e-16])
            %end


            % approximate zodiacal sky background as a function of
            % Heliocentric ecliptic longitude and ecliptic latitude
            % (in Vmag / arcsec^2)
            % from HST WFC3 Handbook; (table 9.4)
            % V-band magnitude of zodi
            % adopted from Table 9.4 in:
            % http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c09_exposuretime08.html#389841
            EcLon=(0:15:180)./RAD;
            EcLat=(0:15:90)./RAD;
            %lat-->
            %v lon
            ZodiVmagTable=[...
            nan ,nan ,nan ,nan ,22.6,23.0,23.3;
            nan ,nan ,nan ,nan ,22.6,23.1,23.3;
            nan ,nan ,nan ,22.3,22.7,23.1,23.3;
            nan ,nan ,22.1,22.5,22.9,23.1,23.3;
            21.3,21.9,22.4,22.7,23.0,23.2,23.3;
            21.7,22.2,22.6,23.9,23.1,23.2,23.3;
            22.0,22.3,22.7,23.0,23.2,23.3,23.3;
            22.2,22.5,22.9,23.1,23.3,23.3,23.3;
            22.4,22.6,22.9,23.2,23.3,23.3,23.3;
            22.4,22.6,22.9,23.2,23.3,23.3,23.3;
            22.4,22.6,22.9,23.1,23.3,23.3,23.3;
            22.3,22.5,22.8,23.0,23.2,23.4,23.3;
            22.1,22.4,22.7,23.0,23.2,23.4,23.3];

            % Interpolate zodi Vmag to requested positions
            ZodiVmag = interp2(EcLat,EcLon,ZodiVmagTable,abs(Lat),HelioEcLong,InterpMethod);
            ZodiVmag(isnan(ZodiVmag)) = ReplaceNaNVal;

            Vmag0 = 22.1;   % synphot(Spec,'Johnson','V','Vega')
            DeltaMag = ZodiVmag-Vmag0;

            ZodiMag = DeltaMag + synphot(Spec,Filter_family,Filter_name,FilterSys);

            %Filter = get_astfilter(Filter_family,Filter_name);
            Filter = AstFilter.get(Filter_family,Filter_name);
            FilterTr = [Filter.nT(:,1), Filter.nT(:,2)./max(Filter.nT(:,2))];
            %[Flux,Counts] = spec_photon_counts([Spec(:,1),Spec(:,2).*10.^(-0.4.*DeltaMag)],FilterTr,[],1,1./3.08e18);
            [Flux,Counts] = spec_photon_counts([Spec.Wave,Spec.Int.*10.^(-0.4.*DeltaMag)],FilterTr,[],1,1./3.08e18);

        end
        
        function AstS=get_pickles(SpC,SpL)
            % load stellar spectrum from the Pickles stellar spectra library
            % Package: @AstSpec
            % Description: Load the Pickles stellar spectra library into a AstSpec
            %              class object. By default will load all the normal stars.
            % Input  : - Optional spectral type to load (e.g., 'G'). If empty, load
            %            all. Default is empty.
            %          - Optional luminosity class to load (e.g., 'V'). If empty, load
            %            all. Default is empty.
            % Output : - AstSpec class object will all the spectra.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Feb 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: AstS=AstSpec.get_pickles;
            %          AstS=AstSpec.get_pickles('g','v');
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (nargin==0)
                SpC = [];
                SpL = [];
            elseif (nargin==1)
                SpL = [];
            else
                % do nothing
            end


            DirPic = '~/matlab/data/spec/PicklesStellarSpec';
            PWD    = pwd;
            cd(DirPic);

            [~,List] = Util.files.create_list('uk*.mat',NaN);
            RE  = regexp(List,'uk(?<SpClass>[obafgkm]+)(?<SpNum>\d+)(?<SpLum>\w+)','names');
            Nre = numel(RE);

            Ind = 0;
            for Ire=1:1:Nre
                if (~isempty(RE{Ire}))
                    Ind = Ind + 1;
                    SRE(Ind) = RE{Ire};
                    ListN{Ind} = List{Ire};
                end
            end

            Nre = numel(SRE);
            if (isempty(SpC))
                FlagC = true(1,Nre);
            else
                FlagC = strcmpi({SRE.SpClass},SpC);
            end
            if (isempty(SpL))
                FlagL = true(1,Nre);
            else
                FlagL = strcmpi({SRE.SpLum},SpL);
            end

            Flag = FlagC & FlagL;


            %AstS = spec_read_mat('uk*.mat');
            AstS = AstSpec.spec_read_mat(ListN(Flag));
            
            SRE  = SRE(Flag);

            Ns   = numel(AstS);
            for Is=1:1:Ns
                AstS(Is).z = 0;
                AstS(Is).source = 'Pickles 1998';

                if (numel(SRE(Is).SpNum)>1)
                    SpNum = str2double(SRE(Is).SpNum(1))+0.5;
                else
                    SpNum = str2double(SRE(Is).SpNum);
                end
                
                AstS(Is).ObjName = sprintf('%s %3.1f %s',upper(SRE(Is).SpClass),SpNum,upper(SRE(Is).SpLum));

            end

            cd(PWD);
        end
        
        function AstS=get_galspec(Name,OutType)
            % get galaxy/qso template spectrum
            % Packge: @AstSpec
            % Description: Get Galaxy or QSO spectral template from local
            %              database.
            % Input  : - Spectral template name. If empty, then will show
            %            the list of available spectral templates.
            %          - Output type: 'astspec'|'mat'. Default is
            %            'astspec'.
            % Output : - An AStSpect object or a matrix with the requested
            %            spectrum. If the first argument is empty then this
            %            is a cell array of available file names.
            % Example: AstSpec.get_galspec
            %          A=AstSpec.get_galspec('QSO_NIR');
            % Reliable: 2
            
            if (nargin==0)
                Name = [];
            elseif (nargin==1)
                OutType = 'astspec';
            else
                % do nothing
            end
            
            Dir = Util.files.which_dir(mfilename);
            SpecDir = 'SpecGalQSO';
            DirPath = sprintf('~/matlab/data/spec/%s%s',SpecDir,filesep);
            
            if (isempty(Name))
                % show all spectral template available
                AllFiles = dir(sprintf('%s*.txt',DirPath));
                Nf=numel(AllFiles);
                % show content of all available files
                sprintf('\n Available spectral templates:\n');
                for If=1:1:Nf
                    FID=fopen(sprintf('%s%s',DirPath,AllFiles(If).name));
                    Line = fgetl(FID);
                    fclose(FID);
                    fprintf('%20s  : %s\n',AllFiles(If).name,Line);
                end
                AstS = {AllFiles.name};
            else
                % get a single file
                if (isempty(strfind(Name,'.txt')))
                    % add .txt to Name
                    Name = sprintf('%s.txt',Name);
                end
                Spec = Util.IO.load2(sprintf('%s%s',DirPath,Name));
                switch lower(OutType)
                    case 'mat'
                        AstS = Spec;
                    case 'astspec'
                        AstS = AstSpec;
                        AstS.Wave = Spec(:,1);
                        AstS.Int  = Spec(:,2);
                        AstS.WaveUnits = 'Ang';
                        AstS.ObjName   = Name;
                        AstS.z         = 0;
                        
                    otherwise
                        error('Unknown OutType option');
                end
            end
        end
        
        function AstS=get_atmospheric_extinction(Name,OutType)
            % get atmopsheric extinction curve for various observatories
            % Package: @AstSpec
            % Description: Get atmopsheric extinction curve for various observatories
            %              local database.
            % Input  : - File name. If empty, then will show
            %            the list of available atmospheric extinction files.
            %          - Output type: 'astspec'|'mat'. Default is
            %            'astspec'.
            % Output : - An AStSpect object or a matrix with the requested
            %            curve. If the first argument is empty then this
            %            is a cell array of available file names.
            % Example: AstSpec.get_atmospheric_extinction
            %          A=AstSpec.get_atmospheric_extinction('KPNO');
            % Reliable: 2
            
            if (nargin==0)
                Name = [];
            elseif (nargin==1)
                OutType = 'astspec';
            else
                % do nothing
            end
            
            
            
            if (isempty(Name))
                Dir = Util.files.list_fun_in_package('cats.spec.AtmoExtinction');
                AllFiles = Dir.m;
                % show all spectral template available
                
                Nf=numel(AllFiles);
                % show content of all available files
                sprintf('\n Available spectral templates:\n');
                AstS = cell(Nf,1);
                for If=1:1:Nf
                    [~,FileName] = fileparts(AllFiles{If});
                    
                    fprintf('%20s  \n',FileName);
                    AstS{If} = FileName;
                end
            else
                % get a single file
                Spec = cats.spec.AtmoExtinction.(Name);
                
                
                switch lower(OutType)
                    case 'mat'
                        AstS = Spec;
                    case 'astspec'
                        AstS = AstSpec;
                        AstS.Wave = Spec(:,1);
                        AstS.Int  = Spec(:,2);
                        AstS.WaveUnits = 'Ang';
                        AstS.IntUnits  = 'mag/airmass';
                        AstS.ObjName   = Name;
                        AstS.z         = 0;
                    otherwise
                        error('Unknown OutType option');
                end
            end
            
        end
                
        function [AstS]=blackbody(VecT,VecW,UnitsOut,UnitsWave,OutType)
            % Calculate a blackbody spectrum into an AstSpec object.
            % Description: Calculate black body (planck) spectrum.
            % Input  : - Array of temperatures [K]. Calculate a planck spectrum for
            %            each temperature.
            %          - Vector of wavelength. Default is (1000:100:25000)';
            %          - Output units. See convert.flux for options.
            %            Default is 'cgs/A'.
            %          - Input wavelength units. See convert.units for options.
            %            Default is 'ang'.
            %          - Output type: 'astspec'|'mat'. Default is 'astspec'.
            %            'mat' works on a single AstSpec element.
            % Output : - AstSpec array with black body spectrum.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: AstS=AstSpec.blackbody([5000;10000]);
            % Reliable: 2
            %--------------------------------------------------------------------------

            Def.VecW      = (1000:100:25000).';
            Def.UnitsOut  = 'cgs/A';
            Def.UnitsWave = 'Ang';
            Def.OutType   = 'astspec';
            if (nargin==1)
                VecW      = Def.VecW;
                UnitsOut  = Def.UnitsOut;
                UnitsWave = Def.UnitsWave;
                OutType   = Def.OutType;
            elseif (nargin==2)
                UnitsOut  = Def.UnitsOut;
                UnitsWave = Def.UnitsWave;
                OutType   = Def.OutType;
            elseif (nargin==3)
                UnitsWave = Def.UnitsWave;
                OutType   = Def.OutType;
            elseif (nargin==4)
                OutType   = Def.OutType;
            elseif (nargin==5)
                % do nothing
            else
                error('Illegal number of input arguments');
            end

            h      = constant.h; %6.6261e-27;      % = get_constant('h','cgs');          % Planck constant [cgs] 
            c      = constant.c; %29979245800;     % = get_constant('c','cgs');          % speed of light [cm]
            k      = constant.kB; %1.380648813e-16; % = get_constant('kB','cgs');         % Boltzmann constant [cgs]

            VecW = VecW(:);
            Lam = VecW.*convert.units(UnitsWave,'cm');   % convert wavelength to cm
            %Lam    = W.*1e-8;            % convert Ang to cm
            Nu     = c./Lam;             % convert wavelength to frequency [Hz]

            Nt = numel(VecT);
            AstS = AstSpec(size(VecT));
            for It=1:1:Nt
                % for each temperature
                AstS(It).Wave = VecW;

                %----------------------
                %--- planck formula ---
                %----------------------
                % emittence [erg/sec/cm^2/Ang(lambda)]
                Flambda = 1e-8.*2.*pi.*h.*c.^2.*Lam.^(-5) ./ (exp(h.*c./(Lam.*k.*VecT(It))) - 1);

                AstS(It).Int = convert.flux(Flambda,'cgs/A',UnitsOut,Lam,'cm');
                AstS(It).WaveUnits = UnitsWave;
                AstS(It).IntUnits  = UnitsOut;
                AstS(It).source    = 'AstSpec.blackbody';
                AstS(It).ObjName   = sprintf('Planck spectrum T=%f',VecT(It));
                AstS(It).z         = 0;
            end
            
            switch lower(OutType)
                case 'mat'
                    AstS = astspec2mat(AstS,'mat');
            end
        end  % AstSpec.blackbody function
        
        
        
    end
       
    %--------------------------
    %--- Non-static classes ---
    %--------------------------
    
    % conversion to other classes
    methods
        % convert to matrix
        function Spec=astspec2mat(AS,Type)
            % Convert AstSpec object to matrix or cell of matrices
            % Description: Convert an AstSpec class object to a matrix or a cell
            %              array of matrices.
            % Input  : - AstSpec class object.
            %          - Output type {'cell','mat'}. Default is 'mat'.
            % Output : - Matrix or cell array of spectra.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Feb 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: A=AstSpec.get_galspec('QSO_NIR'); Spec=astspec2mat(A);
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (nargin==1)
                Type = 'mat';
            end

            if (~AstSpec.isastspec(AS))
                error('First input must be of AstSpec class');
            end

            Ns    = numel(AS);

            switch lower(Type)
                case 'mat'
                    if (Ns>1)
                        warning('More than one output spectrum - Changing Type to cell');
                        Type = 'cell';
                        Spec = cell(Ns,1);
                    end
                case 'cell'
                    Spec = cell(Ns,1);
                otherwise
                    % do nothing
            end
            for Is=1:1:Ns
                % for each spectrum

                Mat = [AS(Is).Wave, AS(Is).Int];
                if (~isempty(AS(Is).Err))
                    Mat = [Mat, AS(Is).Err];
                end
                if (~isempty(AS(Is).Back))
                    Mat = [Mat, AS(Is).Back];
                end
                if (~isempty(AS(Is).Mask))
                    Mat = [Mat, AS(Is).Mask];
                end

                switch lower(Type)
                    case 'mat'
                        Spec = Mat;
                    case 'cell'
                        Spec{Is} = Mat;
                    otherwise
                        error('Unknown Type option');
                end
            end
        end
        
        function St=astspec2struct(AstS)
            % Convert an AstSpec object to structure
            % Description: Convert an AstSpec class to a structure array
            % Input  : - AstSpec class object.
            % Output : - A structure array.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: A=AstSpec.get_galspec('QSO_NIR'); St=astspec2struct(A)
            % Reliable: 2
            %--------------------------------------------------------------------------

            Ns = numel(AstS);
            Size = size(AstS);
            St = struct('Wave',cell(Size),...
                        'Int',cell(Size),...
                        'Err',cell(Size),...
                        'Back',cell(Size),...
                        'Mask',cell(Size),...
                        'WaveUnits',cell(Size),...
                        'IntUnits',cell(Size),...
                        'AddCol',cell(Size),...
                        'ObjName',cell(Size),...
                        'Header',cell(Size),...
                        'comments',cell(Size),...
                        'source',cell(Size),...
                        'FileName',cell(Size),...
                        'z',cell(Size),...
                        'UserData',cell(Size));
            Fields = fieldnames(AstS(1));
            Nf     = numel(Fields);
            for Is=1:1:Ns
                for If=1:1:Nf
                    St(Is).(Fields{If})      = AstS(Is).(Fields{If});
                end
            end
        end
        
        
    end
    
    % get units   
    methods
        function WaveUnits=get_wave_units(AstS)
            % Description: Get wavelength units of AstSpec class
            % Input  : - AstSpec class.
            % Output : - A cell array of string containing the AstSpec
            %            wavelength units.
            
            
            Ns = numel(AstS);
            WaveUnits = cell(Ns,1);
            for Is=1:1:Ns
                if (isempty(AstS(Is).WaveUnits))
                    WaveUnits{Is} = '';
                else
                    WaveUnits{Is} = AstS(Is).WaveUnits;
                end
            end
            
        end

        function IntUnits=get_int_units(AstS)
            % Description: Get intensity units of AstSpec class
            % Input  : - AstSpec class.
            % Output : - A cell array of string containing the AstSpec
            %            intensity units.
            
            Ns = numel(AstS);
            IntUnits = cell(Ns,1);
            for Is=1:1:Ns
                if (isempty(AstS(Is).IntUnits))
                    IntUnits{Is} = '';
                else
                    IntUnits{Is} = AstS(Is).IntUnits;
                end
            end
            
        end
    end
        
    % convert units
    methods
        
        function AstS=convert_wave(AstS,OutUnits)
            % Description: Convert wavelength units of AstSpec class object.
            % Input  : - AstSpec class object.
            %          - Output units. See convert.units.m for options.
            % Outout : - AstSpec class object in which the .Wave field
            %            units are converted to the new system.
            % Example: AstS=AstSpec.get_pickles;
            %          AstS=convert_wave(AstS,'micron')
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                if (~isempty(AstS(Is).Wave))
                    ConvFactor   = convert.units(lower(AstS(Is).WaveUnits),lower(OutUnits));
                    AstS(Is).Wave = AstS(Is).Wave .* ConvFactor;
                    AstS(Is).WaveUnits  = OutUnits;
                    
                end
            end
        end
        
        function [F,WaveUnits]=calc_freq(AstS,WaveUnits)
            % Description: Calculate the frequency [Hz] given
            %              AstSpec class object.
            % Input  : - AstSpec class object.
            %          - Optional cell array of Wavelength units in
            %            AstSpec.
            % Output : - Structure array containing a field (.Freq) with
            %            a vector of frquencies [Hz] for each spectra.
            %          - Cell array of wavelength units in AstSpec.
            
            if (nargin==1)
                WaveUnits = get_wave_units(AstS);
            end
            c = constant.c;
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                % for each spectrum
                if (isempty(WaveUnits{Is}))
                    % no wavelength - assume Ang
                    warning('No wavelength units for spectrum %d - Assume ang',Is);
                    Units = 'ang';
                else
                    Units = WaveUnits{Is};
                end
                
                ConvFactor = convert.units(lower(Units),'cm');
                F(Is).Freq = c./(ConvFactor.*AstS(Is).Wave);
            end
                
        end 
        
        function [F,WaveUnits]=calc_energy(AstS,EnergyUnits,WaveUnits)
            % Description: Calculate the energy given
            %              AstSpec class object.
            % Input  : - AstSpec class object.
            %          - Output energy units (e.g.,'erg','J','eV','GeV','T',...)
            %            See convert_energy for options.
            %            Default is 'eV'.
            %          - Optional cell array of Wavelength units in
            %            AstSpec.
            % Output : - Structure array containing a field (.E) with
            %            a vector of energies for each spectra.
            %          - Cell array of wavelength units in AstSpec.
            
            if (nargin<2)
                EnergyUnits = 'eV';
                WaveUnits   = get_wave_units(AstS);
            elseif (nargin<3)
                WaveUnits   = get_wave_units(AstS);
            else
                % do nothing
            end
            c = constant.c;
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                % for each spectrum
                if (isempty(WaveUnits{Is}))
                    % no wavelength - assume Ang
                    warning('No wavelength units for spectrum %d - Assume ang',Is);
                    Units = 'ang';
                else
                    Units = WaveUnits{Is};
                end
                
                F(Is).E = convert.energy(Units,EnergyUnits,AstS(Is).Wave);
                
            end
                
        end 
    end
        
        
    % plots
    methods 
        function [Hx,Hy]=plot_labels(AstS,Type)
            % Description: Add xlabel and ylabel to spectral plot
            % Input  : - AstSpec class object
            %          - 'wave' | 'freq'. Default is 'wave'
            % Output : - xlabel handle
            %          - ylabel handle
            
            if (nargin<2)
                Type = 'wave';
            end
            
            WaveUnits = get_wave_units(AstS);
            if (numel(unique(WaveUnits))==1)
                WaveUnitsStr = WaveUnits{1};
            else
                warning('AstS class contains spectra with different wavelength units');
                WaveUnitsStr = '';
            end
            
            IntUnits = get_int_units(AstS);
            if (numel(unique(IntUnits))==1)
                IntUnitsStr = IntUnits{1};
            else
                warning('AstS class contains spectra with different intensity units');
                IntUnitsStr = '';
            end
            
            switch lower(Type)
                case 'wave'
                    Hx = xlabel(sprintf('Wavelength [$%s$]',WaveUnitsStr));
                case 'freq'
                    Hx = xlabel(sprintf('Frequency [$1/%s$]',WaveUnitsStr));
                otherwise
                    error('Uknown Type option');
            end
            Hx.Interpreter = 'latex';
            Hx.FontSize    = 16;
            
            Hy = ylabel(sprintf('Intensity [$%s$]',IntUnitsStr));
            Hy.Interpreter = 'latex';
            Hy.FontSize    = 16;
        end
        
        function [H,Hx,Hy]=plot(AstS,varargin)
            % Description: Given an AstSpec class, plot all the spectra
            %              as a function of wavelength.
            % Input  : - AstFilter class.
            %          * Additional arguments to pass to plot.m
            % Output : - Plot handle.
            %          - xlabel handle.
            %          - ylabel handle.
            % Example: A=AstSpec.get_galspec('QSO_NIR'); plot(A)
            
            Ns = numel(AstS);
            H  = zeros(Ns,1);
            for Is=1:1:Ns
                H = plot(AstS(Is).Wave,AstS(Is).Int,varargin{:});
                hold on;
            end
            hold off;
            [Hx,Hy]=plot_labels(AstS,'wave');
        end
            
        function [H,Hx,Hy]=plotf(AstS,varargin)
            % Description: Given an AstSpec class, plot all the spectra
            %              as a function of frequency.
            % Input  : - AstFilter class.
            %          * Additional arguments to pass to plot.m
            % Output : - Plot handle.
            %          - xlabel handle.
            %          - ylabel handle.
            
            Ns = numel(AstS);
            WaveUnits = get_wave_units(AstS);
            [F]=calc_freq(AstS,WaveUnits);
            for Is=1:1:Ns
                H = plot(F(Is).Freq,AstS(Is).Int,varargin{:});
                hold on;
            end
            hold off;
            
            [Hx,Hy]=plot_labels(AstS,'freq');
            
        end
            
        function [H,Hx,Hy]=stairs(AstS,varargin)
            % Description: Given an AstSpec class, plot all the spectra
            %              as a function of wavelength using stairs.
            % Input  : - AstFilter class.
            %          * Additional arguments to pass to plot.m
            % Output : - Plot handle.
            %          - xlabel handle.
            %          - ylabel handle.
            
            Ns = numel(AstS);
            H  = zeros(Ns,1);
            for Is=1:1:Ns
                H = stairs(AstS(Is).Wave,AstS(Is).Int,varargin{:});
                hold on;
            end
            hold off;
            
            [Hx,Hy]=plot_labels(AstS,'wave');
            
        end    
        
        function [H,Hx,Hy]=stairsf(AstS,varargin)
            % Description: Given an AstSpec class, plot all the spectra
            %              as a function of frequency using stairs.
            % Input  : - AstFilter class.
            %          * Additional arguments to pass to plot.m
            % Output : - Plot handle.
            %          - xlabel handle.
            %          - ylabel handle.
            
            Ns = numel(AstS);
            WaveUnits = get_wave_units(AstS);
            [F]=calc_freq(AstS,WaveUnits);
            for Is=1:1:Ns
                H = stairs(F(Is).Freq,AstS(Is).Int,varargin{:});
                hold on;
            end
            hold off;
            
            [Hx,Hy]=plot_labels(AstS,'freq');
            
        end
    end
        
    % Statistics   
    methods
        function [Int]=integral(AstS)
            % Description: Calculate the integral of spectra in AstSpec
            %              class.
            % Input  : - AstSpec class.
            % Output : - A vector of integrals of the .Spec field in each
            %            of the AstSpec class elements.
            % Example: [Int]=AstS.integral;
            
            Nf = numel(AstS);
            SizeS = size(AstS);
            
            Int  = zeros(SizeS);
            for If=1:1:Nf
                if (isempty(AstS(If).Int) || isempty(AstS(If).Wave))
                    Int(If)  = NaN;
                else
                    % ignore NaN's
                    Fnn      = ~isnan(AstS(If).Int);
                    Int(If)  = trapz(AstS(If).Wave(Fnn), AstS(If).Int(Fnn));
                    
                end
            end
            
        end
        
        function [MaxVal,MaxWave,MaxInd]=max(AstS)
            % Description: Get the maximum in each spectrum in a AstSpec
            %              class. Return NaN when spectra is empty.
            % Input  : - AstSpec class.
            % Output : - A vector of maximum in each AstSpec spectrum.
            %          - A vector of the wavelength at the maximum value.
            %          - A vector of indices of the maximum value.
            % Example: [MaxVal,MaxWave,MaxInd]=max(AstS);
            
            Ns = numel(AstS);
            SizeS = size(AstS);
            
            MaxVal   = zeros(SizeS).*NaN;
            MaxWave  = zeros(SizeS).*NaN;
            MaxInd   = zeros(SizeS).*NaN;
            for Is=1:1:Ns
                if (~isempty(AstS(Is).Int) && ~isempty(AstS(Is).Wave))
                    
                    [MaxVal(Is),MaxInd(Is)] = max(AstS(Is).Int);
                    if (nargout>1)
                        MaxWave(Is) = AstS(Is).Wave(MaxInd(Is));

                    end
                end
            end
        end
        
        function [MinVal,MinWave,MinInd]=min(AstS)
            % Description: Get the minimum in each spectrum in a AstSpec
            %              class. Return NaN when spectra is empty.
            % Input  : - AstSpec class.
            % Output : - A vector of minimum in each AstSpec spectrum.
            %          - A vector of the wavelength at the minimum value.
            %          - A vector of indices of the minimum value.
            % Example: [MinVal,MinWave,MinInd]=min(AstS);
            
            Ns = numel(AstS);
            SizeS = size(AstS);
            
            MinVal   = zeros(SizeS).*NaN;
            MinWave  = zeros(SizeS).*NaN;
            MinInd   = zeros(SizeS).*NaN;
            for Is=1:1:Ns
                if (~isempty(AstS(Is).Int) && ~isempty(AstS(Is).Wave))
                    
                    [MinVal(Is),MinInd(Is)] = min(AstS(Is).Int);
                    if (nargout>1)
                        MinWave(Is) = AstS(Is).Wave(MinInd(Is));

                    end
                end
            end
        end
        
        function [WaveRange,IntRange,WaveRatio,IntRatio]=range(AstS)
            % Description: Get the rane and ratios of the maximum to 
            %              minimum wavelength and intensity in each
            %              spectrum in a AstSpec class.
            %              Return NaN when spectra is empty.
            % Input  : - AstSpec class.
            % Output : - A vector of wavelength range for each AstSpec
            %            spectrum.
            %          - A vector of intensity range for each AstSpec
            %            spectrum.
            %          - A vector of max to min wavelength ratio for
            %            each AstSpec spectrum.
            %          - A vector of max to min intensity ratio for
            %            each AstSpec spectrum.
            % Example: [WaveRange]=range(AstS);
            
            Ns = numel(AstS);
            SizeS = size(AstS);
            
            WaveRange   = zeros(SizeS).*NaN;
            IntRange    = zeros(SizeS).*NaN;
            WaveRatio   = zeros(SizeS).*NaN;
            IntRatio    = zeros(SizeS).*NaN;
            for Is=1:1:Ns
                if (~isempty(AstS(Is).Int) && ~isempty(AstS(Is).Wave))
                    WaveRange(Is) = range(AstS(Is).Wave);
                    if (nargout>1)
                        IntRange(Is) = range(AstS(Is).Int);
                        
                        if (nargout>2)
                            WaveRatio(Is) = max(AstS(Is).Wave)./min(AstS(Is).Wave);
                    
                            if (nargout>3)
                                IntRatio(Is) = max(AstS(Is).Int)./min(AstS(Is).Int);
                            end
                        end
                    end
                end
            end
        end
        
        function [MeanInt,MeanWave]=mean(AstS)
            % Description: Return the nanmean intensity and mean wavelength
            %              for each spectrum in an AstSpec object.
            % Input  : - An AstSpec class object.
            % Output : - A vector of nanmean intensity in each spectrum.
            %          - A vector of nanmean wavelength in each spectrum.
            % Example: AstS.mean
           
            
            Ns       = numel(AstS);
            SizeS    = size(AstS);
            MeanInt  = zeros(SizeS).*NaN;
            MeanWave = zeros(SizeS).*NaN;
            for Is=1:1:Ns
                if (~isempty(AstS(Is).Int) && ~isempty(AstS(Is).Wave))
                    MeanInt(Is) = nanmean(AstS(Is).Int);
                    if (nargout>1)
                        MeanWave(Is) = nanmean(AstS(Is).Wave);
                    end
                end
            end
            
        end
                    
        function [MedInt,MedWave]=median(AstS)
            % Description: Return the nanmedian intensity and nanmedian wavelength
            %              for each spectrum in an AstSpec object.
            % Input  : - An AstSpec class object.
            % Output : - A vector of nanmedian intensity in each spectrum.
            %          - A vector of nanmedian wavelength in each spectrum.
            % Example: AstS.median
           
            
            Ns       = numel(AstS);
            SizeS    = size(AstS);
            MedInt  = zeros(SizeS).*NaN;
            MedWave = zeros(SizeS).*NaN;
            for Is=1:1:Ns
                if (~isempty(AstS(Is).Int) && ~isempty(AstS(Is).Wave))
                    MedInt(Is) = nanmedian(AstS(Is).Int);
                    if (nargout>1)
                        MedWave(Is) = nanmedian(AstS(Is).Wave);
                    end
                end
            end
            
        end
               
        function [StdInt,StdWave]=std(AstS)
            % Description: Return the nanstd intensity and nanstd wavelength
            %              for each spectrum in an AstSpec object.
            % Input  : - An AstSpec class object.
            % Output : - A vector of nanstd intensity in each spectrum.
            %          - A vector of nanstd wavelength in each spectrum.
            % Example: AstS.std
           
            
            Ns       = numel(AstS);
            SizeS    = size(AstS);
            StdInt  = zeros(SizeS).*NaN;
            StdWave = zeros(SizeS).*NaN;
            for Is=1:1:Ns
                if (~isempty(AstS(Is).Int) && ~isempty(AstS(Is).Wave))
                    StdInt(Is) = nanstd(AstS(Is).Int);
                    if (nargout>1)
                        StdWave(Is) = nanstd(AstS(Is).Wave);
                    end
                end
            end
            
        end
        
        function [WeightedWave]=wwave(AstS)
            % Description: Return the weighted wavelength (weighted by
            %              intensity)
            %              for each spectrum in an AstSpec object.
            % Input  : - An AstSpec class object.
            % Output : - A vector of weighted wavelength in each spectrum.
            % Example: AstS.wwave
           
            
            Ns       = numel(AstS);
            SizeS    = size(AstS);
            WeightedWave  = zeros(SizeS).*NaN;
            for Is=1:1:Ns
                if (~isempty(AstS(Is).Int) && ~isempty(AstS(Is).Wave))
                    WeightedWave(Is) = nansum(AstS(Is).Int.*AstS(Is).Wave)./ ...
                                       nansum(AstS(Is).Int);
                    
                    
                end
            end
            
        end
    end
     
    % filters
    methods
        
        function AS=medfilt1(AS,varargin)
            % Description: Run median filter of order N on spectra
            %              points using the medfilt1.m function.
            % Input  : - AstSpec class object.
            %          * Additional arguments to pass to medfilt1.m.
            %            Note that the defualt order is 3.
            % Output : - AstSpec class object in which the .Int and .Back
            %            fields are filtered.
            %            The .Err field is divided by sqrt(N).
            % Example: A=AstSpec.get_galspec('QSO_NIR'); plot(medfilt1(A))
            
            if (isempty(varargin))
                N = 3;
            end
            Ns = numel(AS);
            for Is=1:1:Ns
                if (~isempty(AS(Is).Int))
                    AS(Is).Int = medfilt1(AS(Is).Int,varargin{:});
                end
                if (~isempty(AS(Is).Back))
                    AS(Is).Back = medfilt1(AS(Is).Back,varargin{:});
                end
                if (~isempty(AS(Is).Err))
                    AS(Is).Err = AS(Is).Err./sqrt(N);
                end
            end
        end
        
        function AS=sgolayfilt(AS,varargin)
            % Description: Run Savitzky-Golay (polynomial) smoothing filter
            %              on points using the sgolayfilt.m function.
            % Input  : - AstSpec class object.
            %          * Additional arguments to pass to sgolayfilt.m.
            %            Note that the defualt order is 3, and the default
            %            window size is 11 (window must be odd).
            % Output : - AstSpec class object in which the .Int and .Back
            %            fields are filtered.
            %            The .Err field is divided by sqrt(window-order).
            
            if (numel(varargin)==0)
                varargin{1} = 3;    % polynomial order
                varargin{2} = 11;   % window size
            elseif (numel(varargin)==1)
                varargin{2} = 11;   % window size
            else
                % do nothing
            end
            
            Ns = numel(AS);
            for Is=1:1:Ns
                if (~isempty(AS(Is).Int))
                    AS(Is).Int = sgolayfilt(AS(Is).Int,varargin{:});
                end
                if (~isempty(AS(Is).Back))
                    AS(Is).Back = sgolayfilt(AS(Is).Back,varargin{:});
                end
                if (~isempty(AS(Is).Err))
                    AS(Is).Err = AS(Is).Err./sqrt(varargin{2}-varargin{1});
                end
            end
        end
        
        function AS=hampel(AS,varargin)
            % Description: Run the hampel filter (Outlier removal via
            %              Hampel identifier) on points using the hampel.m
            %              function.
            % Input  : - AstSpec class object.
            %          * Additional arguments to pass to hampel.m.
            %            Note that the defualt half window is 3,
            %            and the default nsigma is 3.
            % Output : - AstSpec class object in which the .Int and .Back
            %            fields are filtered.
            %            The .Err field is unchanged.
            
            
            Ns = numel(AS);
            for Is=1:1:Ns
                if (~isempty(AS(Is).Int))
                    AS(Is).Int = hampel(AS(Is).Int,varargin{:});
                end
                if (~isempty(AS(Is).Back))
                    AS(Is).Back = hampel(AS(Is).Back,varargin{:});
                end
               
            end
        end
        
        function AS=convolve(AS,KernelType,KernelPar)
            % Description: Convolve a spectrum with a kernel.
            % Input  : - AstSpec class object.
            %          - A vector representing a kernel.
            %            Alternatively, a few built in options
            %            are possible (see construct_matched_filter1d.m):
            %            'tophat' - top hat filter.
            %            'gauss' - gaussian (default).
            %            'triangle' - triangle.
            %            The kernel is generated using
            %            construct_matched_filter1d.m. Note that the
            %            units are sampling rather than wavelength.
            %          - A vector of parameters to pass to the:
            %            construct_matched_filter1d.m function.
            %            Default is [5 2].
            % Output : - AstSpec class object in which the .Int and .Back
            %            fields are convolved with the kernel.
            %            The .Err field is unchanged.
            % Example: S=convolve(S,'twohat',[16 8]);
            
            if (nargin==1)
                KernelType = 'gauss';
                KernelPar  = [5 2];
            elseif (nargin==2)
                KernelPar  = [5 2];
            else
                % do nothing
            end
            
            Kernel = construct_matched_filter1d(KernelType,KernelPar,1);
            
            Ns = numel(AS);
            for Is=1:1:Ns
                if (~isempty(AS(Is).Int))
                    AS(Is).Int = conv(AS(Is).Int,Kernel,'same');
                end
                if (~isempty(AS(Is).Back))
                    AS(Is).Back = conv(AS(Is).Back,Kernel,'same');
                end
               
            end
        end
        
        function [FiltAstS,AstS]=filter_lines(AstS)
            % Description: Attempt to filter spectral lines from a spectrum
            %              using two step, two-top-hat median filter.
            %              Each spectrum is resampled logarithmically.
            %              A median filter with width of MedFiltVel1
            %              and MedFiltVel2 (measured in km/s) are
            %              calculated. For each median filter, the
            %              mean between points sperated by  TwoHatVel1 and
            %              TwoHatVel2 (km/s) are calculated.
            %              The best of the two is used to estimate the
            %              local continuum.
            % Input  : -
            % Output : - AstSpec class array in which the resampled spectra 
            %            are filtered for spectral lines (i.e., spectral 
            %            lines are removed).
            %          - The original AstSpec class object with the
            %            logarithmic resampling.
            % Example: [FiltAstS,AstS]=filter_lines(AstS)
            
            Samp = 'log';
            OverSamp = 1;
            MedFiltVel1 = 3000;
            TwoHatVel1 = 10000;
            MedFiltVel2 = 5000;
            TwoHatVel2 = 20000;
            
            c = constant.c./1e5; % [km/s]
            
            Nsigma    = 3;
            
            % resample the spectrum to a uniform logarithmic grid
            AstS = resample(AstS,Samp,OverSamp);
                
            FiltAstS = AstSpec(size(AstS));
            Ns = numel(AstS);
            for Is=1:1:Ns
                % for each spectrum
                
                FiltAstS(Is) = AstS(Is);   
                
                % resolution of spectrum
                Tmp = diff(AstS(Is).Wave)./AstS(Is).Wave(1:end-1);
                ResolutionZ = Tmp(1);
                ResolutionV = ResolutionZ.*c;  % [km/s]
                
                %--- Filtering with filter #1 ---
                Nmf1 = ceil(MedFiltVel1./ResolutionV);  % block size of median filter
                MedFilt1 = medfilt1(AstS(Is).Int,Nmf1);
                
                Nth1 = ceil(TwoHatVel1./ResolutionV);  % dist of two top hats
                % make sure Nth is even number
                if (~Util.array.is_evenint(Nth1))
                    Nth1 = Nth1 + 1;
                end
                Delta1 = MedFilt1(Nth1+1:end) - MedFilt1(1:end-Nth1);
                Delta1 = padarray(Delta1,Nth1.*0.5,'both');
                RStd1  = diff(Util.stat.err_cl(Delta1,0.68)).*0.5;  % robust std
                FlagBad1  = abs(Delta1)>Nsigma.*RStd1;  % flag for regions with bad background estimate
                
                Mean1 = 0.5.*(MedFilt1(Nth1+1:end) + MedFilt1(1:end-Nth1));
                Mean1  = padarray(Mean1,Nth1.*0.5,'both');
                
                
                %--- Filtering with filter #2 ---
                Nmf2 = ceil(MedFiltVel2./ResolutionV);  % block size of median filter
                MedFilt2 = medfilt1(AstS(Is).Int,Nmf2);
                
                Nth2 = ceil(TwoHatVel2./ResolutionV);  % dist of two top hats
                % make sure Nth is even number
                if (~Util.array.is_evenint(Nth2))
                    Nth2 = Nth2 + 1;
                end
                Delta2 = MedFilt2(Nth2+1:end) - MedFilt2(1:end-Nth2);
                Delta2 = padarray(Delta2,Nth2.*0.5,'both');
                RStd2  = diff(Util.stat.err_cl(Delta2,0.68)).*0.5;  % robust std
                FlagBad2  = abs(Delta2)>Nsigma.*RStd2;  % flag for regions with bad background estimate
                
                Mean2 = 0.5.*(MedFilt2(Nth2+1:end) + MedFilt2(1:end-Nth2));
                Mean2  = padarray(Mean2,Nth2.*0.5,'both');
                
                Mean   = Mean1;
                Flag2isBetter = abs(Delta1)>abs(Delta2);
                Mean(Flag2isBetter) = Mean2(Flag2isBetter);
                
                FiltAstS(Is).Int = Mean;
            end
            
        end
        
        function AstS=filter(AstS,Filter,Pars)
            % Filter spectra using an fft filter.
            % Description: Filter spectra using an fft filter.
            % Input  : - AstSpec class object.
            %          - Filter. Either a vector or one of the following
            %            options:
            %            'smooth' - Smooth spectra (high pass).
            %                       Paramaeter is the number of high
            %                       frequencies to remove.
            %            'flat'   - Flatten spectra (low pass).
            %                       Paramaeter is the number of low
            %                       frequencies to remove.
            %            'flat-smooth' - Smooth and flatten (band filter)
            %                       Parameters are the number of low,
            %                       and high frequencies to remove.
            %          - Filter parameters.
            % Output : - AstSpec class object filtered.
            % Example: AstS=filter(AstS,'flat-smooth',[10 10])
            
            error('Under construction');
            Ns = numel(AstS);
            for Is=1:1:Ns
                % for each spectrum
                % filter construction
                if (~isempty(AstS(Is).Wave))
                    Len = numel(AstS(Is).Wave);
                    'need to check!!!!'
                    if (ischar(Filter))
                        switch lower(Filter)
                            case 'smooth'
                                Filter = ones(Len,1);
                                Filter(1:Pars(1)) = 0;
                            case 'flat'
                                Filter = ones(Len,1);
                                Filter(end-Pars(1):end) = 0;
                            case 'flat-smooth'
                                Filter = ones(Len,1);
                                Filter(end-Pars(1):end) = 0;
                                Filter(1:Pars(2)) = 0;
                            otherwise
                                error('Unknwon Filter option');
                        end
                        
                    else
                        % filter is given by user
                    end
                    % apply filter
                    if (~isempty(AstS(Is).Int))
                        AstS(Is).Int = ifft(fft(AstS(Is).Int).*Filter);
                    end
                    if (~isempty(AstS(Is).Back))
                        AstS(Is).Back = ifft(fft(AstS(Is).Back).*Filter);
                    end
                end
            end
            
        end
        
        function AstS=stdfilt1(AstS,varargin)
            % Calculate the std filter of an AstSpec object
            % Description: Calculate the one dimensional standard deviation
            %              filter of the intensity field in an AstSpec
            %              object. The blck size is on the number of
            %              points.
            %              If you want the std due to the noise its
            %              recomended to first subtract a smoothed version
            %              of the spectra.
            %              To calculate the error on the mean divide by
            %              sqrt(FilterOrder), or use errfilt1.m.
            % Input  : - An AstSpec object.
            %          - Filter order. Default is 3.
            %            See stdfilt1.m for details.
            %          - BlockSize. Default is the length of the input vector.
            %            See stdfilt1.m for details.
            %          - If this parameter is provided then use quantile instaed
            %            of std, where the parameter specify the fraction of data
            %            to be in the returned range. For example 0.6834 is analog
            %            to one sigma.
            %            If two elements vector is given than these are the lower
            %            and upper quantile range.
            %            Default is empty matrix (i.e., use std).
            % Output : - An AstSpec object with the std filtered spectrum
            %            in the .Int field.
            % Example: A=stdfilt1(AstS);
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                AstS(Is).Int = timeseries.stdfilt1(AstS(Is).Int,varargin{:});
            end
            
        end
                
        % Calculate the error on the mean 
        function AstS=errfilt1(AstS,varargin)
            % Calculate the error on the mean (from std filtering) of an AstSpec object
            % Description: Calculate the one dimensional standard deviation
            %              filter of the intensity field in an AstSpec
            %              object divided by the sqrt(Order).
            %              The blck size is on the number of
            %              points.
            %              If you want the std due to the noise its
            %              recomended to first subtract a smoothed version
            %              of the spectra.
            %              To calculate the std use stdfilt1.m.
            % Input  : - An AstSpec object.
            %          - Filter order. Default is 5.
            %            See stdfilt1.m for details.
            %          - BlockSize. Default is the length of the input vector.
            %            See stdfilt1.m for details.
            %          - If this parameter is provided then use quantile instaed
            %            of std, where the parameter specify the fraction of data
            %            to be in the returned range. For example 0.6834 is analog
            %            to one sigma.
            %            If two elements vector is given than these are the lower
            %            and upper quantile range.
            %            Default is empty matrix (i.e., use std).
            % Output : - An AstSpec object with the std filtered spectrum
            %            in the .Int field.
            % Example: A=errfilt1(AstS);
            
            
            if (nargin<2)
                varargin{1} = 5;
            end
            AstS = stdfilt1(AstS,varargin{:});
            AstS = AstS./sqrt(varargin{1});
        end
        
    end

    % fft and related functions
    methods
        
        function [AstS,F]=fft(AstS,Field)
            % Calculate the fft of A field of an AstSpec class.
            % Description: Calculate the fft of A field of an AstSpec
            %              class.
            % Input  : - AstSpec class object.
            %          - Field for which to calculate the FFT.
            %            Default is 'Int'.
            % Output : - The AstSpec array with the fft of the requested
            %            field stored in the new object.
            %            The other fields are not changed.
            %          - Vector of fft of the last AstSpec element.
            
            if (nargin<2)
                Field = 'Int';
            end
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                if (~isempty(AstS(Is).(Field)))
                    AstS(Is).(Field) = fft(AstS(Is).(Field));
                end
            end
            F = AstS(Is).(Field);
            
        
        end
        
        function [AstS,F]=ifft(AstS,Field)
            % Calculate the ifft of A field of an AstSpec class.
            % Description: Calculate the ifft of A field of an AstSpec
            %              class.
            % Input  : - AstSpec class object.
            %          - Field for which to calculate the IFFT.
            %            Default is 'Int'.
            % Output : - The AstSpec array with the ifft of the requested
            %            field stored in the new object.
            %            The other fields are not changed.
            %          - Vector of ifft of the last AstSpec element.
            
            if (nargin<2)
                Field = 'Int';
            end
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                if (~isempty(AstS(Is).(Field)))
                    AstS(Is).(Field) = ifft(AstS(Is).(Field));
                end
            end
            F = AstS(Is).(Field);
        
        end
        
        function [AstS,F]=fftshift(AstS,Field)
            % Calculate the fftshift of A field of an AstSpec class.
            % Description: Calculate the fftshift of A field of an AstSpec
            %              class.
            % Input  : - AstSpec class object.
            %          - Field for which to calculate the fftshift.
            %            Default is 'Int'.
            % Output : - The AstSpec array with the fftshift of the requested
            %            field stored in the new object.
            %            The other fields are not changed.
            %          - Vector of fftshift of the last AstSpec element.
            
            if (nargin<2)
                Field = 'Int';
            end
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                if (~isempty(AstS(Is).(Field)))
                    AstS(Is).(Field) = fftshift(AstS(Is).(Field));
                end
            end
            F = AstS(Is).(Field);
        
        end
        
        function [AstS,F]=ifftshift(AstS,Field)
            % Calculate the ifftshift of A field of an AstSpec class.
            % Description: Calculate the ifftshift of A field of an AstSpec
            %              class.
            % Input  : - AstSpec class object.
            %          - Field for which to calculate the ifftshift.
            %            Default is 'Int'.
            % Output : - The AstSpec array with the ifftshift of the requested
            %            field stored in the new object.
            %            The other fields are not changed.
            %          - Vector of ifftshift of the last AstSpec element.
            
            if (nargin<2)
                Field = 'Int';
            end
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                if (~isempty(AstS(Is).(Field)))
                    AstS(Is).(Field) = ifftshift(AstS(Is).(Field));
                end
            end
            F = AstS(Is).(Field);
        
        end
        
        function [AstS,F]=conj(AstS,Field)
            % Calculate the conj of A field of an AstSpec class.
            % Description: Calculate the conj of A field of an AstSpec
            %              class.
            % Input  : - AstSpec class object.
            %          - Field for which to calculate the conj.
            %            Default is 'Int'.
            % Output : - The AstSpec array with the conj of the requested
            %            field stored in the new object.
            %            The other fields are not changed.
            %          - Vector of conj of the last AstSpec element.
            
            if (nargin<2)
                Field = 'Int';
            end
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                if (~isempty(AstS(Is).(Field)))
                    AstS(Is).(Field) = conj(AstS(Is).(Field));
                end
            end
            F = AstS(Is).(Field);
        
        end
        
        function [AstS,F]=abs(AstS,Field)
            % Calculate the abs of a field of an AstSpec class.
            % Description: Calculate the abs of a field of an AstSpec
            %              class.
            % Input  : - AstSpec class object.
            %          - Field for which to calculate the abs.
            %            Default is 'Int'.
            % Output : - The AstSpec array with the abs of the requested
            %            field stored in the new object.
            %            The other fields are not changed.
            %          - Vector of abs of the last AstSpec element.
            
            if (nargin<2)
                Field = 'Int';
            end
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                if (~isempty(AstS(Is).(Field)))
                    AstS(Is).(Field) = abs(AstS(Is).(Field));
                end
            end
            F = AstS(Is).(Field);
        
        end
        
        function [AstS,F]=real(AstS,Field)
            % Calculate the real of a field of an AstSpec class.
            % Description: Calculate the real of a field of an AstSpec
            %              class.
            % Input  : - AstSpec class object.
            %          - Field for which to calculate the real.
            %            Default is 'Int'.
            % Output : - The AstSpec array with the real of the requested
            %            field stored in the new object.
            %            The other fields are not changed.
            %          - Vector of real of the last AstSpec element.
            
            if (nargin<2)
                Field = 'Int';
            end
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                if (~isempty(AstS(Is).(Field)))
                    AstS(Is).(Field) = real(AstS(Is).(Field));
                end
            end
            F = AstS(Is).(Field);
        
        end
        
        function [AstS,F]=imag(AstS,Field)
            % Calculate the imag of a field of an AstSpec class.
            % Description: Calculate the imag of a field of an AstSpec
            %              class.
            % Input  : - AstSpec class object.
            %          - Field for which to calculate the imag.
            %            Default is 'Int'.
            % Output : - The AstSpec array with the imag of the requested
            %            field stored in the new object.
            %            The other fields are not changed.
            %          - Vector of imag of the last AstSpec element.
            
            if (nargin<2)
                Field = 'Int';
            end
            
            Ns = numel(AstS);
            for Is=1:1:Ns
                if (~isempty(AstS(Is).(Field)))
                    AstS(Is).(Field) = imag(AstS(Is).(Field));
                end
            end
            F = AstS(Is).(Field);
        
        end
        
    end
    
        %---------------------------
        %--- Intensity mataching ---
        %---------------------------
        
        
    % correlations
    methods
        % correlation between two AstSpec objects
        function [Corr,Noverlap,Par,ParErr]=corr(AstS1,AstS2,Type)
            % Calculate the correlation between the spectra in AstSpec object.
            % Description: Calculate the correlation between the spectra
            %              in two AstSpec class array.
            %              The spectra are first resampled to their
            %              overlap region, flux normalized using lscov.m
            %              and the correlation coef. is calculated.
            % Input  : - AstSpec class array.
            %          - AstSpec class array. Single spectra, or the same
            %            length as the first AstSpec object.
            %          - Correlation type. See corr.m for options.
            %            Default is 'Pearson'.
            % Output : - Correlation coef. per spectrum.
            %          - Number of overlap pixels per correlation.
            %          - Best fit scale factor between spectra.
            %            Multiply first to get the second.
            %          - Error in best fit scale factor between spectra.
            % Example: S=AstSpec.get_pickles('g','v');
            %          [C,N,P,Pe] = corr(S,S(1));
            
            if (nargin<3)
                Type = 'Pearson';
            end
            
            % equalize sampling
            [AstS1,AstS2]=equalize_sampling_overlap(AstS1,AstS2);
            
            Ns1 = numel(AstS1);
            Ns2 = numel(AstS2);
            Ns  = max(Ns1,Ns2);
            Corr     = zeros(Ns,1).*NaN;
            Noverlap = zeros(Ns,1).*NaN;
            Par      = zeros(Ns,1).*NaN;
            ParErr   = zeros(Ns,1).*NaN;
            for Is=1:1:Ns
                % for each spectrum
                Is1 = min(Is,Ns1);
                Is2 = min(Is,Ns2);
                if (~isempty(AstS1(Is1).Int) && ~isempty(AstS2(Is2).Int))
                    % calc errors
                    if (isempty(AstS1(Is1).Err))
                        if (isempty(AstS2(Is2).Err))
                            Err = ones(size(AstS1(Is1).Int));
                        else
                            Err = AstS2(Is2).Err;
                        end
                    else
                        if (isempty(AstS2(Is2).Err))
                            Err = AstS1(Is1).Err;
                        else
                            Err = sqrt(AstS1(Is1).Err.^2 + AstS2(Is2).Err.^2);
                        end
                    end

                    % remove NaNs
                    Fnn = ~isnan(AstS1(Is1).Int) & ~isnan(AstS2(Is2).Int);
                    % for flux level
                    [Par(Is),ParErr(Is)] = lscov(AstS1(Is1).Int(Fnn),AstS2(Is2).Int(Fnn),1./Err(Fnn).^2);
                    Corr(Is) = corr(AstS1(Is1).Int(Fnn).*Par(Is),AstS2(Is2).Int(Fnn),'type',Type);
                    Noverlap(Is) = sum(Fnn);
                end
            end
        end
             
        % rms between spectra and photometry
        function Res=compare_phot(AstS,Data)
            % Compare an AstSpec object with photometric observations
            % Package: @AstSpec
            % Description: Given an AstSpec object and photometric data
            %              points, for each AstSpec elemnt perform
            %              synthetic photometry and compare the photometric
            %              data points with the synthetic photometry. For
            %              each spectrum, the function returns its rms and
            %              chi^2 with the photometric data points.
            % Input  : - An AstSpec object
            %          - Photometric data points:
            %            Either a cell array of
            %            {Mag, MagErr, Family, Band, MagSys}
            %            or a structure array with these fields.
            % Outout : - A structure array with an element per spectrum
            %            and the following fields:
            %            .RMS
            %            .Chi2
            %            .Dof
            % Example: 
            
            if (iscell(Data))
                Mag    = [Data{:,1}]';
                MagErr = [Data{:,2}]';
                Family = Data(:,3);
                Band   = Data(:,4);
                MagSys = Data(:,5);
            elseif (isstruct(Data))
                Mag    = [Data.Mag]';
                MagErr = [Data.MagErr]';
                Family = {Data.Family}';
                Band   = {Data.Band}';
                MagSys = {Data.MagSys}';
            else
                error('Unkown Data format');
            end
            Nband = numel(Band);
                
            Ns   = numel(AstS);
            SynMag  = zeros(Nband,Ns);
            Flag    = zeros(Nband,Ns);
            EffW    = zeros(Nband,Ns);
            for Is=1:1:Ns
                % for each spectrum
                for Iband=1:1:Nband
                    [SynMag(Iband,Is),Flag(Iband,Is),EffW(Iband,Is)] = synphot(AstS(Is),Family{Iband},Band{Iband},MagSys{Iband});
                end
                
                Res(Is).RMS  = std(Mag - SynMag(:,Is),[],1);
                Res(Is).Chi2 = sum( ((Mag - SynMag(:,Is))./MagErr).^2 );
                Res(Is).Dof  = numel(Mag) - 1;
            end
        end
                
        function [SG,Mat]=grid_z_ext(S,Vecz,VecEbv,VecR)
            % Duplicate a spectrum into a grid of redshifts and extinctions
            % Package: @AstSpec
            % Description: Given an AstSpec object with a single element
            %              create multiple versions of the spectrum in a
            %              grid of redshifts and extinctions.
            % Input  : - A single element AstSpec object (restframe).
            %          - A vector of redshifts.
            %          - A vector of E_{B-V} extinctions.
            %          - A vector of selective extinction R_{V}.
            % Output : - An AstSpect object with spectra with redshifts and
            %            extinctions in thr requested grid.
            %          - A three column matrix with [z, E_{B-V}, R_{V}].
            %            where the line correspondinf to the element index
            %            in the first output.
            
            
            Nz   = numel(Vecz);
            NEbv = numel(VecEbv);
            NR   = numel(VecR);
            N    = Nz.*NEbv.*NR;
            SG   = AstSpec(N,1);
            Mat  = zeros(N,3);
            K    = 0;
            for Iz=1:1:Nz
                Sz    = shift(S,-Vecz(Iz),'f');
                for IEbv=1:1:NEbv
                    for IR=1:1:NR
                        K        = K + 1;
                        SG(K)    = extinction(Sz,VecEbv(IEbv),VecR(IR));
                        Mat(K,:) = [Vecz(Iz), VecEbv(IEbv),VecR(IR)];
                    end
                end
            end
                    
        end
        
%         function []=xcorr(AstS1,AstS2,VecZ,Type)
%             %
%           
%         end
        
        
        %% << got here
%         function [Z,XC]=xcorr(AstS1,AstS2)
%             % Description: 
%             
%             % equalize sampling
%             [AstS1,AstS2]=equalize_sampling_overlap(AstS1,AstS2,'log');
%             
%             Ns1 = numel(AstS1);
%             Ns2 = numel(AstS2);
%             Ns  = max(Ns1,Ns2);
%             Corr     = zeros(Ns,1).*NaN;
%             Noverlap = zeros(Ns,1).*NaN;
%             Par      = zeros(Ns,1).*NaN;
%             ParErr   = zeros(Ns,1).*NaN;    
%             for Is=1:1:Ns,
%                 % for each spectrum
%                 
%                 Is1 = min(Is,Ns1);
%                 Is2 = min(Is,Ns2);
%                 if (~isempty(AstS1(Is1).Int) && ~isempty(AstS2(Is2).Int)),
%                     % remove NaNs
%                     Fnn = ~isnan(AstS1(Is1).Int) & ~isnan(AstS2(Is2).Int);
%                     % for flux level
%                     
%                     [Par(Is),ParErr(Is)] = lscov(AstS1(Is1).Int(Fnn),AstS2(Is2).Int(Fnn),1./Err(Fnn).^2);
%                                         
%                     Res(Is).XC  = ifft(fft(AstS1(Is1).Int.*Par(Is)).*conj(fft(AstS2(Is2).Int)));
%                     Res(Is).Lag = [] %  %<<<---
%                     
%                     [MaxXC,MaxInd] = max(Res(Is).XC);
%                     MaxLag = Res(Is).Lag(MaxInd);
%                 end
%             end
%         end
                     
    end
    
    % model fitting
    methods
        
        function [Res,SpecBB]=fit_bbT(AstS,T,varargin)
            % Fit the normalization of black-body spectra with a known temperature
            % Package: @AstSpec
            % Description: Fit the normalization of black-body spectra with a known temperature
            % Input  : - An AstSpec object
            %          - Temperature. If array, then each element
            %            correspond to an AstSpec element.
            % Output : - Structure with best fit information and scaling
            %          - AstSpec object containing best fit black-body
            %            spectra.
            % Example: S=AstSpec.get_pickles('g','v');
            %          S3=filter_lines(S(3));
            %          [R,Sb]=fit_bbT(S3,5700)
            
           
            WaveField = 'Wave';
            IntField  = 'Int';
            
            Ns = numel(AstS);
            WaveUnits = get_wave_units(AstS);
            IntUnits  = get_int_units(AstS);
            for Is=1:1:Ns
                It = min(Is,numel(T));
                if (isempty(IntUnits{Is}))
                    IntU = 'cgs/A';
                else
                    IntU = IntUnits{Is};
                end
                SpecBB(Is) = AstSpec.blackbody(T(It),AstS(Is).(WaveField),IntU,WaveUnits{Is});
            
                Res(Is).Factor = SpecBB(Is).(IntField)\AstS(Is).(IntField);
                %Res(Is).Factor = mean(SpecBB(Is).(IntField)./AstS(Is).(IntField))
                Res(Is).AngRad = sqrt(Res(Is).Factor);
                Res(Is).Resid   = SpecBB(Is).(IntField).* Res(Is).Factor - AstS(Is).(IntField);
                
                Res(Is).RMS     = std(Res(Is).Resid);
                Res(Is).RelRMS  = std(Res(Is).Resid./(SpecBB(Is).(IntField).* Res(Is).Factor));
               
                SpecBB(Is).(IntField) = SpecBB(Is).(IntField).*Res(Is).Factor;
            end
            
        end
        
        function [Res,SpecBB]=fit_bb(AstS,varargin)
            % Fit a black-body spectrum to all spectra in an AstSpec object.
            % Description: Fit a black-body spectrum to all spectra in
            %              an AstSpec class object.
            % Input  : - AstSpec class object.
            %          * Additional arguments to pass to fit_bb.m
            %            See fit_bb.m for options.
            % Output : - Structure array of best fit results for each
            %            spectrum.
            %          - The best fit black-body spectrum for each
            %            input spectrum.
            % Example: S=AstSpec.get_pickles('g','v');
            %          S3=filter_lines(S(3));
            %          [R,Sb]=fit_bb(S3);
            
            Ns = numel(AstS);
            SpecBB = AstSpec(Ns,1);
            for Is=1:1:Ns
                % for each spectrum
                if (~isempty(AstS(Is).Wave) && ~isempty(AstS(Is).Int))
                    if (isempty(AstS(Is).Err))
                        % no errors provided
                        [Res(Is),SpecBB(Is)] = AstroUtil.spec.fit_bb([AstS(Is).Wave,AstS(Is).Int],varargin{:});
                    else
                        % error are provided 
                        [Res(Is),SpecBB(Is)] = AstroUtil.spec.fit_bb([AstS(Is).Wave,AstS(Is).Int,AstS(Is).Err],varargin{:});
                    end
                end
            end
        end
        
    end
    
    % wavelength shifts
    methods
            
        function AS=shift(AS,Z,MethodF)
            % redshift/blue shift spectrum
            % Description: redshift/blue shift spectrum
            % Input  : - AstSpec class object.
            %          - Redshift. If >0 then transform from observed
            %            frame to rest frame, and vice versa of <0.
            %          - Method for flux conversion:
            %            'w' - f_lambda: (1+z)^3
            %            'f' - f_nu: (1+z)
            %            'n' - do nothing
            % Output : AstSpec in which the wavelength is shifted
            %          and the .Int and .Err fields are corrected.
            % Example: A(1)=shift(A(1),0.1,'w')
            
            Ns = numel(AS);
            for Is=1:1:Ns
                % for each spectrum
                % update redshift
                if (~isempty(AS(Is).z))
                    if (Z>0)
                        AS(Is).z = AS(Is).z./(1+Z);
                    else
                        AS(Is).z = AS(Is).z.*(1+abs(Z));
                    end
                end
                
                % wavelength
                if (~isempty(AS(Is).Wave))
                    if (Z>0)
                        AS(Is).Wave = AS(Is).Wave./(1+Z);
                    else
                        AS(Is).Wave = AS(Is).Wave.*(1+abs(Z));
                    end
                end
                
                % intensity
                if (~isempty(AS(Is).Int))
                    switch lower(MethodF)
                        case 'w'
                            if (Z>0)
                                AS(Is).Int = AS(Is).Int.*(1+Z).^3;
                            else
                                AS(Is).Int = AS(Is).Int.*(1+Z).^-3;
                            end
                        case 'f'
                            if (Z>0)
                                AS(Is).Int = AS(Is).Int.*(1+Z);
                            else
                                AS(Is).Int = AS(Is).Int./(1+Z);
                            end
                        case 'n'
                            % do nothing
                        otherwise
                            error('Unknown MethodF option');
                    end
                end
                
                % error
                if (~isempty(AS(Is).Err))
                    switch lower(MethodF)
                        case 'w'
                            if (Z>0)
                                AS(Is).Err = AS(Is).Err.*(1+Z).^3;
                            else
                                AS(Is).Err = AS(Is).Err.*(1+Z).^-3;
                            end
                        case 'f'
                            if (Z>0)
                                AS(Is).Err = AS(Is).Err.*(1+Z);
                            else
                                AS(Is).Err = AS(Is).Err./(1+Z);
                            end
                        case 'n'
                            % do nothing
                        otherwise
                            error('Unknown MethodF option');
                    end
                end
                
            end
        end
        
        function AS=shift_vel(AS,Vel,MethodF)
            % redshift/blue shift spectrum given velocity.
            % Description: redshift/blue shift spectrum given velocity.
            % Input  : - AstSpec class object.
            %          - Recession velocity. If >0 then transform from
            %            observed frame to rest frame, and vice versa of <0.
            %          - Method for flux conversion:
            %            'w' - f_lambda: (1+z)^3
            %            'f' - f_nu: (1+z)
            %            'n' - do nothing
            % Output : AstSpec in which the wavelength is shifted
            %          and the .Int and .Err fields are corrected.
            % Example: AS=shift_vel(A,1000,'w')
            
            Z = AstroUtil.spec.vel2shift(Vel);
            
            AS=shift(AS,Z,MethodF);
        end
            
    end
    
      
    % spectral lines
    methods
        function AS=region_nan(AS,Ranges)
            % Set the .Int field in AstSpec object to NaN
            % Description: Set the .Int field in AstSpec object
            %              to NaN when their corresponding .Wave
            %              is in some specific wavelength ranges.
            % Input  : - AstSpec class object.
            %          - Two column matrix of ranges [Low, High].
            %            For each wavelength in each of the low-to-high
            %            ranges, the .Int will be set to NaN.
            %            If this is one of the following strings:
            %            'telluric' - set Ranges to []
            % Output : - AstSpec class object with the required NaN
            %            in the .Int field only.
            % Example: AS=region_nan(AS,[5000 5500]);
            %          plot(region_nan(A(1),'telluric'))
            
            if (ischar(Ranges))
                switch lower(Ranges)
                    case 'telluric'
                        Ranges = [7450 7620; 6461 6711];
                    otherwise
                        error('Unknown Ranges option');
                end
            end
            
            Nr = size(Ranges,1);
            Ns = numel(AS);
            for Is=1:1:Ns
                % for each spectrum
                Ind = Util.array.find_ranges(AS(Is).Wave,Ranges);
                AS(Is).Int(Ind) = NaN;
            end
        end
                
    end
    
       
    % syntheic photometry
    methods 
        function [Mag,Flag,EffW]=synphot(AS,varargin)
            % Synthetic photometry on AstSpec class spectra.
            % Description: Synthetic photometry on AstSpec class spectra.
            %              OBSOLETE: use synthetic_phot instead.
            % Input  : - AstSpec class object.
            %          * Additional arguments to pass to synphot.m
            % Output : - Vector of synthetic magnitude for each spectrum.
            %          - Vector of the fraction of flux that was
            %            extrapolated in case of partial coverage between
            %            spectrum and filter. 0 - means no extrapolation.
            %          - Vector of filter effective wavelength [Ang].
            % Example: [Mag,Flag,EffW]=synphot(AS,'SDSS','r','AB');

            Ns   = numel(AS);
            Mag  = zeros(Ns,1);
            Flag = zeros(Ns,1);
            EffW = zeros(Ns,1);
            % convert Wavelength to ang:
            AS = convert_wave(AS,'Ang');
            
            for Is=1:1:Ns
                % for each spectrum
                [Mag(Is),Flag(Is),EffW(Is)]=AstroUtil.spec.synphot([AS(Is).Wave, AS(Is).Int],varargin{:});
                
            end
        end
        
        
        function [Mag,Cover,Flux]=synthetic_phot(AS,varargin)
            % Synthetic photometry on AstSpec class spectra.
            % Input  : - AstSpec class object.
            %          * Additional arguments to pass to synthetic_phot.m
            % Output : - Vector of synthetic magnitude for each spectrum.
            %          - Vector of the fraction of flux that was
            %            extrapolated in case of partial coverage between
            %            spectrum and filter. 0 - means no extrapolation.
            %          - Vector of filter effective wavelength [Ang].
            % Example: [Mag,Flag,EffW]=synthetic_phot(AS,'SDSS','r','AB');
            
            Ns   = numel(AS);
            Mag  = zeros(Ns,1);
            Flag = zeros(Ns,1);
            EffW = zeros(Ns,1);
            % convert Wavelength to ang:
            AS = convert_wave(AS,'Ang');
            
            for Is=1:1:Ns
                % for each spectrum
                [Mag(Is),Flag(Is),EffW(Is)]=AstroUtil.spec.synthetic_phot([AS(Is).Wave, AS(Is).Int],varargin{:});
                
            end
            
            
        end
        
        
        function Spec=scale2mag(Spec,Mag,Family,Name,System)
            % Scale an AstCat object spectrum to have a specific synthetic magnitude
            % Package: @AstCat
            % Description: Scale an AstCat object spectrum to have a
            %              specific synthetic magnitude.
            % Input  : - An AstCat object.
            %          - Synthetic magnitude of the output spectrum.
            %          - Filter family. Default is 'SDSS'.
            %          - Filter name. Default is 'r'.
            %          - Mag system. Default is 'AB'.
            % Output : - An AstCat magnitude with scaled flux.
            
            if (nargin<5)
                System = 'AB';
                if (nargin<4)
                    Name = 'r';
                    if (nargin<3)
                        Family = 'SDSS';
                    end
                end
            end
            
            Ns = numel(Spec);
            for Is=1:1:Ns
                Mag0 = synphot(Spec(Is),Family,Name,System);
                Spec(Is).Int = Spec(Is).Int.*10.^(-0.4.*(Mag-Mag0));
            end
            
            
        end
        
    end
    
    
        
    % Interpolation, resampling, stretching, normalization
    methods
        function AstS=interp(AstS,W,varargin)
            % Interpolate an astronomical spectra class into a new wavelngth grid.
            % Description: Interpolate an astronomical spectra class
            %              into a new wavelngth grid.
            % Input  : - AstSpec class.
            %          - Column vector of wavelength grid, in the same
            %            units as the wavelength in the spectra class.
            % Output : - AstSpec class with a new wavelength grid.
            %            The .Wave, .Int, .Err, .Back, .Mask and .AddCol
            %            fields are resampled.
            % Example: interp(A,[5000:1:5100]);


            if (nargin<2)
                error('Wavelength column vector must be provided');
            end

            Ns = numel(AstS);

            for Is=1:1:Ns
                if (isempty(AstS(Is).Wave) || isempty(AstS(Is).Int))
                    warning('AstSpec number %d is empty - not interpolating',Is);
                else
                    AstS(Is).Int  = interp1(AstS(Is).Wave, AstS(Is).Int, W, varargin{:});
                    if (~isempty(AstS(Is).Err))
                        AstS(Is).Err  = interp1(AstS(Is).Wave, AstS(Is).Err, W, varargin{:});
                    end
                    if (~isempty(AstS(Is).Back))
                        AstS(Is).Back  = interp1(AstS(Is).Wave, AstS(Is).Back, W, varargin{:});
                    end
                    if (~isempty(AstS(Is).Mask))
                        AstS(Is).Mask  = interp1(AstS(Is).Wave, AstS(Is).Mask, W, 'nearest');
                    end
                    if (~isempty(AstS(Is).AddCol))
                        AstS(Is).AddCol  = interp1(AstS(Is).Wave, AstS(Is).AddCol, W, varargin{:});
                    end
                    % only now populate also Wavelength
                    AstS(Is).Wave = W;

                end
            end
        end

        function AstS=resample(AstS,Samp,OverSamp)
            % Description: Resample a spectra in AstSpec class linearly
            %              or logarithmically.
            % Input  : - AstSpec class object.
            %          - Sampling type: 'linear' | 'log'. Default is
            %            'linear'.
            %          - Over sampling. Default is 1.
            % Output : - AstSpec class object resampled.
            % Example: S=resample(S,'log')

            if (nargin<2)
                Samp = 'linear';
                OverSamp = 1;
            elseif (nargin<3)
                OverSamp = 1;
            else
                % do nothing
            end

            N = numel(AstS);
            for I=1:1:N
                if (isempty(AstS(I).Wave) || isempty(AstS(I).Int))
                    warning('AstSpec number %d is empty - not interpolating',I);
                else
                    % define the range between the objects
                    StartRange = min(AstS(I).Wave);
                    EndRange   = max(AstS(I).Wave);
                    Step       = min(diff(AstS(I).Wave));
                    % linear or log sampling
                    switch lower(Samp)
                        case 'linear'
                            WaveVec       = (StartRange:Step./OverSamp:EndRange).';
                        case 'log'
                            Nlog          = ceil(OverSamp.*(EndRange-StartRange)./Step);
                            WaveVec       = logspace(log10(StartRange),log10(EndRange),Nlog).';
                        otherwise
                            error('Illegal Samp option');
                    end
                    if (isempty(WaveVec))
                        warning('AstSpec number %d - no overlap',I);
                    else
                        AstS(I) = interp(AstS(I),WaveVec);

                    end
                end
            end



        end

        function AstS=equalize_sampling(AstS1,AstS2)
            % Description: Equalize the wavelength sampling of of two
            %              astronomical spectra class objects.
            %              Change the sampling of the first object to
            %              be equal to that of the second object.
            % Input  : - AstSpec class 
            %          - AstSpec class
            % Output : - The first AstSpec class resampled.
            %            The .Wave, .Int, .Err, .Back, .Mask and .AddCol
            %            fields are resampled.


            N1 = numel(AstS1);
            N2 = numel(AstS2);
            N  = max(N1,N2);
            AstS = AstS1;
            for I=1:1:N
                Is1 = min(I,N1);
                Is2 = min(I,N2);
                if (isempty(AstS1(Is1).Wave) || isempty(AstS1(Is1).Int) || ...
                    isempty(AstS2(Is2).Wave) || isempty(AstS2(Is2).Int) )
                    warning('AstSpec number %d is empty - not interpolating',I);
                else
                    AstS(I) = interp(AstS1(Is1),AstS2(Is2).Wave);

                end

            end
        end

        function [AstS1,AstS2]=equalize_sampling_overlap(AstS1,AstS2,Samp,OverSamp)
            % Description: Equalize the wavelength sampling of of two
            %              astronomical spectra class objects.
            %              Change the sampling of the first object to
            %              be equal to that of the overlap range between
            %              the two objects.
            %              The resampled spectra have uniform sampling
            %              in linear space or log space.
            % Input  : - AstSpec class. 
            %          - AstSpec class.
            %          - Samling is uniform in 'linear' space or 'log'
            %            space. Default is 'linear'.
            %          - Oversampling factor. Default is 1.
            % Output : - The first AstSpec class resampled.
            %            The .Wave, .Int, .Err, .Back, .Mask and .AddCol
            %            fields are resampled.
            %          - The second AstSpec class resampled.
            % Example: [AstS1,AstS2]=equalize_sampling_overlap(AstS1,AstS2,'log')

            if (nargin<3)
                Samp = 'linear';
                OverSamp = 1;
            elseif (nargin<4)
                OverSamp = 1;
            else
                % do nothing
            end

            N1 = numel(AstS1);
            N2 = numel(AstS2);
            N  = max(N1,N2);
            for I=1:1:N
                Is1 = min(I,N1);
                Is2 = min(I,N2);
                if (isempty(AstS1(Is1).Wave) || isempty(AstS1(Is1).Int) || ...
                    isempty(AstS2(Is2).Wave) || isempty(AstS2(Is2).Int) )
                    warning('AstSpec number %d is empty - not interpolating',I);
                else
                    % define the range between the objects
                    StartRange = max(min(AstS1(Is1).Wave), min(AstS2(Is2).Wave));
                    EndRange   = min(max(AstS1(Is1).Wave), max(AstS2(Is2).Wave));
                    Step       = min(min(diff(AstS1(Is1).Wave)), min(diff(AstS2(Is2).Wave)));
                    % linear or log sampling
                    switch lower(Samp)
                        case 'linear'
                            WaveVec       = (StartRange:Step./OverSamp:EndRange).';
                        case 'log'
                            Nlog          = ceil(OverSamp.*(EndRange-StartRange)./Step);
                            WaveVec       = logspace(log10(StartRange),log10(EndRange),Nlog).';
                        otherwise
                            error('Illegal Samp option');
                    end
                    if (isempty(WaveVec))
                        warning('AstSpec number %d - no overlap',I);
                    else
                        AstS1(Is1) = interp(AstS1(Is1),WaveVec);
                        AstS2(Is2) = interp(AstS2(Is2),WaveVec);
                    end
                end
            end

        end

        function AstS=norm(AstS,Norm)
            % Description: normalize the integral of the .Int field to be
            %              be some number. Also scales the .Err and .Back
            %              fields.
            % Input  : - AstSpec class.
            %          - Optional normalization. Default is 1.
            % Output : - AstSpec class in which the .Int, .Err and .Back
            %            fields are normalized to have an integral equal
            %            Norm.

            if (nargin==1)
                Norm = 1;
            end
            Ns = numel(AstS);
            for Is=1:1:Ns
                if (~isempty(AstS(Is).Wave) && ~isempty(AstS(Is).Int))
                    Integral = trapz(AstS(Is).Wave,AstS(Is).Int);
                    AstS(Is).Int = AstS(Is).Int.*Norm./Integral;
                end
                if (~isempty(AstS(Is).Err))
                    AstS(Is).Err = AstS(Is).Err.*Norm./Integral;
                end
                if (~isempty(AstS(Is).Back))
                    AstS(Is).Back = AstS(Is).Back.*Norm./Integral;
                end
            end
        end

        function [AstS,Factor]=stretch(AstS,Factor)
            % Description: Stretch spectra intensity in an AstSpec class
            % Input  : - AstSpec class object.
            %          - Stretch factor or a function handle
            %            to use for the stretch calculation.
            %            This can be a scalar or vector (element per
            %            spectrum).
            % Output : - AstSpec class object, in which the .Int,
            %            .Back, and .Err fields are multiplies by
            %            a stretch factor.
            %          - The stretch factor.
            % Example: AstS=stretch(A(2:3),2);
            %          [A1,F]=stretch(A(2),@mean)

            Ns = numel(AstS);
            if (isnumeric(Factor))
                Factor = Factor.*ones(size(AstS));
            elseif (isa(Factor,'function_handle'))
                Factor = Factor(AstS);
            else
                error('Unknwon Factor option');
            end


            for Is=1:1:Ns
                % for each spectrum
                AstS(Is).Int  = AstS(Is).Int.*Factor(Is);
                AstS(Is).Back = AstS(Is).Back.*Factor(Is);
                AstS(Is).Err  = AstS(Is).Err.*Factor(Is);
                AstS(Is).IntUnits = sprintf('%s * %f',AstS(Is).IntUnits,Factor(Is));
            end
        end

        function [AstS,Factor]=addcon(AstS,Factor)
            % Description: Add a value to spectra intensity in an AstSpec 
            % Input  : - AstSpec class object.
            %          - A value to add/subtract or a function handle
            %            to use for subtraction calculation.
            %            This can be a scalar or vector (element per
            %            spectrum).
            % Output : - AstSpec class object, in which the .Int,
            %            and .Back fields are added to a constant.
            %          - The addition factor (minus sign for function
            %            handle).
            % Example: AstS=addcon(A(2:3),2);
            %          [A1,F]=addcon(A(2),@mean)

            Ns = numel(AstS);
            if (isnumeric(Factor))
                Factor = Factor.*ones(size(AstS));
            elseif (isa(Factor,'function_handle'))
                Factor = -Factor(AstS);
            else
                error('Unknwon Factor option');
            end


            for Is=1:1:Ns
                % for each spectrum
                AstS(Is).Int  = AstS(Is).Int + Factor(Is);
                AstS(Is).Back = AstS(Is).Back + Factor(Is);
                %AstS(Is).Err  = AstS(Is).Err + Factor(Is);
                AstS(Is).IntUnits = sprintf('%s + %f',AstS(Is).IntUnits,Factor(Is));
            end
        end

    end
    
    % apply extinction curves
    methods
        
        function AstS=atmospheric_extinction(AstS,AirMass,Name)
            % Apply atmospheric extinction to an AstSpec object
            % Description: Get atmospheric extinction curve and apply it an
            %              an AstSpec object.
            %              If the airmass is positive than apply the
            %              extinction to the spectra. If negative, then
            %              assume that the spectra are extincted and
            %              convert it to a non-extincted spectra.
            % Input  : - An AstSpec object.
            %          - Airmass to apply. Default is 1.
            %            If the airmass is positive than apply the
            %            extinction to the spectra. If negative, then
            %            assume that the spectra are extincted and
            %            convert it to a non-extincted spectra.
            %          - An atmospheric extinction file name, or an AstSpec
            %            object containing the atmospheric extinction, or a
            %            matrix with the atmospheric extinction
            %            [Wave(ang) Ext(mag/airmass)].
            %            See AstSpec.get_atmospheric_extinction for
            %            details.
            %            Default is 'VLT'.
            % Output : - An AstSpec object with the extincted/corrected
            %            spectra.
            % Example: S=AstSpec.get_pickles('g','v');
            %          atmospheric_extinction(S,-2,'VLT')
            % Reliable: 2
            
            Def.Name    = 'VLT';
            Def.AirMass = 1;
            if (nargin==1)
                AirMass = Def.AirMass;
                Name    = Def.Name;
            elseif (nargin==2)
                Name    = Def.Name;
            elseif (nargin==3)
                % do nothing
            else
                error('Illegal number of input arguments: AstS=AstSpec.atmospheric_extinction(AstS,Name,AirMass)');
            end
            
            if (ischar(Name))
                Ext = AstSpec.get_atmospheric_extinction(Name);
            elseif (AstSpec.isastspec(Name))
                Ext = Name;
            elseif (isnumeric(Name))
                Ext = AstSpec.array2astspec(Name);
            else
                error('Unknown atmospheric extinction type');
            end
            
            %Ext.Wave(end+1)=12000;
            %Ext.Int(end+1) = 0.02;
            
            Ns = numel(AstS);
            AirMass = AirMass(:).*ones(Ns,1);
            AstS    = convert_wave(AstS,'Ang');
            for Is=1:1:Ns
                ExtI = interp(Ext,AstS(Is).Wave);
                AstS(Is).Int = AstS(Is).Int.*10.^(-0.4.*ExtI.Int.*AirMass(Is));
            end
               
            
        end
        
        function AS=extinction(AS,Ebv,R)
            % Applay extinction to AstSpec class
            % Description: Applay extinction to AstSpec class
            % Input  : - AstSpec class object.
            %          - Vector or scalar E_{B-V} [mag].
            %            If >0 then extinct the spectra,
            %            If <0 then correct for extinction.
            %            If scalar, then apply the same extinction to
            %            all spectra. If vector, apply each elemnt to
            %            to one spectrum.
            %          - R_{V}. Default is 3.08.
            % Output : - AstSpec for which extinction is applied to
            %            the .Int, .Back and .Err fields.
            % Example: S=AstSpec.get_pickles('g','v');
            %          plot(extinction(S(1),0.2,3.08)); hold on; plot(S(1))
            
            if (nargin<3)
                R = 3.08;
            end
            
            Ns = numel(AS);
            Ne = numel(Ebv);
            if (Ne>1 && Ne~=Ns)
                error('Numeber of extinction element does not corresponds to number of spectra');
            end
            Ebv = Ebv.*ones(Ns,1);
            
            WaveUnits = get_wave_units(AS);
            
            for Is=1:1:Ns
                if (~isempty(WaveUnits{Is}))
                    Conv = convert.units(lower(WaveUnits{Is}),'micron',1);
                else
                    Conv = convert.units('ang','micron',1);
                    warning('Assume wavelength for AstSpec %d is in ang',Is);
                end
                Ext = AstroUtil.spec.extinction(abs(Ebv(Is)),AS(Is).Wave.*Conv,[],R);
                if (~isempty(AS(Is).Int))
                    % apply extinction to .Int
                    AS(Is).Int = AS(Is).Int.*10.^(-sign(Ebv(Is)).*0.4.*Ext);
                end
                
                if (~isempty(AS(Is).Err))
                    % apply extinction to .Err
                    AS(Is).Err = AS(Is).Err.*10.^(-sign(Ebv(Is)).*0.4.*Ext);
                end
                
                if (~isempty(AS(Is).Back))
                    % apply extinction to .Back
                    AS(Is).Back = AS(Is).Back.*10.^(-sign(Ebv(Is)).*0.4.*Ext);
                end
            end        
        end

    end
       
    % General operator functions
    methods
        function AstS=astspec_fun1(AstS,Field,Fun,varargin)
            % Description: Run a unary function on one/many fields of
            %              an AstSpec class object, and write the result
            %              to this field.
            % Input  : - An AstSpec class object.
            %          - Field name
            %            {'Wave','Int','Err','Back','Mask','AddCol'}
            %            or a cell array of fields, or 'all', or 'notwave'.
            %            If empty use 'Int'.
            %          - Function handle (e.g., @sin).
            %          * Additional arguments to pass to the function.
            % Output : - An AstSpec class object in which the field
            %            is modified by the function operation.
            % Example: AstS=astspec_fun1(S,'all',@sin);
            
            if (isempty(Field))
                Field = 'Int';
            end
            if (~iscell(Field))
                if (strcmp(Field,'all'))
                    % all fields
                    Field = {'Wave','Int','Err','Back','Mask','AddCol'};
                elseif (strcmp(Field,'notwave'))
                    % all fields except Wave
                    Field = {'Int','Err','Back','Mask','AddCol'};
                else
                    Field = {Field};
                end
            end
            
            Nf = numel(Field);
            Ns = numel(AstS);
            for Is=1:1:Ns
                for If=1:1:Nf
                    if (~isempty(AstS(Is).(Field{If})))
                        AstS(Is).(Field{If}) = Fun(AstS(Is).(Field{If}),varargin{:});
                    end
                end
            end
        end
        
        function AstS=astspec_fun2(AstS1,AstS2,Field,Fun,varargin)
            % Description: Run a binary function on one/many fields of two
            %              AstSpec class object, and write the result
            %              to this field. The rest of the fields will
            %              be copied from the first AstSpec input.
            %              The wavelength grid of the two AstSpec object
            %              will be equalized using
            %              equalize_sampling_overlap.m
            % Input  : - First AstSpec class object.
            %          - Second AstSpec class object, or a a vector,
            %            or a scalar.
            %          - Field, or a cell array of fields:
            %            {'Wave','Int','Err','Back','Mask','AddCol'}
            %            or a cell array of fields, or 'all', or 'notwave'.
            %            If empty use 'Int'.
            %          - Function handle (e.g., @sin).
            %          * Additional arguments to pass to the function.
            % Output : - An AstSpec class object in which the field
            %            is modified by the function operation.
            % Example: AstS=astspec_fun2(S,S(1),'Int',@plus);
            %          AstS=astspec_fun2(AS1,2,'Int',@times);
            
            if (isempty(Field))
                Field = 'Int';
            end
            if (~iscell(Field))
                if (strcmp(Field,'all'))
                    % all fields
                    Field = {'Wave','Int','Err','Back','Mask','AddCol'};
                elseif (strcmp(Field,'notwave'))
                    % all fields except Wave
                    Field = {'Int','Err','Back','Mask','AddCol'};
                else
                    Field = {Field};
                end
            end
            
            Nf = numel(Field);
            Ns1 = numel(AstS1);
            if (~AstSpec.isastspec(AstS2))
                Ns2 = 1;
                % assume that the sampling is the same
            else
                Ns2 = numel(AstS2);
                % equalize the sampling to the overlap region
                [AstS1,AstS2]=equalize_sampling_overlap(AstS1,AstS2);
            end
            Ns  = max(Ns1,Ns2);
            AstS = AstS1;
            for Is=1:1:Ns
                Is1 = min(Is,Ns1);
                Is2 = min(Is,Ns2);
                for If=1:1:Nf
                    if (AstSpec.isastspec(AstS2))
                        % binary operation between two AstSpec objects
                        if (~isempty(AstS1(Is1).(Field{If})) && ~isempty(AstS2(Is2).(Field{If}))),
                            AstS(Is).(Field{If}) = Fun(AstS1(Is1).(Field{If}),AstS2(Is2).(Field{If}),varargin{:});
                        end
                    else
                        % second object is a scalar or a vector
                        AstS(Is).(Field{If}) = Fun(AstS1(Is1).(Field{If}),AstS2,varargin{:});
                    end
                end
                   
            end
        end
        
        function Vec=vec_fun1(AstS,Field,Fun,varargin)
            % Description: Run a unary function on a single field of
            %              an AstSpec class object, and return the result
            %              in a vector.
            % Input  : - An AstSpec class object.
            %          - Field name:
            %            'Wave','Int','Err','Back','Mask','AddCol'
            %            If empty then set to 'Int'.
            %          - Function handle (e.g., @sin).
            %          * Additional arguments to pass to the function.
            % Output : - An AstSpec class object in which the field
            %            is modified by the function operation.
            % Example: Vec=vec_fun1(AS,'Int',@sin);
            
            
            Ns = numel(AstS);
            if (Ns>1)
                error('vec_fun1 can work on a single AstSpec element, use astspec_fun1 for multiple elements');
            end
            Vec = Fun(AstS.(Field),varargin{:});
            
        end
        
        function Vec=vec_fun2(AstS1,Field1,AstS2,Field2,Fun,varargin)
            % Description: Run a binary function between one field in
            %              one AstSpec and another field in the second
            %              AstSpec, and return the result in a vector.
            %              The wavelength grid of the two AstSpec object
            %              will be equalized using
            %              equalize_sampling_overlap.m
            % Input  : - First AstSpec class object.
            %          - Field name in first AstSpec object.
            %            'Wave','Int','Err','Back','Mask','AddCol'
            %            If empty then set to 'Int'.
            %          - Second AstSpec class object, or a vector or a
            %            scalar.
            %          - Field name in second AstSpec object.
            %            'Wave','Int','Err','Back','Mask','AddCol'
            %            If empty then set to 'Int'.
            %          - Function handle (e.g., @plus).
            %          * Additional arguments to pass to the function.
            % Output : - An AstSpec class object in which the field
            %            is modified by the function operation.
            % Example: Vec=vec_fun2(S(1),'Int',S(2),'Wave',@plus);
            %          Vec=vec_fun2(AS1,'Int',2,[],@times);
            
            

            % equalize the sampling to the overlap region
            if (AstSpec.isastspec(AstS2))
                [AstS1,AstS2]=equalize_sampling_overlap(AstS1,AstS2);
                Ns = max(numel(AstS1),numel(AstS2));
            else
                Ns = numel(AstS1);
            end
            
            if (Ns>1)
                error('vec_fun2 can work on a single AstSpec element, use astspec_fun1 for multiple elements');
            end
            if (AstSpec.isastspec(AstS2))
                if (~isempty(AstS1.(Field1)) && ~isempty(AstS2.(Field2)))
                    Vec = Fun(AstS1.(Field1),AstS2.(Field2),varargin{:});
                end
            else
                if (~isempty(AstS1.(Field1)))
                    Vec = Fun(AstS1.(Field1),AstS2,varargin{:});
                end
            end
                    
            
        end
    
        function [FunInt,FunWave]=fun_scalar(AstS,Fun,varargin)
            % Description: Run a function that get a vector and return a
            %              a scalar, on the wavelength and intensity
            %              columns of each spectra in an AstSpect object.
            % Input  : - An AstSpec class object.
            %          - A function handle (e.g., @mean).
            %            The function must return a scalar value.
            %          * Additional arguments to apss to the function.
            % Output : - A vector of function evaluations on the intensity
            %            column in each spectrum.
            %          - A vector of function evaluations on the wavelength
            %            column in each spectrum.
            % Example: AstS.fun_scalar(@mean)
           
            
            Ns       = numel(AstS);
            SizeS    = size(AstS);
            FunInt  = zeros(SizeS).*NaN;
            FunWave = zeros(SizeS).*NaN;
            for Is=1:1:Ns
                if (~isempty(AstS(Is).Int) && ~isempty(AstS(Is).Wave))
                    FunInt(Is) = Fun(AstS(Is).Int,varargin{:});
                    if (nargout>1)
                        FunWave(Is) = Fun(AstS(Is).Wave,varargin{:});
                    end
                end
            end
            
        end

    end
    
    % Specific Operators
    methods
        function [AS]=astspec_arith(AS1,AS2,Operator)
            %--------------------------------------------------------------------------
            % astspec_arith function                                         AstroSpec
            % Description: Basic arithmetics on AstSpec class objects.
            %              Binary operations on two AstSpec class objects, or
            %              a AstSpec class object and a scalar.
            % Input  : - AstSpec class object.
            %          - AstSpec class object, a scalar or a vector.
            %            If vector, then the vector size should be equal to the
            %            size of the first argument.
            %          - Operator. E.g, @plus.
            % Output : - The result in an AstSpec class object.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Feb 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: B = astspec_arith(S(1:2),S(1),@times);
            % Reliable: 2
            %--------------------------------------------------------------------------

%             N1 = numel(AS1);
%             N2 = numel(AS2);
%             if (~(N2==1 || N2==N1))
%                 error('Size of second argument should be 1 or equal to that of the first argument');
%             end
            S1 = size(AS1);
            S2 = size(AS2);
            [S, S1, S2] = Util.array.bsx_size(S1,S2);
            
%             ndim = length(S);
%             pow1 = [1 cumprod(S1(1:end-1))];
%             pow2 = [1 cumprod(S2(1:end-1))];
            
            AS = AstSpec(S);
            if (~AstSpec.isastspec(AS2))
                % assume AS2 is numeric
%                 AS2 = AS2.*ones(size(AS1));
                
%                 for I1=1:1:N1
%                     I2 = I1;
                for I=1:prod(S)
                    sub = Util.array.ind2sub_array(S,I);
                    sub1 = sub;
                    sub1(S1==1) = 1;
                    I1 = Util.array.sub2ind_array(S1,sub1);
                    sub2 = sub;
                    sub2(S2==1) = 1;
                    I2 = Util.array.sub2ind_array(S2,sub2);
                        
                    % copy from AS1 fields which are not affected by the
                    % operation.
                    AS(I).Wave      = AS1(I1).Wave;
                    AS(I).Mask      = AS1(I1).Mask;
                    AS(I).WaveUnits = AS1(I1).WaveUnits;
                    AS(I).IntUnits  = AS1(I1).IntUnits;
                    AS(I).AddCol    = AS1(I1).AddCol;
                    AS(I).ObjName   = AS1(I1).ObjName;
                    AS(I).comments  = AS1(I1).comments;
                    AS(I).source    = AS1(I1).source;
                    AS(I).FileName  = AS1(I1).FileName;
                    AS(I).z         = AS1(I1).z;

                    if (~isempty(AS1(I1).Int))
                        AS(I).Int  = Operator(AS1(I1).Int,AS2(I2));
                    end
                    if (~isempty(AS1(I1).Back))
                        AS(I).Back = Operator(AS1(I1).Back,AS2(I2));
                    end
                    if (~isempty(AS1(I1).Err))
                        if (any(strcmpi(func2str(Operator),{'times','rdivide'})))
                            AS(I).Err  = Operator(AS1(I1).Err,AS2(I2));
                        else
                            % do nothing to .Err
                        end
                    end
                end
            else
                % AS2 is AstSpec
%                 N = max(N1,N2);
%                 for I=1:1:N
%                     I1 = min(I,N1);
%                     I2 = min(I,N2);
                for I=1:prod(S)
                    sub = Util.array.ind2sub_array(S,I);
                    sub1 = sub;
                    sub1(S1==1) = 1;
                    I1 = Util.array.sub2ind_array(S1,sub1);
                    sub2 = sub;
                    sub2(S2==1) = 1;
                    I2 = Util.array.sub2ind_array(S2,sub2);
                    
                    if (AS1(I1).Wave ~= AS2(I2).Wave)
                        error(AstSpec:astspec_arith:WaveNotAgreed,...
                            'On AstSpec arithmetic operation both AstSpec waves should be agreed.');
                    end
                    % copy from AS1 fields which are not affected by the
                    % operation.
                    AS(I).Wave      = AS1(I1).Wave;
                    AS(I).Mask      = AS1(I1).Mask;
                    AS(I).WaveUnits = AS1(I1).WaveUnits;
                    AS(I).IntUnits  = AS1(I1).IntUnits;
                    AS(I).AddCol    = AS1(I1).AddCol;
                    AS(I).ObjName   = AS1(I1).ObjName;
                    AS(I).comments  = AS1(I1).comments;
                    AS(I).source    = AS1(I1).source;
                    AS(I).FileName  = AS1(I1).FileName;
                    AS(I).z         = AS1(I1).z;

                    if (~isempty(AS1(I1).Int))
                        AS(I).Int = Operator(AS1(I1).Int,AS2(I2).Int);
                    end
                    if (~isempty(AS1(I1).Back))
                        if (~isempty(AS2(I2).Back))
                            AS(I).Back = Operator(AS1(I1).Back,AS2(I2).Back);
                        else
                            % Back of AS2 is not available - use AS1 back
                            % do nothing
                        end
                    end
                    if (~isempty(AS1(I1).Err))
                        if (~isempty(AS2(I2).Err))
                            if (any(strcmpi(func2str(Operator),{'times'})))
                                [~,AS(I).Err] = times_err(AS1(I1).Int,AS1(I1).Err, AS2(I2).Int,AS2(I2).Err);
                            elseif (any(strcmpi(func2str(Operator),{'rdivide'})))
                                [~,AS(I).Err] = rdivide_err(AS1(I1).Int,AS1(I1).Err, AS2(I2).Int,AS2(I2).Err);
                            elseif (any(strcmpi(func2str(Operator),{'plus','minus'})))
                                AS(I).Err = sqrt(AS1(I1).Err.^2 + AS2(I2).Err.^2);
                            else
                                % do nothing 
                            end
                        else
                            % Err of AS2 is not available - use AS1 Err
                            % do nothing
                        end
                    end
                end

            end
        end
       
        function AstS1=plus(AstS1,AstS2)
            % Description: Add AstSpec arrays (+).
            %              See astspec_arith for details.
            % Input   : - AstSpec array
            %           - AstSpec array, scalar or vector.
            % Outoput : - AstSpec class in which the sum of the
            %             .Int, .Back fields are given, an the
            %             error is propogated into .Err
            % Example : S+S(1)
            
            AstS1 = astspec_arith(AstS1,AstS2,@plus);
            
        end
            
        function AstS1=minus(AstS1,AstS2)
            % Description: Subtract AstSpec arrays (-).
            %              See astspec_arith for details.
            % Input   : - AstSpec array
            %           - AstSpec array, scalar or vector.
            % Outoput : - AstSpec class in which the subtraction of the
            %             .Int, .Back fields are given, an the
            %             error is propogated into .Err
            % Example : S-S
            
            AstS1 = astspec_arith(AstS1,AstS2,@minus);
            
        end
        
        function AstS1=rdivide(AstS1,AstS2)
            % Description: Divide AstSpec arrays (./).
            %              See astspec_arith for details.
            % Input   : - AstSpec array
            %           - AstSpec array, scalar or vector.
            % Outoput : - AstSpec class in which the division of the
            %             .Int, .Back fields are given, an the
            %             error is propogated into .Err
            % Example : A./A
            
            AstS1 = astspec_arith(AstS1,AstS2,@rdivide);
            
        end
        
        function AstS1=times(AstS1,AstS2)
            % Description: Multiply AstSpec arrays (.*).
            %              See astspec_arith for details.
            % Input   : - AstSpec array
            %           - AstSpec array, scalar or vector.
            % Outoput : - AstSpec class in which the multiplication of the
            %             .Int, .Back fields are given, an the
            %             error is propogated into .Err
            % Example : A.*A
            
            AstS1 = astspec_arith(AstS1,AstS2,@times);
            
        end
                
    end
        
end

