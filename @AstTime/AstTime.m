% AstTime class: Astronomical time class
% Package: class/@AstTime
% Description: A class to store and manipulate time formats and time
%              systems. The class is a structure array like object that
%              contains the following properties:
%                'Time' - An array of the time value.
%                'Type' - Time format.
%                'System' - Time system.
%              The class may contain multiple elements, and each element
%              may be an array of times.
%              The following time types are supported:
%                'jd' - Julian day (JD).
%                'mjd' - Modified JD (JD-2400000.5).
%                'rjd' - Reduced JD (JD-2400000).
%                'jy'|'jyear' - Julian year: 2000 + (JD-2451545)./365.25
%                'jcy' - Juilan century since J2000: (JD-2451545)./36525
%                'by'|'byear' - Bessilian year: 1900 + (JD-2415020.3135)/365.2421988
%                'iso'|'isot' - ISO time string 'YYYY-MM-DDTHH:MM:SS.FFF'
%                'dmyhms' - [D M Y H M S]
%                'ymdhms' - [Y M D H M S]
%                'dmyf'   - [D M Y Frac]
%                'ymdf'   - [Y M D Frac]
%                'unix'   - Unix time (seconds since UTC JD=2440587.5)
%              The following time systems are supported:
%                'UTC' - Coordinate Universal Time.
%                'UT1' - Universal time corrected for the polar motion.
%                'TAI' - International Atomic Time.
%                'TT'  - Terrestrial Time.
%                'TDB' - Dynamical Barycentric Time.
%              The time systems and Earth orientation parameters are
%              retrieved from the IERS website.
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef AstTime
    properties (SetAccess = public)
        Time     % 'now'|...         
        Type     % JDref|'jd'|'mjd'|'rjd'|'jyear'|'jcy'|'byear'|'iso'|'dmyhms'|'ymdhms'|'unix'|dmyf'|'ymdf'
        System   % 'UTC'|'UT1'|'TAI'|'TT'|'TDB'
    end
  
    % static methods for fields name and properties content
    methods (Static)
        function F=TimeField
           % AstTime class Time field name
           % Package: @AstTime
           F = 'Time';
        end
        function F=TypeField
           % AstTime class Type field name
           % Package: @AstTime
           F = 'Type';
        end
        function F=SystemField
           % AstTime class System field name
           % Package: @AstTime
           F = 'System';
        end
        function [Cell,Desc,Units]=timetype
            % Return all possible time types for AstTime object
            % Package: @AstTime
            % Input  : null
            % Output : - Cell array of time type options.
            %          - Cell array of time type description.
            %          - Cell array of units
            % Example: [Cell,Desc,Units]=AstTime.timetype
            Cell = {0,'jd','mjd','rjd','jy','jcy','by','iso','dmyhms','ymdhms','unix','dmyf','ymdf'};
            Desc = {'Ref Julian Day','Modified Julian Day','Truncated Julian Day',...
                    'Julian Year','Julian century','Besselian Year','ISO time string',...
                    '[D M Y H M S]','[Y M D H M S]','[D M Y Frac]','[Y M D Frac]'};
            Units = {'day','day','day','day','Julian year','Julian Century','Besselian year','ISO date','DMYHMS','YMDHMS','DMYF','YMDF'};
        end
        function [Cell,Desc]=timesystem
            % Return all possible time systems for AstTime object
            % Package: @AstTime
            % Input  : null
            % Output : - Cell array of time system options.
            %          - Cell array of time system description.
            % Example: [Cell,Desc]=AstTime.timesystem
            Cell = {'UTC','UT1','TAI','TT','TDB'};
            Desc = {'Coordinated Universal Time',...
                    'Universal Time 1 (corrected for polar motion)',...
                    'International Atomic Time',...
                    'Terrestrial Time',...
                    'Barycentric Dynamical Time'};
        end
    end
    
        
    % AstTime constructor
    methods
        
        function AstT=AstTime(Time,Type,System)
            % AstTime constructor method
            % Package: @AstTime
            % Description: AstTime constructor method
            % Input  : - Time. Default is now.
            %          - Time type (see AstTime.timetype)
            %            Default is 'JD'.
            %          - Time system (see AstTime.timesystem)
            %            Default is 'UTC'.
            
            Def.Time   = celestial.time.julday;
            Def.Type   = 'JD';
            Def.System = 'UTC';
            if (nargin==0)
                Time   = Def.Time;
                Type   = Def.Type;
                System = Def.System;
            elseif (nargin==1)
                Type   = Def.Type;
                System = Def.System;
            elseif (nargin==2)
                System = Def.System;
            else
                % do nothing
            end
            
            AstT.(AstTime.TimeField)   = Time;
            AstT.(AstTime.TypeField)   = Type;
            AstT.(AstTime.SystemField) = System;
             
        end

    end
    
    % set/get
    methods
        
        
    end
    
    % conversion
    methods
        function Out=convert_type2array(AstT,OutUnits)
            % convert type of an AstTime object and return an array
            % Package: @AstTime
            % Input  : - An AstTime object
            %          - Output time type (see AstTime.timetype)
            % Output : - Array of time with the requested type
           
            TimeField  = AstTime.TimeField;
            TypeField  = AstTime.TypeField;
            Nel = numel(AstT);
%             if (Nel>1)
%                 error('julday method work on a single element AstTime object');
%             end
            
            Value   = AstT.(TimeField);
            if (isnumeric(AstT.(TypeField)))
                % user supplied units in days
                RefJD  = AstT.(TypeField);
                Units  = 1;
            else
                % string format units
                % JDref|'jd'|'mjd'|'tjd'|'jyear'|'byear'|'iso'|'dmyhms'|'ymdhms'
                switch lower(AstT.(TypeField))
                    case 'jd'
                        RefJD = 0;
                        Units = 1;  % day
                        JD = Value.*Units + RefJD;
                    case 'mjd'
                        RefJD = 2400000.5;
                        Units = 1;  % day
                        JD = Value.*Units + RefJD;
                    case 'rjd'
                        % reduced JD
                        RefJD = 2400000;
                        Units = 1;
                        JD = Value.*Units + RefJD;
                    case 'unix'
                        % UNIX time: 1 Jan 1970
                        RefJD = 2440587.5;
                        Units = 1./86400;
                        JD = Value.*Units + RefJD;
                    case {'jy','jyear','julianyear'}
                        RefJD = 2451545.0;
                        Units = 365.25;
                        JD = (Value-2000).*Units + RefJD;
                    case 'jcy'
                        RefJD = 2451545.0;
                        Units = 36525;
                        JD = Value.*Units + RefJD;
                    case {'by','byear','besselianyear'}
                        RefJD = 2415020.3135;
                        Units = 365.2421988;
                        JD = (Value-1900).*Units + RefJD;
                    case {'iso','isodate','isot'}
                        % ISO date string
                        % special treatment
                        Value  = celestial.time.julday(AstT.(TimeField));
                        % Time already in JD
                        RefJD  = 0;
                        Units = 1;
                        JD = Value.*Units + RefJD;
                    case 'dmyhms'
                        % 3-6 column matrix of [D M Y H M S]
                        Value = celestial.time.julday(AstT.(TimeField));
                        % Time already in JD
                        RefJD = 0;
                        Units = 1;
                        JD = Value.*Units + RefJD;
                    case 'ymdhms'
                        % 3-6 column matrix of [Y M D H M S]
                        Value = celestial.time.julday(AstT.(TimeField)(:,[3 2 1 4 5 6]));
                        % Time already in JD
                        RefJD = 0;
                        Units = 1;
                        JD = Value.*Units + RefJD;
                    case 'dmyf'
                        % 4 column matrix [D M Y Frac]
                        Value = convert.julday(AstT.(TimeField));
                        RefJD = 0;
                        Units = 1;
                        JD = Value.*Units + RefJD;
                    case 'ymdf'
                        % 4 column matrix [Y M D Frac]
                        Value = celestial.time.julday(AstT.(TimeField)(:,[3 2 1 4]));
                        RefJD = 0;
                        Units = 1;
                        JD = Value.*Units + RefJD;
                    otherwise
                        error('Unknown AstTime class time type');
                end
            end
            

            % convert JD to output type:
            if (isnumeric(AstT.(TypeField)))
                % user supplied units in days
                RefJD  = AstT.(TypeField);
                Units  = 1;
            else
                % string format units
                % JDref|'jd'|'mjd'|'tjd'|'jyear'|'byear'|'iso'|'dmyhms'|'ymdhms'
                switch lower(OutUnits)
                    case 'jd'
                        RefJD = 0;
                        Units = 1;  % day
                        Out   = JD;
                    case 'mjd'
                        RefJD = 2400000.5;
                        Units = 1;  % day
                        Out   = JD - RefJD;
                    case 'rjd'
                        % reduced JD
                        RefJD = 2400000;
                        Units = 1;
                        Out   = JD - RefJD;
                    case 'unix'
                        % UNIX time: 1 Jan 1970
                        RefJD = 2440587.5;
                        Units = 86400;
                        Out   = (JD - RefJD).*Units;
                    case {'jy','jyear','julianyear'}
                        RefJD = 2451545.0;
                        Units = 365.25;
                        Out   = 2000 + (JD - RefJD)./Units;
                    case 'jcy'
                        RefJD = 2451545.0;
                        Units = 36525;
                        Out   = (JD - RefJD)./Units;
                    case {'by','byear','besselianyear'}
                        RefJD = 2415020.3135;
                        Units = 365.2421988;
                        Out   = 1900 + (JD - RefJD)./Units;
                    case {'iso','isodate','isot'}
                        % ISO date string
                        % special treatment
                        Out = convert.date2str(celestial.time.jd2date(JD(:),'H'));
                    case 'dmyhms'
                        % 3-6 column matrix of [D M Y H M S]
                        Out = celestial.time.jd2date(JD(:),'H');
                    case 'ymdhms'
                        % 3-6 column matrix of [Y M D H M S]
                        Out = celestial.time.jd2date(JD(:),'H');
                        Out = Out(:,[3 2 1 4 5 6]);
                    case 'dmyf'
                        % [D M Y Frac]
                        Out = celestial.time.jd2date(JD(:),'f');
                    case 'ymdf'
                        % [Y M D Frac]
                        Out = celestial.time.jd2date(JD(:),'f');
                        Out = Out(:,[3 2 1 4]);
                    otherwise
                        error('Unknown AstTime class time type');
                end
                        
            end
        end
        
        function AstT=convert_type(AstT,OutUnits)
            % convert type of an AstTime object and return an AstTime object
            % Package: @AstTime
            % Input  : - An AstTime object
            %          - Output time type (see AstTime.timetype)
            % Output : - An AstTime object of the requested type
            % Example: T=AstTime(2451545);
            %          MJD=convert_type(T,'mjd')
           
            TimeField = AstTime.TimeField;
            TypeField = AstTime.TypeField;
            
            Nel = numel(AstT);
            for Iel=1:1:Nel
                AstT(Iel).(TimeField) = convert_type2array(AstT(Iel),OutUnits);
                AstT(Iel).(TypeField) = OutUnits;
            end
             
        end
        
        function Out=convert_sys2array(AstT,OutSys,Mode)
            % Convert AstTime time system and return an array
            % Package: @AstTime
            % Input  : - A single element AstTime object 
            %          - Output system (see AstTime.timesystem)
            %          - Mode 'use'|'get' for retrieve time systems tables
            %            Default is 'use'.
            % Output : - An array of the time in the requested system
            % Example: T=AstTime(2451545);
            %          TAI=convert_sys2array(T,'TAI');
            
            if (nargin<3)
                Mode = 'use';
            end
            
            DAY_SECOND = 1./86400;
            if (nargin<2)
                Mode = 'use';
            end
            SystemField = AstTime.SystemField;
            TimeField   = AstTime.TimeField;
            TypeField   = AstTime.TypeField;
            
            % convert AstTime type to JD:
            AstT_JD = convert_type(AstT,'JD');
            
            InSystem = AstT_JD.(SystemField);
            % convert the time system to UTC
            AstT_JDUTC = AstT_JD;
            % 'UTC'|'UT1'|'TAI'|'TT'|'TDB'
            switch lower(InSystem)
                case 'utc'
                    % do nothing - already in UTC
                case 'ut1'
                    % convert UT1 to UTC
                    % get UT1-UTC
                    [UT1mUTC] = celestial.time.ut1_utc(AstT_JDUTC.(TimeField),Mode);
                    AstT_JDUTC.(TimeField) = AstT_JDUTC.(TimeField) - UT1mUTC.*DAY_SECOND;
                case 'tai'
                    % convert TAI to UTC
                    % get TAI-UTC
                    [TAImUTC] = celestial.time.tai_utc(AstT_JDUTC.(TimeField),Mode);
                    AstT_JDUTC.(TimeField) = AstT_JDUTC.(TimeField) - TAImUTC.*DAY_SECOND;
                case {'tt','tdt'}
                    % convert TT to UTC
                    % get TT-UTC
                    [~,TTmUTC] = celestial.time.tai_utc(AstT_JDUTC.(TimeField),Mode);
                    AstT_JDUTC.(TimeField) = AstT_JDUTC.(TimeField) - TTmUTC.*DAY_SECOND;
                case {'tdb'}
                    % convert TDB to UTC
                    % get TT-UTC
                    [~,TTmUTC] = celestial.time.tai_utc(AstT_JDUTC.(TimeField),Mode);
                    % get TDB-TT (approximate)
                    TDBmTT     = celestial.time.tdb_tdt(AstT_JDUTC.(TimeField));
                    AstT_JDUTC.(TimeField) = AstT_JDUTC.(TimeField) - (TDBmTT + TTmUTC).*DAY_SECOND;
                otherwise
                    error('Unknown time system option');
            end
            
            % convert UTC to output system
            AstT_JDout = AstT_JDUTC;
            switch lower(OutSys)
                case 'utc'
                    % do nothing
                case 'ut1'
                    % convert UTC to UT1
                    % get UT1-UTC
                    [UT1mUTC] = celestial.time.ut1_utc(AstT_JDout.(TimeField),Mode);
                    AstT_JDout.(TimeField) = AstT_JDout.(TimeField) + UT1mUTC.*DAY_SECOND;
                case 'tai'
                    % convert UTC to TAI
                    % get TAI-UTC
                    [TAImUTC] = celestial.time.tai_utc(AstT_JDout.(TimeField),Mode);
                    AstT_JDout.(TimeField) = AstT_JDout.(TimeField) + TAImUTC.*DAY_SECOND;
                case {'tt','tdt'}
                    % convert UTC to TT 
                    % get TT-UTC
                    [~,TTmUTC] = celestial.time.tai_utc(AstT_JDout.(TimeField),Mode);
                    AstT_JDout.(TimeField) = AstT_JDout.(TimeField) + TTmUTC.*DAY_SECOND;
                case {'tdb'}
                    % convert UTC to TDB
                    % get TT-UTC
                    [~,TTmUTC] = celestial.time.tai_utc(AstT_JDout.(TimeField),Mode);
                    % get TDB-TT (approximate)
                    TDBmTT     = celestial.time.tdb_tdt(AstT_JDout.(TimeField));
                    AstT_JDout.(TimeField) = AstT_JDout.(TimeField) + (TDBmTT + TTmUTC).*DAY_SECOND;    
                otherwise
                    error('Unknown time system option');
            end    
            
            AstT_JDout.(SystemField) = OutSys;
            % convert back to original time type
            AstT = convert_type(AstT_JDout,AstT.(TypeField));
            Out  = AstT.(TimeField);
            
        end
        
        function AstT=convert_sys(AstT,OutSys,Mode)
            % Convert AstTime time system and return an AstTime object
            % Package: @AstTime
            % Input  : - An AstTime object 
            %          - Output system (see AstTime.timesystem)
            %          - Mode 'use'|'get' for retrieve time systems tables
            %            Default is 'use'.
            % Output : - An AstTime object of the time in the requested system
            % Example: T=AstTime(2451545);
            %          TAI=convert_sys(T,'TAI');
            
            if (nargin<3)
                Mode = 'use';
            end
           
            TimeField   = AstTime.TimeField;
            SystemField = AstTime.SystemField; 
            Nel = numel(AstT);
            for Iel=1:1:Nel
                AstT(Iel).(TimeField)   = convert_sys2array(AstT(Iel),OutSys,Mode);
                AstT(Iel).(SystemField) = OutSys;
            end
            
        end
        
        function AstT=convert(AstT,OutType,OutSys,Mode)
            % convert AstTime object time type and time system
            % Package: @AstTime
            % Input  : - An AstTime object
            %          - Time type to convert into. See AstTime.timetype
            %          - Time system to convert into.
            %            See AstTime.timesystem.
            %            Default is 'UTC'.
            %          - Mode for getting time system information:
            %            'use'|'get'. Default is 'use'.
            % Output : - An AstTime object with the requested time type and
            %            system.
            % Example: T=AstTime(2451545);
            %          AstT=convert(T,'iso','UT1');
            
            if (nargin<4)
                Mode = 'use';
                if (nargin<3)
                    OutSys = 'UTC';
                end
            end
            
            AstT = convert_type(AstT,OutType);
            AstT = convert_sys(AstT,OutSys,Mode);
            
            
        end
        
        function AstT=convert2array(AstT,OutType,OutSys,Mode)
            % convert AstTime object time type and time system to array
            % Package: @AstTime
            % Input  : - A single element AstTime object
            %          - Time type to convert into. See AstTime.timetype
            %          - Time system to convert into.
            %            See AstTime.timesystem.
            %            Default is 'UTC'.
            %          - Mode for getting time system information:
            %            'use'|'get'. Default is 'use'.
            % Output : - An AstTime object with the requested time type and
            %            system.
            % Example: T=AstTime(2451545);
            %          AstT=convert2array(T,'iso','UT1');
            
            if (numel(AstT)>1)
                error('AstTime object must contain a single element');
            end
            
            if (nargin<4)
                Mode = 'use';
                if (nargin<3)
                    OutSys = 'UTC';
                end
            end
            
            AstT = convert_type(AstT,OutType);
            AstT = convert_sys2array(AstT,OutSys,Mode);
            
            
        end
        
    end
    
    % operators implementation
    methods
        function AstT=bfun2asttime(AstT1,AstT2,Operator)
            % Binary operator on AstTime object
            % Package: @AstTime
            % Input  : - An AstTime object
            %          - Either a scalar, array of the size of the first input
            %            AstTime object, or another AstTime object.
            %            If a scalar or vector than this is a time in units
            %            of days (i.e., 86400 SI seconds).
            %          - Binary operator (e.g., @plus).
            % Output : - An AstTime object
            % Example: T=AstTime(2451545); bfun2asttime(T,T,@plus)
            %          T=convert(T,'jy'); bfun2asttime(T,365,@plus)
            TimeField   = AstTime.TimeField;
            TypeField   = AstTime.TypeField;
            SystemField = AstTime.SystemField;
         
            AstT = convert(AstT1,'jd','utc');
            N1   = numel(AstT1);
            N2   = numel(AstT2);
            if ~(N2==1 || N2==N1)
                error('Size of AstT2 must be 1 or as the size of AstT1');
            end
            if (isnumeric(AstT2))
                % AstT2 is an array
                % each element in AstT2 is pairs with element in AstT1
                for I1=1:1:N1
                    I2 = min(N2,I1);
                    % Assume AstT2 is in units of days
                    AstT(I1).Time = Operator(AstT(I1).(TimeField),AstT2(I2));
                end
            else
                % assume AstT2 is an AstTime object
                AstT2 = convert(AstT2,'jd','utc');
                for I1=1:1:N1
                    I2 = min(N2,I1);
                    AstT(I1).(TimeField) = Operator(AstT(I1).(TimeField), AstT2(I2).(TimeField));
                end
            end
            % return AstT to original units
            AstT = convert(AstT,AstT1(1).(TypeField),AstT1(1).(SystemField));
                    
                
            
        end
        
        function AstTV=ufun2asttime(AstT,Operator,Pars)
            % Apply a unary operator to an AstTime and return AstTime
            % Package: @AstTime
            % Input  : - An AstTime object
            %          - Operator function handle.
            %          - Cell array of additional parameters to pass to the 
            %            operator function. Default is {}.
            % Output : - An AstTime object
            % Example: T = AstTime(2451545+(0:1:10)');
            %          AstTV=ufun2asttime(T,@mean)
            %          AstTV=ufun2asttime(T,@sin)
           
            if (nargin<3)
                Pars = {};
            end
            TimeField  = AstTime.TimeField;
            TypeField  = AstTime.TypeField;
            AstTV      = convert_type(AstT,'JD');
            
            Nel = numel(AstT);
            for Iel=1:1:Nel
                AstTV(Iel).(TimeField) = Operator(AstTV(Iel).(TimeField),Pars{:});
                AstTV(Iel) = convert_type(AstTV(Iel),AstT(Iel).(TypeField));
            end
            
        end
        
        function Array=ufun2scalar(AstT,Operator,varargin)
            % Apply an operator that return a scalar to each elements in AstTime
            % Package: @AstTime
            % Description: Apply an operator on each element in AstTime and
            %              return the operator scalar result in an array 
            %              (element in the array per elementin the AstTime object).
            %              The times in each AstTime element are treated as
            %              a column vector (e.g., the mean will be
            %              calculated for all elements).
            % Input  : - An AstTime object
            %          - Opertaor function handle.
            %          * Additional arguments to pass to the operator
            % Output : - An array in which each element is the result of
            %            the operator on the corresponding AstTime object.
            %            Note that regardless the input type the result is
            %            always represented in days (864100 SI seconds).
            % Example: T(1) = AstTime(2451545+(1:1:10)');
            %          T(2) = AstTime((1:1:10)','jy');
            %          A = ufun2scalar(T,@std);
            
            Nel   = numel(AstT);
            Array = zeros(size(AstT));
            for Iel=1:1:Nel
                Tmp = convert_type2array(AstT(Iel),'JD');
                Array(Iel) = Operator(Tmp(:),varargin{:});
            end
           
        end
    end
    
    % operators / overload
    methods
        
        % plus
        function AstT=plus(AstT1,AstT2)
            % Plus operator on AstTime objects
            % Package: @AstTime
            % Input  : - An AstTime object
            %        : - A scalar, array, or an AstTime object.
            %            Scalar will be added to all elements, each array 
            %            element is added to each AstTime element.
            %            If a scalar or vector than this is a time in units
            %            of days (i.e., 86400 SI seconds).
            % Output : - An AstTime object contain the addition of the two
            %            input objects, the format and time system of the
            %            first AstTime objects are assumed.
            AstT = bfun2asttime(AstT1,AstT2,@plus);
            
        end
        
        % minus
        function AstT=minus(AstT1,AstT2)
            % Minus operator on AstTime objects
            % Package: @AstTime
            % Input  : - An AstTime object
            %        : - A scalar, array, or an AstTime object.
            %            Scalar will be subtracted from all elements, each array 
            %            element is subtracted from each AstTime element.
            %            If a scalar or vector than this is a time in units
            %            of days (i.e., 86400 SI seconds).
            % Output : - An AstTime object contain the subtraction of the two
            %            input objects, the format and time system of the
            %            first AstTime objects are assumed.
            AstT = bfun2asttime(AstT1,AstT2,@minus);
            
        end
        
        % mean time
        function MeanAstT=mean(AstT,varargin)
            % Calculate the mean for each element in an AstTime object
            % Package: @AstTime
            % Input  : - An AstTime object
            %          * Additional arguments to pass to the mean function
            % Output : - An AstTime object with the mean time in each
            %            element.
            % Example: T=AstTime(2451545+(1:1:10)');
            %          convert(mean(T),'iso')
            
            MeanAstT = ufun2asttime(AstT,@mean,varargin{:});
            
        end
        
        % median time
        function MedAstT=median(AstT,varargin)
            % Calculate the median for each element in an AstTime object
            % Package: @AstTime
            % Input  : - An AstTime object
            %          * Additional arguments to pass to the median function
            % Output : - An AstTime object with the median time in each
            %            element.
            % Example: T=AstTime(2451545+(1:1:10)');
            %          convert(median(T),'iso')
            
            MedAstT = ufun2asttime(AstT,@median,varargin{:});
            
        end
        
        % mean time to array
        function Array=mean2array(AstT)
            % Array of the mean of each AstTime object element
            % Package: @AstTime
            % Input  : - An AstTime object
            % Output : - An array with the mean of each element of an
            %            AstTime object.
            % Example: T=AstTime(2451545+(1:1:10)');
            %          mean2array(T)
            
            Array = ufun2scalar(AstT,@mean);
        end
        
        % median time to array
        function Array=median2array(AstT)
            % Array of the median of each AstTime object element
            % Package: @AstTime
            % Input  : - An AstTime object
            % Output : - An array with the median of each element of an
            %            AstTime object.
            % Example: T=AstTime(2451545+(1:1:10)');
            %          median2array(T)
            
            Array = ufun2scalar(AstT,@median);
        end
        
        % std time to array
        function Array=std2array(AstT,varargin)
            % Array of the std of each AstTime object element
            % Package: @AstTime
            % Input  : - An AstTime object
            %          * Additional arguments to pass to the std function.
            % Output : - An array with the std of each element of an
            %            AstTime object.
            % Example: T=AstTime(2451545+(1:1:10)');
            %          std2array(T)
            
            Array = ufun2scalar(AstT,@std,varargin{:});
        end
        
        % std/sqrt(N) time to array
        function Array=errormean2array(AstT,varargin)
            % Array of the error on the mean of each AstTime object element
            % Package: @AstTime
            % Input  : - An AstTime object
            %          * Additional arguments to pass to the std function.
            % Output : - An array with the error on the mean (std/sqrt(N))
            %            of each element of an AstTime object.
            % Example: T=AstTime(2451545+(1:1:10)');
            %          errormean2array(T)
            
            Array = ufun2scalar(AstT,@mean_error,varargin{:});
        end
        
        % range time to array
        function Array=range2array(AstT)
            % Array of the range of each AstTime object element
            % Package: @AstTime
            % Input  : - An AstTime object
            % Output : - An array with the range of each element of an
            %            AstTime object.
            % Example: T=AstTime(2451545+(1:1:10)');
            %          range2array(T)
            
            Array = ufun2scalar(AstT,@range);
        end
        
        % histogram of time
        function [N,Edges]=histcounts(AstT,varargin)
            % Calculate the histogram of time in an AstTime object
            % Package: @AstTime
            % Input  : - A single element AstTime object.
            %            The times will be converted to JD.
            %          * Additional aeguments to pass to the histcounts
            %            function. Default is {}.
            % Output : - Number of events in bin.
            %          - Bins edges.
            % Example: T=AstTime(2451545+(1:1:100));
            %          [N,E]=histcounts(T);
            
            Array = julday(AstT);
            [N,Edges] = histcounts(Array(:),varargin{:});
        end
    end
    
    % spatial information
    methods
        % lst
        function [LST,JD_UT1]=lst(AstT,EastLong,STtype)
            % Calculate the local sidereal time
            % Package: @AstTime
            % Input  : - A single element AstTime object
            %          - East longitude [radians].
            %            Default is 0.
            %          - Sidereal time type:
            %            'm' - Mean. Default.
            %            'a' - Apparent.
            % Output : - An array of LST in fraction of day.
            %          - An array of JD in the UT1 time system.
            % Example: T=AstTime(2451545); [LST,JD_UT1]=lst(T)
           
            if (numel(AstT)>1)
                error('AstTime object must contain a single element');
            end
            
            if (nargin<3)
                STtype = 'm';
                if (nargin<2)
                    EastLong = 0;
                end
            end
            JD_UT1  = convert2array(AstT,'JD','UT1');
            
            LST = celestial.time.lst(JD_UT1,EastLong,STtype);
        end
        
        % month name
        % JD
        function JD=julday(AstT,TimeSystem)
            % Convert an AstTime object to an array of JD
            % Package: @AstTime
            % Input  : - A single element AstTime object
            %          - Time system. Default is the AstTime object time
            %            system. See AStTime.timesystem for options.
            % Output : - An array of JD.
            % Example: julday(AstTime)
        
            if (numel(AstT)>1)
                error('AstTime should contain a single element');
            end
            SystemField = AstTime.SystemField;
            
            if (nargin<2)
                TimeSystem = AstT.(SystemField);
            end
            JD = convert2array(AstT,'JD',TimeSystem);
        end
        
        % delta_t
        % TAI-
        % Sidereal month
        
        % Synodic month
        function SM=synodic_month(AstT)
            % The mean synodic month at the time in an AstTime object
            % Package: @AstTime
            % Input  : - A single element AstTime object
            % Output : - An array of the mean synodic month [days]
            % Example: synodic_month(AstTime)
           
            T  = convert2array(AstT,'jcy','UTC');
            SM = 29.5305888531 + 0.00000021621.*T - 3.64.*1e-10.*T.^2;
            
        end
        % Anomalistic month
        % Draconic month
        % Tropical year
        % Sidereal year
        % anomalistic year
        % Draconic year
        % Hebrew year
        
        % obliquity
        function Obl=obliquity(AstT,varargin)
            % The mean obliquity of the Earth at an epochs given in AstTime object
            % Package: @AstTime
            % Input  : - A single element AstTime object
            %          - Calculation type:
            %            'L' - IAU 1976, good from 1000-3000 AD. Default.
            %            'H' - Laskar expression, more accurate.
            % Output : - An array of obliquity (radians).
            % Example: obliquity(AstTime)
           
            if (numel(AstT)>1)
                error('AstTime should contain a single element');
            end
            JD = convert2array(AstT,'JD','UTC');
            Obl = celestial.coo.obliquity(JD,varargin{:});
        end
 
        % Nutation
        % Precession
        
    end
    
    % spatial information (static)
    methods (Static)
        % MeanDay
        function L=day
            % Return the length of the day = 86400 SI seconds
            % Package: @AstTime
            % Input  : null
            % Output : 86400 s.
            L = 86400;
        end
            
        % Mean sidereal day at J2000.0
        function L=sidereal_day
            % Return the length of the sidereal day at J2000
            % Package: @AstTime
            % Input  : null
            % Output : 86164.091 s.
            L = 86164.091;
        end
        
        % Mean sidearl day short term
        
        % Mean sidereal day
%         function LOD=sidereal_day_long(T)
%             % Length of the sideral day at the +/-0.25 billion years.
%             % Input  : Vector of time in years since J2000.0
%             % Output : Length of sidereal day in [days=86400 SI s]
%             % Reference: Laskar et al. 2004
%             % Example: LOD=AstTime.sidereal_day_long(-1e8)
%             T = T./1e9;
%             if (any(abs(T)>0.25))
%                 error('Time in years must be between -0.25 to 0.25 billion years');
%             end
%             
%             LOD = T;
%             Ipast   = T<0;
%             Ifuture = T>0;
%             LOD(Ipast)   = 23.934468 + 7.432167.*T(Ipast)   - 0.727046.*T(Ipast).^2   + 0.409572.*T(Ipast).^3 - 0.589692.*T(Ipast).^4;
%             LOD(Ifuture) = 23.934468 + 7.444649.*T(Ifuture) - 0.715049.*T(Ifuture).^2 + 0.458097.*T(Ifuture).^3;
%             LOD = LOD./24;
%             
%         end
        
        % ...
        
    end
end

            
