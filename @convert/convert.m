%--------------------------------------------------------------------------
% convert class                                                      class
% Description: A static class for conversion functions.
%              This class include many static function for units
%              conversion.
%              Highlights include:
%              convert.units - General units conversion (e.g., 'cm/s' to
%                              km/h').
%              convert.angular - Convert between angular units.
%              convert.energy  - Convert between energy units.
%              convert.flux    - Convert fluxes (e.g., AB mag to ph/s/Ang).
%              convert.time    - Years, to JD etc.
%              convert.date2jd - Date to JD.
%              convert.jd2date - JD to date.
%              Type "convert." followed by <tab> to see the full list of
%              functions.
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Aug 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef convert 
             
    
    methods (Static)
    
        % General units convertor
        function [Factor,Out,String]=units(InUnits,OutUnits,In)
            % Unit conversion
            % Package: @convert
            % Description: Unit conversion function. Given an input and output strings
            %              containing unit names, return the conversion multiplication 
            %              factor needed for converting the input units to the output
            %              units. The user is responsible for the balance of the
            %              transformation.
            %              Type of units:
            %               Length:
            %                  'ang','mm' ; 'cm' ; 'inch' ; 'feet' ;
            %                  'micron', 'm' - meter ; 'km' ; 'mile' ;
            %                  'earth_rad' - Earth radius ; 'au' ;
            %                  'ly' - light year ; 'pc'; 'yard'
            %                Time:
            %                  's' ; 'min' ; 'hour' ; 'day';
            %                  'sday' - sidereal day ; week ;
            %                  'year'; 'cen' - century
            %                Mass:
            %                  'gr'/'g'; 'kg'; 'earth_mass' - Earth mass;
            %                  'jupiter_mass' - Jupiter mass;
            %                  'sun_mass' - Solar mass;
            %                  'mp' - proton mass;
            %                  'me' - electron mass;
            %                  'libra';'pound'
            %                Energy: (see also convert_energy.m)
            %                  'erg'; 'J'
            %                Angle:
            %                  'rad' ; 'deg' ;
            %                  'amin' (or 'arcmin') - arcmin ; 'asec' (or 'arcsec') -
            %                  arcsec, 'mas'
            %                Solid Angle:
            %                  'ster' ; 'sdeg' - square degree ;
            %                  'smin' - square arcmin ;
            %                  'ssec' - square arcsec
            % Input  : - String containing the input units.
            %          - String containing the output units.
            %          - Optional value to convert from input units to output
            %            units. Default is 1.
            % Output : - Multiplication factor for converting input units to
            %            output units.
            %          - The input value given in the output units.
            %          - A string describing the transformation (e.g., 
            % Example : convert.units('m^3 * kg^-1 * s^-2','cm^3 * gr^-1 * s^-2')
            % Tested : Matlab 6.5
            %     By : Eran O. Ofek                    Jul 2003
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [Factor,Out,String]=convert.units('erg','J',1)
            % Reliable: 1
            %--------------------------------------------------------------------------
            RAD  = 180./pi;

            if (nargin==2)
                In = 1;
            end

            InUnits  = lower(InUnits);
            OutUnits = lower(OutUnits);

            %--- Length units ---
            ang  = 1e-10;                        % ang
            mm   = 1e-3;                         % mm
            cm   = 1e-2;                         % cm
            micron = 1e-6;                       % micron
            inch = 0.0254;                       % inch
            feet = 0.30480;                      % foot
            yard = 0.9144;                       % yard
            m    = 1;                            % meter
            km   = 1000;                         % km
            mile = 1609;                         % mile
            earth_rad = celestial.Earth.refellipsoid('WGS84');
            earth_rad = earth_rad(1);                      % Earth equatorial radius (WGS84)
            au   = constant.au('SI');      % au
            ly   = constant.ly('SI');      % light-year
            pc   = constant.pc('SI');      % pc


            %--- Time units ---
            s    = 1;                            % sec
            min  = 60;                           % minute
            hour = 3600;                         % hour
            sday = 86164.09053;                  % sidereal day
            day  = 86400;                        % day
            week = 7.*day;                       % week
            year = 365.25.*day;                  % tropical year
            yr   = year;
            cen  = year*100;                     % tropical century

            %--- Mass units ---
            gr   = 1e-3;                         % gram
            g    = 1e-3;                         % gram
            kg   = 1;                            % kg
            earth_mass= constant.EarthM('SI');    % Earth Mass
            jupiter_mass= constant.JupiterM('SI');    % Jupiter Mass
            sun_mass= constant.SunM('SI');    % Solar Mass
            me   = constant.me('SI');      % electron mass
            mp   = constant.mp('SI');      % proton mass
            libra= 0.32736409;                   % libra
            pound= 0.45359237;                   % pound

            %--- Energy Units ---
            erg  = 1e-7;                         % erg
            j    = 1;                            % joul

            %--- Angle units ---
            rad  = 1;                            % radian
            deg  = 1./RAD;                       % degree
            amin = deg./60;                      % arcmin
            asec = amin./60;                     % arcsec
            arcmin = deg./60;
            arcsec = arcmin./60;
            mas    = arcsec/1000;

            %--- Solid Angle units ---
            ster = 1;                            % steradian
            sdeg = 1./(RAD.^2);                  % square degree
            smin = 1./((60.*RAD).^2);            % square arcmin
            ssec = 1./((3600.*RAD).^2);          % square arcsec


            %--- Find conversion factor ---
            F1   = eval(InUnits);
            F2   = eval(OutUnits);

            Factor = F1./F2;

            if (nargout>1)
               Out = Factor.*In;
            end

            if (nargout>2)
               String = sprintf('1 %s = %g %s',InUnits,Factor,OutUnits);
            end

        end % convert.units function

    end % Static


    % conversion functions
    % (angular, energy, flux)
    methods (Static)
        % angular conversion
        function ConvVal=angular(In,Out,Val,SpecialString)
            % Convert angular units
            % Package: @convert
            % Description: Convert angular units.
            % Input  : - Input angular units name:
            %            'rad'|'deg'|'arcmin'|'arcsec'|'mas'|'frac'|'hour'
            %            Note that the 'pix' option generate a conversion factor of 1.
            %          - Output angular units (options like first input argument).
            %          - Value to convert. Default is 1.
            %          - A string or a cell array of string that if the first input
            %            argument (In) is equal to one of them then the function will
            %            set the conversion factor to 1.
            % Output : - Converted value (by default this is the conversion factor).
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Jun 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: convert.angular('rad','mas');
            %          convert.angular('rad','mas',1,'pix');
            % Reliable: 2
            %--------------------------------------------------------------------------

            arguments
                In
                Out
                Val             = 1;
                SpecialString   = [];
            end
            
            if strcmp(In,Out)
                % input=output
                ConvVal = Val;
            else
                
                if any(strcmpi(In,SpecialString))
                    % user used Special string - set output to 1
                    ConvVal = Val;
                else

                    % Input: convert to deg
                    switch lower(In)
                        case {'rad','radian','radians'}
                            ConvVal = 180./pi;
                        case {'deg','degree','degrees'}
                            ConvVal = 1;
                        case 'arcmin'
                            ConvVal = 1./60;
                        case 'arcsec'
                            ConvVal = 1./3600;
                        case 'mas'
                            ConvVal = 1./(3600.*1000);
                        case 'frac'
                            % fraction to deg
                            ConvVal = 360;
                        case 'hour'
                            % hours to deg
                            ConvVal = 15;
                        otherwise
                            error('Unknown input angular units');
                    end
                    ConvVal = Val.*ConvVal;

                    % Output: Convert from deg to requested output units
                    switch lower(Out)
                        case {'rad','radian','radians'}
                            ConvVal = ConvVal.* pi./180;
                        case {'deg','degree','degrees'}
                            % do nothing
                            %ConvVal = ConvVal;
                        case 'arcmin'
                            ConvVal = ConvVal.*60;
                        case 'arcsec'
                            ConvVal = ConvVal.*3600;
                        case 'mas'
                            ConvVal = ConvVal.*3600.*1000;
                        case 'frac'
                            ConvVal = ConvVal./360;
                        case 'hour'
                            ConvVal = ConvVal./15;
                        otherwise
                            error('Unknown input angular units');
                    end
                end
            end

        end % convert.angular function
        
        function Val = timeUnits(In, Out, Val)
            % convert time units
            % Input  : - In time units
            %            's'|'min'|'hr'|'day'|'sday'|'week'|'jyr'|'jcy'
            %          - Out Units (like input).
            %          - Value to convert. Default is 1.
            % Output : - Converted time
            % Author : Eran Ofek (May 2021)
            % Example: convert.timeUnits('jcy','min',1)
           
            arguments
                In
                Out
                Val    = 1;
            end
            
            % convert to [s]
            switch lower(In)
                case 's'
                    % do nothing
                case {'min'}
                    Val = Val.*60;
                case {'hr','hour'}
                    Val = Val.*3600;
                case {'day'}
                    Val = Val.*86400;
                case {'sday'}
                    Val = Val.*86164.091;
                case {'week'}
                    Val = Val.*86400.*7;
                case {'jyear','yr'}
                    Val = Val.*86400.*365.25;
                case {'jcy'}
                    Val = Val.*86400.*365.25.*100;
                otherwise
                    error('Unknown time units');
            end
            
            % convert from 's' to output
            switch lower(Out)
                case 's'
                    % do nothing
                case {'min'}
                    Val = Val./60;
                case {'hr','hour'}
                    Val = Val./3600;
                case {'day'}
                    Val = Val./86400;
                case {'sday'}
                    Val = Val./86164.091;
                case {'week'}
                    Val = Val./(86400.*7);
                case {'jyear','yr'}
                    Val = Val./(86400.*365.25);
                case {'jcy'}
                    Val = Val./(86400.*365.25.*100);
                otherwise
                    error('Unknown time units');
            end
                    
            
        end
        
        function Val = length(In, Out, Val)
            % Convert length units
            % Input  : - Input units:
            %            'cm'|'m'|'mm'|'ang'|...
            %          - Output units (like input).
            %          - Input value. Default is 1.
            % Output : - Converted output value.
            % AUthor : Eran Ofek (May 2021)
            % Example: convert.length('yard','m')
           
            arguments
                In
                Out
                Val    = 1;
            end
            
            % convert to 'cm'
            switch lower(In)
                case 'cm'
                    % do nothing
                case 'm'
                    Val = Val.*100;
                case 'mm'
                    Val = Val./10;
                case 'micrometer'
                    Val = Val.*1e-4;
                case {'a','ang'}
                    Val = Val.*1e-8;
                case 'km'
                    Val = Val.*1e5;
                case 'au'
                    Val = Val.*constant.au;
                case 'pc'
                    Val = Val.*constant.pc;
                case 'earthr'
                    Val = Val.*constant.EarthR;
                case 'sunr'
                    Val = Val.*constant.SunR;
                case 'inch'
                    Val = Val.*2.54;
                case 'mile'
                    Val = Val.*1609.344.*100;
                case 'nmile'
                    Val = Val.*1852.*100;
                case 'yard'
                    Val = Val.*91.44;
                case 'feet'
                    Val = Val.*30.48;
                otherwise
                    error('Unknown length units');
            end
            
              % convert to 'cm'
            switch lower(Out)
                case 'cm'
                    % do nothing
                case 'm'
                    Val = Val./100;
                case 'mm'
                    Val = Val.*10;
                case 'micrometer'
                    Val = Val.*1e4;
                case {'a','ang'}
                    Val = Val.*1e8;
                case 'km'
                    Val = Val.*1e-5;
                case 'au'
                    Val = Val./constant.au;
                case 'pc'
                    Val = Val./constant.pc;
                case 'earthr'
                    Val = Val./constant.EarthR;
                case 'sunr'
                    Val = Val./constant.SunR;
                case 'inch'
                    Val = Val./2.54;
                case 'mile'
                    Val = Val./(1609.344.*100);
                case 'nmile'
                    Val = Val./(1852.*100);
                case 'yard'
                    Val = Val./91.44;
                case 'feet'
                    Val = Val./30.48;
                otherwise
                    error('Unknown length units');
            end
            
                        
                    
            
            
        end
        
        function Val = velocity(In, Out, Val)
            % Convert velocity units of the form 'length/time'
            % Input  : - Input velocity units.
            %          - Output velocity units.
            %          - Value to convert. Default is 1.
            % Output : - Converted value.
            % Author : Eran Ofek (May 2021)
            % Example: convert.velocity('km/s','au/day',1)
           
            arguments
                In
                Out
                Val    = 1;
            end
            
            % identify length and time units
            InSplit = regexp(In,'/','split');
            if numel(InSplit)~=2
                error('In units does not look like velocity');
            end
            OutSplit = regexp(Out,'/','split');
            if numel(OutSplit)~=2
                error('Out units does not look like velocity');
            end
            LengthFactor = convert.length(InSplit{1}, OutSplit{1}, 1);
            TimeFactor   = convert.timeUnits(InSplit{2}, OutSplit{2}, 1);
            
            Val = Val.*LengthFactor./TimeFactor;
            
        end
        
        function Val = proper_motion(In, Out, Val)
            % Convert proper motion units of the form 'angle/time'
            % Input  : - Input velocity units.
            %          - Output velocity units.
            %          - Value to convert. Default is 1.
            % Output : - Converted value.
            % Author : Eran Ofek (May 2021)
            % Example: convert.proper_motion('arcsec/yr','mas/yr',1)
           
            arguments
                In
                Out
                Val    = 1;
            end
            
            % identify length and time units
            InSplit = regexp(In,'/','split');
            if numel(InSplit)~=2
                error('In units does not look like proper motion');
            end
            OutSplit = regexp(Out,'/','split');
            if numel(OutSplit)~=2
                error('Out units does not look like proper motion');
            end
            AngleFactor = convert.angular(InSplit{1}, OutSplit{1}, 1);
            TimeFactor   = convert.timeUnits(InSplit{2}, OutSplit{2}, 1);
            
            Val = Val.*AngleFactor./TimeFactor;
            
        end
        
        
        
        function Out=minusPi2Pi(In,Units)
            % convert angular units to -pi to pi range
            % Input  : - Angular values.
            %        : - Output and input units: 'rad' | 'deg'. Default is
            %            'deg'
            % Output : - Angular values in the range of -pi to pi.
           
            if nargin<2
                Units = 'deg';
            end
            
            switch lower(Units)
                case 'deg'
                    Out = mod(In,360);
                    Out(Out>180) = Out(Out>180) - 360;
                case 'rad'
                    Out = mod(In,2.*pi);
                    Out(Out>pi) = Out(Out>pi) - 2.*pi;
                otherwise
                    error('Unknown Units option');
            end
            
        end
        
        % mass conversion
        function OutM=mass(InUnit,OutUnit,InM)
            % Convert between different mass units
            % Package: @convert
            % Description: Convert between different mass units
            % Input  : - Input system:
            %            'gr'
            %            'kg'
            %            'SunM'
            %            'EarthM'
            %            'mp'
            %            'me'
            %          - Output system.
            %          - Input Mass. Default is 1.
            % Output : - Output mass.
            
            Def.InM = 1;
            if (nargin<3)
                InM = Def.InM;
            end
            
            % convert to gr:
            switch lower(InUnit)
                case 'gr'
                    Conv1 = 1;
                case 'kg'
                    Conv1 = 1000;
                case 'sunm'
                    Conv1 = constant.SunM;
                case 'earthm'
                    Conv1 = constant.EarthM;
                case 'mp'
                    Conv1 = constant.mp;
                case 'me'
                    Conv1 = constant.me;
                otherwise
                    error('Unknown InUnit option');
            end
            switch lower(OutUnit)
                case 'gr'
                    Conv2 = 1;
                case 'kg'
                    Conv2 = 1./1000;
                case 'sunm'
                    Conv2 = 1./constant.SunM;
                case 'earthm'
                    Conv2 = 1./constant.EarthM;
                case 'mp'
                    Conv2 = 1./constant.mp;
                case 'me'
                    Conv2 = 1./constant.me;
                otherwise
                    error('Unknown InUnit option');
            end
            
            OutM = InM.*Conv1.*Conv2;
            
        end
        
        % energy conversion
        function OutE=energy(InUnit,OutUnit,InE,SpecialString)
            % Convert between different energy units
            % Package: @convert
            % Description: Convert between different energy units.
            % Input  : - Input system:
            %            'erg'   - ergs
            %            'J'     - Jouls
            %            'Hz'    - Frequency [1/s]
            %            'A','Ang'- Wavelength [Ang]
            %            'cm'    - Wavelength [cm]
            %            'nm'    - Wavelength [nm]
            %            'm'     - Wavelength [m]
            %            'eV'    - Electron volts [h nu/q]
            %            'keV'   - kilo Electron volts [h nu/q]
            %            'MeV'   - Mega Electron volts [h nu/q]
            %            'GeV'   - Giga Electron volts [h nu/q]
            %            'T'     - Temperature [K]
            %            'me'    - Electron mass [E/m_e]
            %            'mp'    - Proton mass [E/m_p]
            %            'cal'   - calorie (4.184 J)
            %            'Btu'   - (1.055x10^3 J)
            %            'kWh'   - kilowatt-hour (3.6x10^6 J)
            %            'TNT'   - one ton of TNT (4.2x10^9 J) 
            %            'gr'    - Energy equivalent of 1 gram of matter (9x10^13 J)
            %          - Output system (see input system for options).
            %          - Input energy value to convert. Default is 1.
            %          - A string or a cell array of string that if the first input
            %            argument (In) is equal to one of them then the function will
            %            set the conversion factor to 1.
            % Output : - Energy in output system.
            % Tested : Matlab 7.0
            %     By : Eran O. Ofek                    Feb 2006
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: OutE=convert.energy('Hz','erg',1); % Convert 1Hz to ergs
            % Reliable: 2
            %--------------------------------------------------------------------------
            
            Def.InE           = 1;
            Def.SpecialString = [];
            if (nargin==2)
                InE           = Def.InE;
                SpecialString = Def.SpecialString;
            elseif (nargin==3)
                SpecialString = Def.SpecialString;
            elseif (nargin==4)
                % do nothing
            else
                error('Illegal number of input arguments: convert.energy(InUnits,OutUnits,Val,SpecialString)');
            end
                
            if any(strcmpi(InUnit,SpecialString))
                % user used Special string - set output to 1
                OutE = InE;
            else

                Erg2Erg = 1;
                Erg2J   = 1e-7;
                Erg2Hz  = 1.5092e26;
                Erg2A   = 1.9864e-8; 
                Erg2eV  = 6.2415e11;
                Erg2T   = 7.2430e15;
                Erg2me  = 1.2214e6;
                Erg2mp  = 665.214577;
                Erg2cal = Erg2J./4.184;
                Erg2Btu = Erg2J./1.055e3;
                Erg2kWh = Erg2J./3.6e6;
                Erg2TNT = Erg2J./4.2e9;
                Erg2gr  = constant.c.^-2; %('c','cgs').^-2;

                Relation = 'lin';

                switch lower(InUnit)
                 case 'erg'
                    ConvFactor = Erg2Erg;
                 case 'j'
                    ConvFactor = Erg2J;
                 case 'hz'
                    ConvFactor = Erg2Hz;
                 case {'a','ang'}
                    Relation   = 'inv';
                    ConvFactor = Erg2A;
                 case 'cm'
                    Relation   = 'inv';
                    ConvFactor = Erg2A.*1e-8;
                 case 'nm'
                    Relation   = 'inv';
                    ConvFactor = Erg2A.*1e-1;
                 case 'm'
                    Relation   = 'inv';
                    ConvFactor = Erg2A.*1e-10;
                 case 'micron'
                    Relation   = 'inv';
                    ConvFactor = Erg2A.*1e-4;
                 case 'ev'
                    ConvFactor = Erg2eV;
                 case 'kev'
                    ConvFactor = Erg2eV.*1e-3;
                 case 'mev'
                    ConvFactor = Erg2eV.*1e-6;
                 case 'gev'
                    ConvFactor = Erg2eV.*1e-9;
                 case 't'
                    ConvFactor = Erg2T;
                 case 'me'
                    ConvFactor = Erg2me;
                 case 'mp'
                    ConvFactor = Erg2mp;
                 case 'cal'
                    ConvFactor = Erg2cal;
                 case 'btu'
                    ConvFactor = Erg2Btu;
                 case 'kwh'
                    ConvFactor = Erg2kWh;
                 case 'tnt'
                    ConvFactor = Erg2TNT;
                 case 'gr'
                    ConvFactor = Erg2gr;
                 otherwise
                    error('Unknown InUnit Option');
                end

                switch Relation
                 case 'lin'
                    ErgE = InE./ConvFactor;
                 case 'inv'
                    ErgE = ConvFactor./InE;
                 otherwise
                    error('Unknown Relation Option');
                end


                Relation = 'lin';
                switch lower(OutUnit)
                 case 'erg'
                    ConvFactor = Erg2Erg;
                 case 'j'
                    ConvFactor = Erg2J;
                 case 'hz'
                    ConvFactor = Erg2Hz;
                 case {'a','ang'}
                    Relation   = 'inv';
                    ConvFactor = Erg2A;
                 case 'nm'
                    Relation   = 'inv';
                    ConvFactor = Erg2A.*1e-1;
                 case 'cm'
                    Relation   = 'inv';
                    ConvFactor = Erg2A.*1e-8;
                 case 'm'
                    Relation   = 'inv';
                    ConvFactor = Erg2A.*1e-10;
                 case 'micron'
                    Relation   = 'inv';
                    ConvFactor = Erg2A.*1e-4;
                 case 'ev'
                    ConvFactor = Erg2eV;
                 case 'kev'
                    ConvFactor = Erg2eV.*1e-3;
                 case 'mev'
                    ConvFactor = Erg2eV.*1e-6;
                 case 'gev'
                    ConvFactor = Erg2eV.*1e-9;
                 case 't'
                    ConvFactor = Erg2T;
                 case 'me'
                    ConvFactor = Erg2me;
                 case 'mp'
                    ConvFactor = Erg2mp;
                 case 'cal'
                    ConvFactor = Erg2cal;
                 case 'btu'
                    ConvFactor = Erg2Btu;
                 case 'kwh'
                    ConvFactor = Erg2kWh;
                 case 'tnt'
                    ConvFactor = Erg2TNT;
                 case 'gr'
                    ConvFactor = Erg2gr;
                 otherwise
                    error('Unknown InUnit Option');
                end


                switch Relation
                 case 'lin'
                    OutE = ErgE.*ConvFactor;
                 case 'inv'
                    OutE = ConvFactor./ErgE;
                 otherwise
                    error('Unknown Relation Option');
                end

             end
        end % convert.energy function
        
        % flux conversion
        function Out=flux(In,InUnits,OutUnits,Freq,FreqUnits)
            % Convert between different flux units
            % Package: @convert
            % Description: Convert between different flux units
            % Input  : - Flux
            %          - Input units - options are:
            %            'cgs/A'    [erg/s/cm^2/A]
            %            'cgs/Hz'   [erg/s/cm^2/Hz]
            %            'mJy'      [1e-26 erg/s/cm^2/Hz]
            %            'AB'       [AB magnitude]
            %            'STmag'    [ST magnitude]
            %            'ph/A'     [ph/s/cm^2/A]
            %            'ph/Hz'    [ph/s/cm^2/Hz]
            %          - Output units (see input units).
            %          - Frequency at which flux is measured.
            %          - Units of frequency - options are:
            %            'Hz'
            %            'A','Ang'
            %            'cm'
            %            'nm'
            %            'eV'
            %            'keV'...
            %            see convert.energy.m for more options
            % Output : - Flux in output units
            % Tested : Matlab 7.11
            %     By : Eran O. Ofek                    Mar 2011
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Updates: STmag option added on 2016-02-15
            % Example: Out=convert.flux(1,'mJy','AB',5000,'A')
            % Reliable: 2
            %---------------------------------------------------------------------------

            c = constant.c; %('c');  % speed of light [cm/s]
            h = constant.h; %('h');  % planck constant [cgs]


            Freq_Hz = convert.energy(FreqUnits,'Hz',Freq);

            switch lower(InUnits)
                case {'cgs/a','erg*cm^-2*s^-1*ang^-1','erg*cm^{-2}*s^{-1}*ang^{-1}'}
                   % convert to mJy
                   Flux_mJy = In.*(c.*1e8./(Freq_Hz.^2))./1e-26;
                case 'cgs/hz'
                   % convert to mJy
                   Flux_mJy = In./1e-26;
                case 'mjy'
                   % already in mJy
                   Flux_mJy = In;
                case 'ab'
                   % convert to mJy
                   Flux_mJy = 1e3.*10.^(23-(In+48.6)./2.5);
                case 'stmag'
                    % STmag to mJy
                    Flam = 10.^(-0.4.*(In+21.10));
                    Flux_mJy = Flam.*(c.*1e8./(Freq_Hz.^2))./1e-26;
                case 'ph/a'
                   Flux_mJy = In.*h.*Freq_Hz  .*(c.*1e8./(Freq_Hz.^2))./1e-26;
                case 'ph/hz'
                   Flux_mJy = In.*h.*Freq_Hz    ./1e-26;
                otherwise
                    error('Unknown InUnits option');
            end

            switch lower(OutUnits)
                case {'cgs/a','erg*cm^-2*s^-1*ang^-1','erg*cm^{-2}*s^{-1}*ang^{-1}'}
                   % convert mJy to cgs/A
                   Out = Flux_mJy./((c.*1e8./(Freq_Hz.^2))./1e-26);
                case 'cgs/hz'
                   % convert mJy to cgs/Hz
                   Out = Flux_mJy.*1e-26;
                case 'mjy'
                   Out = Flux_mJy;
                case 'ab'
                   Out = 2.5.*(23-log10(Flux_mJy.*1e-3))-48.6;
                case 'stmag'
                    Flam = Flux_mJy./((c.*1e8./(Freq_Hz.^2))./1e-26);
                    Out = -21.10 -2.5.*log10(Flam);
                case 'ph/a'
                   Out = Flux_mJy./(h.*Freq_Hz  .*(c.*1e8./(Freq_Hz.^2))./1e-26);
                case 'ph/hz'
                   Out = Flux_mJy./(h.*Freq_Hz    ./1e-26);
                otherwise
                    error('Unknown OutUnits option');
            end 

        end % convert.flux function
        
        % temperature conversion
        function OutTemp=temp(InTemp,InUnits,OutUnits)
            % Temperature conversion
            % Package: @convert
            % Description: Convert between temperature systems.
            % Input  : - Input temperature.
            %          - Units of the input temperature:
            %            'C' - degree Celsius
            %            'F' - degree Fahrenheit
            %            'K' - Kelvin
            %            'R' - Rankine
            %          - Units of output temperature.
            % Output : - Output temperature.
            % Tested : Matlab 7.6
            %     By : Eran O. Ofek                    Oct 2008
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: OutTemp=convert.temp(0,'C','K');
            % Reliable: 2
            %--------------------------------------------------------------------------


            switch InUnits
             case 'K'
                TempK = InTemp;
             case 'C'
                TempK = InTemp + 273.15;
             case 'F'
                TempK = (InTemp + 459.67).*5./9;
             case 'R'
                TempK = InTemp.*5./9;
             otherwise
                error('Unknown temperature units');
            end

            % from Kelvin to OutUnits

            switch OutUnits
             case 'K'
                OutTemp = TempK;
             case 'C'
                OutTemp = TempK - 273.15;
             case 'F'
                OutTemp = TempK.*9./5 - 459.67;
             case 'R'
                OutTemp = TempK.*9./5;
             otherwise
                error('Unknown temperature units');
            end

        end
        
    end
    
    methods (Static)
        
        % flux to luptitude
        function Lup=luptitude(Flux,Flux0,B)
            % Convert flux to luptitudes (asinh magnitudes)
            % Package: @convert
            % Description: Convert flux to luptitudes (asinh magnitudes).
            % Input  : - Flux.
            %          - Reference flux (Flux0), default is 1.
            %            Note that this parameter should equal to 10.^(0.4.*ZP);
            %          - Softening parameter (B), default is 1e-10.
            % Output : - Luptitude.
            % Tested : Matlab 7.13
            %     By : Eran O. Ofek                    Jul 2012
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Reliable: 2
            %--------------------------------------------------------------------------

            Def.Flux0 = 1;
            Def.B     = 1e-10;
            if (nargin==1)
               Flux0 = Def.Flux0;
               B     = Def.B;
            elseif (nargin==2)
               B     = Def.B;
            elseif (nargin==3)
               % do nothing
            else
               error('Illegal number of input arguments');
            end

            Lup = -2.5./log(10).*(asinh((Flux/Flux0)./(2.*B))+log(B));

        end % convert.luptitude function
        
        function Mag=flux2mag(Flux,ZP,Luptitude,Soft)
            % Convert flux to magnitude or luptitude
            % Package: @convert
            % Description: Convert flux to magnitude or luptitude. 
            % Input  : - Flux
            %          - ZP
            %          - A flag indicating if to return Luptitude (true) or
            %            magnitude (false). Default is true.
            %          - Luptitude softening parameter (B), default is 1e-10.
            % Output : - Magnitude or luptitude.
            % License: GNU general public license version 3
            % See also: luptitude.m
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    May 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Mag=convert.flux2mag(1,1)
            % Reliable: 2
            %--------------------------------------------------------------------------

            Def.ZP        = 0;
            Def.Luptitude = true;
            Def.Soft      = 1e-10;
            if (nargin==1)
                ZP        = Def.ZP;
                Luptitude = Def.Luptitude;
                Soft      = Def.Soft;
            elseif (nargin==2)
                Luptitude = Def.Luptitude;
                Soft      = Def.Soft;
            elseif (nargin==3)
                Soft      = Def.Soft;
            elseif (nargin==4)
                % do nothing
            else
                error('Illegal number of input arguments');
            end


            if (Luptitude)
                Mag = convert.luptitude(Flux,10.^(0.4.*ZP),Soft);
            else
                Mag = ZP - 2.5.*log10(Flux);
            end

        end
        
        function [SumMag,FracFlux]=sum_mag(Mag,Dim)      
            % Sum magnitudes
            % Package: @convert
            % Description: Sum a set of magnitudes.
            % Input  : - Matrix or vector of magnitudes to sum.
            %          - Dimesnions along to sum. Default is 1.
            % Output : - Sum of magnitude in mag.
            %          - Fraction of flux for each entry out of the total flux.
            % Tested : Matlab R2011b
            %     By : Eran O. Ofek                    Apr 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [SumMag,FracFlux]=convert.sum_mag([15;17;19])
            % Reliable: 2
            %--------------------------------------------------------------------------

            Def.Dim = 1;
            if (nargin==1)
                Dim = Def.Dim;
            end

            Flux     = 10.^(-0.4.*Mag);
            SumFlux  = nansum(Flux,Dim);
            SumMag   = -2.5.*log10(SumFlux);
            FracFlux = bsxfun(@rdivide,Flux,SumFlux);

        end
        
        % length conversion
%         function Out=length(In,Out,Val,SpecialString)
%             %--------------------------------------------------------------------------
%             % convert.length function                                   class/@convert
%             % Description: Convert length units.
%             % Input  : - Input length units name:
%             %            'm'|'cm'|'mm'|'micron'|'nm'|'ang'|'mile'|
%             %            Note that the 'pix' option generate a conversion factor of 1.
%             %          - Output angular units (options like first input argument).
%             %          - Value to convert. Default is 1.
%             %          - A string or a cell array of string that if the first input
%             %            argument (In) is equal to one of them then the function will
%             %            set the conversion factor to 1.
%             % Output : - Converted value (by default this is the conversion factor).
%             % License: GNU general public license version 3
%             % Tested : Matlab R2015b
%             %     By : Eran O. Ofek                    Jun 2016
%             %    URL : http://weizmann.ac.il/home/eofek/matlab/
%             % Example: convert.angular('rad','mas');
%             %          convert.angular('rad','mas',1,'pix');
%             % Reliable: 2
%             %--------------------------------------------------------------------------
% 
%             Def.Val = 1;
%             Def.SpecialString = [];
%             if (nargin==2),
%                 Val = Def.Val;
%                 SpecialString = Def.SpecialString;
%             elseif (nargin==3),
%                 SpecialString = Def.SpecialString;
%             elseif (nargin==4),
%                 % do nothing
%             else
%                 error('Illegal number of input arguments');
%             end
% 
% 
%             if any(strcmpi(In,SpecialString)),
%                 % user used Special string - set output to 1
%                 ConvVal = Val;
%             else
% 
%                 % Input: convert to deg
%                 switch lower(In)
%                     case {'rad','radian','radians'}
%                         ConvVal = 180./pi;
%                     case {'deg','degree','degrees'}
%                         ConvVal = 1;
%                     case 'arcmin'
%                         ConvVal = 1./60;
%                     case 'arcsec'
%                         ConvVal = 1./3600;
%                     case 'mas'
%                         ConvVal = 1./(3600.*1000);
%                     case 'frac'
%                         % fraction to deg
%                         ConvVal = 360;
%                     case 'hour'
%                         % hours to deg
%                         ConvVal = 15;
%                     otherwise
%                         error('Unknown input angular units');
%                 end
%                 ConvVal = Val.*ConvVal;
% 
%                 % Output: Convert from deg to requested output units
%                 switch lower(Out)
%                     case {'rad','radian','radians'}
%                         ConvVal = ConvVal.* pi./180;
%                     case {'deg','degree','degrees'}
%                         % do nothing
%                         %ConvVal = ConvVal;
%                     case 'arcmin'
%                         ConvVal = ConvVal.*60;
%                     case 'arcsec'
%                         ConvVal = ConvVal.*3600;
%                     case 'mas'
%                         ConvVal = ConvVal.*3600.*1000;
%                     case 'frac'
%                         ConvVal = ConvVal./360;
%                     case 'hour'
%                         ConvVal = ConvVal./15;
%                     otherwise
%                         error('Unknown input angular units');
%                 end
%             end         
%        end % convert.length function
        

    end % static
    
    % specific conversion functions
    methods (Static)
        
        % HMS to angle
        function OutHMS=hms2angle(InHMS,OutUnits)
            %  Convert Hour,Minutes,Seconds to angle.
            % Package: @convert
            % Description: Convert Hour,Minutes,Seconds to angle.
            % Input  : - A three column matrix of [H M S],
            %            or a sexagesimal string ('HH MM SS.SSS') with
            %            arbitrary delimiters or a cell array of
            %            sexagesimal strings.
            %          - Output units: 'rad'|'deg'|'hour'|frac'.
            %            Default is 'rad'
            % Output : - The output coordinates.
            % Example: convert.hms2angle([10 10 10])
            %          convert.hms2angle('10:10:10','deg')
            
            InvRAD = pi./180;
            
            if (nargin==1)
                OutUnits = 'rad';
            end
            
            if (isnumeric(InHMS))
                OutRad = (InHMS(:,1) + InHMS(:,2)./60 + InHMS(:,3)./3600).*15.*InvRAD;
            else
                if (ischar(InHMS))
                    InHMS = {InHMS};
                end
                Ncoo   = numel(InHMS);
                OutRad = zeros(size(InHMS)); 
                for Icoo=1:1:Ncoo
                    C = textscan(InHMS{Icoo},'%2d%*1c%2d%*1c%f');
                    OutRad(Icoo) = (double(C{1}) + double(C{2})./60 + C{3}./3600).*15.*InvRAD;
                end
            end
            
            % convert from radians to requested output units
            OutHMS = convert.angular('rad',OutUnits,OutRad);
            
        end % convert.hms2angle function
            
        % DMS to angle
        function OutDMS=dms2angle(InDMS,OutUnits)
            % Convert Sign,Deg,Minutes,Seconds to angle.
            % Package: @convert
            % Description: Convert Sign,Deg,Minutes,Seconds to angle.
            % Input  : - A four column matrix of [Sign D M S],
            %            or a sexagesimal string ('+DD MM SS.SSS') with
            %            arbitrary delimiters or a cell array of
            %            sexagesimal strings.
            %          - Output units: 'rad'|'deg'|'hour'|frac'.
            %            Default is 'rad'
            % Output : - The output coordinates.
            % Example: convert.hms2angle([-1 10 10 10])
            %          convert.hms2angle('-10:10:10','deg')
            
            InvRAD = pi./180;
            
            if (nargin==1)
                OutUnits = 'rad';
            end
            
            if (isnumeric(InDMS))
                OutRad = InDMS(:,1).*(InDMS(:,2) + InDMS(:,3)./60 + InDMS(:,4)./3600).*InvRAD;
            else
                if (ischar(InDMS))
                    InDMS = {InDMS};
                end
                Ncoo   = numel(InDMS);
                OutRad = zeros(size(InDMS)); 
                for Icoo=1:1:Ncoo
                    C = textscan(InDMS{Icoo},'%c%2d%*1c%2d%*1c%f');
                    Sign = 2.*(0.5-double(strcmp(C{1},'-')));  % 1 for -; 0 for other
                    OutRad(Icoo) = Sign.*(double(C{2}) + double(C{3})./60 + C{4}./3600).*InvRAD;
                end
            end
            
            % convert from radians to requested output units
            OutDMS = convert.angular('rad',OutUnits,OutRad);
            
        end % convert.dms2angle function
                
        
    end % Static
    
    % date and time conversion
    methods (Static)
        
        % date to ISO string of dates
        function DateStr=date2str(Date,SepSym)
            % Convert date [D M Y H M S] to an ISO date string
            % Package: @convert
            % Description: Convert date [D M Y H M S] to an ISO date string: YYYY-MM-DDTHH:MM:SS.FFF
            % Input  : - Matrix of dates [D M Y H M S] (H, M, S are optional).
            %            Each line corresponds to a single date.
            %          - Seperator character between the output date and
            %            time string. Default is 'T'.
            % Output : - A cell array of date string: YYYY-MM-DDTHH:MM:SS.FFF
            % Example: DateSte=convert.date2str([1 1 2000]);
            % Reliable: 2
            
            if (nargin==1)
                SepSym = 'T';
            end
            
            % add time to Date matrix
            [Ndate,Ncol] = size(Date);
            Date = [Date, zeros(Ndate,6-Ncol)];
            
            DateStr = cell(Ndate,1);
            for Idate=1:1:Ndate
                DateStr{Idate} = sprintf('%04d-%02d-%02d%c%02d:%02d:%06.3f',Date(Idate,[3 2 1]), SepSym, Date(Idate,[4 5 6]));
            end
            
        end % convert.date2str function
        
        function DateVec=str2date(String)
            % Convert a ISO time string  to matrix of times
            % Package: @convert
            % Description: Convert a string or a cell array of string containing date
            %              and time in the format 'YYYY-MM-DD HH:MM:SS.frac'
            %              or 'YYYY-MM-DD', to a matrix of dates with the following
            %              columns [Y M D H M S].
            % Input  : - A string or a cell array of string containing date
            %            and time in the format 'YYYY-MM-DD HH:MM:SS.frac'
            %            or 'YYYY-MM-DD'
            % Output : - A matrix of dates with the following columns [Y M D H M S].
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Dec 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: DateVec=convert.str2date({'2010-10-11T10:10:10.11','2014-12-12 10:23:59.1'});
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (ischar(String))
                String = mat2cell(String,ones(size(String,1),1));
            end
            if (size(String,1)==1)
                String = String.';
            end
            Time = regexp(String,'T|:|\s|-','split');
            
            DateVec = cell2mat(cellfun(@str2double,Time,'UniformOutput',false));

        end % convert.str2date function

	    function Frac=hour_str2frac(String)
            % Convert hour string to fraction of day
            % Package: @convert
            % Description: Convert a string or cell array of strings containing
            %              the hour in format HH:MM:SS.frac' to fraction of day.                                  % Input  : - String or cell array of strings containing the hour in                                   %            format HH:MM:SS.frac'.
            % Output : - Vector of fraction of day for each hour in the input.
            % Tested : Matlab R2014a
            %     By : Eran O. Ofek                    Dec 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/                                                   % Example: [Frac]=convert.hour_str2frac('10:10:10.111')
            %          [Frac]=convert.hour_str2frac({'10:10:10.111','11:11:11.12'})
            % Reliable: 2

            if (ischar(String))
                String = {String};
            end

            Time = regexp(String,':','split');
            FunH = @(C) (str2double(C{1}).*3600 + str2double(C{2}).*60 + str2double(C{3}))./86400;
            Frac = cellfun(FunH,Time);
        end

        % Convert between time systems
        function Output = time(Input, InType, OutType)
            % Convert between different types of time systems and years
            % Package: @convert
            % Description: Convert between different types of time systems and years.
            %              For example, this
            %              program can convert Julian years to Besselian years or JD
            %              and visa versa.
            % Input  : - Input to convert;
            %          - Input type - options are:
            %            'J'    - Julian year.
            %            'B'    - Besselian year.
            %            '1'    - Jan 1 of the year. (only input)
            %            'JD'   - Julian days.
            %            'MJD'  - MJD.
            %            'Date' - date matrix [D M Y] or [D M Y frac] or [D M Y H M S].
            %            'StrDate'- String date 'YYYY-MM-DDTHH:MM:SS.FFF',
            %                     or cell array of string dates.
            %            'StrDateO'- String date 'YYYY-MM-DD', or cell array of string dates.
            %          - Output type, options are as the input type.
            % Output : - Output.
            % Tested : Matlab 7.10
            %     By : Eran O. Ofek                    Oct 2010
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: convert.time(2000,'J','B')
            % Reliable: 2
            %--------------------------------------------------------------------------

            switch lower(InType)
                case 'j'
                    % convert Julian years to JD
                    JD = (Input - 2000).*365.25 + 2451545.0;
                case 'b'
                    % convert Besselian years to JD
                    JD = (Input - 1900).*365.2421988 + 2415020.3135;
                case '1'
                    % convert Jan 1st Year to JD
                    N = numel(Input);
                    JD = convert.date2jd([ones(N,2), Input(:)]);
                case 'jd'
                    JD = Input;
                case 'mjd'
                    JD = Input + 2400000.5;
                case 'date'
                     JD = convert.date2jd(Input);
                case 'strdate'
                     JD = convert.date2jd(Input);
                case 'strdateo'
                     JD = convert.date2jd(Input);     
                otherwise
                 error('Unknown InType option');
            end


            switch lower(OutType)
                case 'j'
                    % convert JD to Julian years
                    Output = 2000 + (JD-2451545.0)./365.25;
                case 'b'
                    % convert JD to Besselian years
                    Output = 1900 + (JD-2415020.3135)./365.2421988;
                case 'jd'
                    Output = JD;
                case 'mjd'
                    Output = JD - 2400000.5;
                case 'date'
                    Output = convert.jd2date(JD);
                case 'strdate'
                    Output = convert.date2str(convert.jd2date(JD,'H'));
                case 'strdateo'
                    Output = convert.date2str(convert.jd2date(JD,'H'));  
                    Nout= numel(Output);
                    for Iout=1:1:Nout
                        Output{Iout} = Output{Iout}(1:10);
                    end
             otherwise
                error('Unknown OutType option');
            end

        end % convert.time function
        
        % Date to JD
        function JD=date2jd(Date,Output,TreatOnlyDate)
            % Convert Julian/Gregorian date to Julian Day
            % Package: @convert
            % Description: Convert Julian/Gregorian date to Julian Day.
            % Input  : - Gregorian of Julian date in one of the following formats
            %            [Day, Month, Year, Day_Fraction]
            %            or [Day, Month, Year, Hour, Min, Sec]
            %            or [Day, Month, Year] - in this case set Day_Fraction to 0.
            %            Alternatively this can be a string or a cell array of strings
            %            in which each string contains the date in the format:
            %            'yyyy-mm-ddTHH:MM:SS' (e.g., '2010-08-11T15:01:56') or:
            %            'yyyy-mm-dd HH:MM:SS'.
            %            If argument is not provided then the program will calculate
            %            the JD for now using the clock UTC computer (time zone
            %            included).
            %          - Output type. Options are:
            %            'JD'  - Julian days (default).
            %            'MJD' - Modified JD (JD-2400000.5).
            %          - Flag indicating if possible to treat strings
            %            containing only dates (no H:M:S). Default is false.
            % Output : - Row vector of Julian days.
            % Tested : Matlab 3.5
            %     By : Eran O. Ofek                    Jan 1994
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: convert.date2jd([1 1 2000 10 30 0]);
            %          convert.date2jd([1 1 2000; 2 2 2000]);
            %          convert.date2jd('2010-10-12 10:10:10.111');
            %          convert.date2jd({'2010-10-12 10:10:10.111'});
            %          convert.date2jd;  % JD of now
            % Reliable: 1
            %--------------------------------------------------------------------------

            if nargin<3
                TreatOnlyDate = true;
            end
            
            IsStr = false;
            if (nargin==0)
               %Date = clock;
               Date = datevec(datetime('now', 'TimeZone', 'UTC'));
               IsStr = true;
               
            end
            if (isempty(Date))
               %Date = clock;
               Date = datevec(datetime('now', 'TimeZone', 'UTC'));
               IsStr = true;
               
            end
            if (ischar(Date) || iscell(Date))
                Date=convert.str2date(Date);
                IsStr = true;
                
            end
            if any(isnan(Date))
                Date = nan(1,6);
            end
            
            if IsStr
                if size(Date,2)==6
                    Date = Date(:,[3 2 1 4 5 6]);
                else
                    if TreatOnlyDate
                        Date = [Date(:,[3 2 1]), zeros(size(Date,1),3)];
                    else
                        Date = nan(size(Date,1),6);
                    end
                end
            end
            
            Y = Date(:,3);
            M = Date(:,2);
            D = Date(:,1);

            [Lines, Rows] = size(Date);
            switch Rows
             case 3
                F = zeros(Lines,1);
             case 4
                F = Date(:,4);
             case 6
                F = convert.hms2frac(Date(:,4:6));
             otherwise
                error('Illegal number of column in Date matrix');
            end

            B  = zeros(Lines,1);
            Im3 = find(M<3);
            Y(Im3) = Y(Im3) - 1;
            M(Im3) = M(Im3) + 12;

            Iy = find(Y>1582 | (Y==1582 & M>10) | (Y==1582 & M==10 & D>=15));
            A = floor(Y(Iy).*0.01);
            B(Iy) = 2 - A + floor(0.25.*A);
            JD = floor(365.25.*(Y + 4716)) + floor(30.6001.*(M + 1)) + D + B - 1524.5 + F;

            if (nargin>1)
               switch lower(Output)
                case 'jd'
                   % do nothing
                case 'mjd'
                   JD = JD - 2400000.5;
                otherwise
                   error('Unknown Output option');
               end
            end

        end % jd function
        
        % JD to date
        function Date=jd2date(JD,Format)
            % Convert Julian days to Gregorian/Julian date
            % Package: @convert
            % Description: Convert Julian days to Gregorian/Julian date.
            % Input  : - Row vector of (positive) Julian Days.
            %          - Output format:
            %            'f'  - [Day Month Year, Day_Fraction(UT)] (default).
            %            'H'  - [Day Month Year, H M S]
            % Output : - Matrix of dates.
            %            e.g., [Day, Month, Year, Day_Fraction(UT)].
            % Tested : Matlab 5.2
            %     By : Eran O. Ofek                    Sep 1999
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: convert.jd2date(convert.date2jd([1 1 2000; 2 2 2000]).')
            % Reliable: 1
            %--------------------------------------------------------------------------
            if (min(JD)<0)
               error('The method is valid only for poitive JDs');
            end

            if (nargin==1)
               Format = 'f';
            end

            Z = floor(JD+0.5);
            F = (JD+0.5) - floor(JD+0.5);

            A     = zeros(size(JD));
            %Day   = zeros(size(JD));
            Month = zeros(size(JD));
            Year  = zeros(size(JD));

            IZ_s = find(Z<2299161);
            IZ_l = find(Z>=2299161);
            A(IZ_s) = Z(IZ_s);

            Alpha = fix((Z(IZ_l) - 1867216.25)./36524.25);
            A(IZ_l) = Z(IZ_l) + 1 + Alpha - fix(Alpha./4);

            B = A + 1524;
            C = fix((B - 122.1)./365.25);
            D = fix(365.25.*C);
            E = fix((B-D)./30.6001);

            Day   = B - D - fix(30.6001.*E) + F;
            IE_s  = find(E<14);
            IE_l  = find(E>=14);
            Month(IE_s) = E(IE_s) - 1;
            Month(IE_l) = E(IE_l) - 13;

            IM_l  = find(Month>2);
            IM_s  = find(Month==1 | Month==2);
            Year(IM_l) = C(IM_l) - 4716;
            Year(IM_s) = C(IM_s) - 4715;

            Date = [floor(Day), Month, Year, F];

            switch lower(Format)
             case 'f'
                % do nothing
             case 'h'
                %Date = [Date(:,1:3), convertdms(Date(:,4),'f','H')];
                Date = [Date(:,1:3), convert.frac2hms(Date(:,4))];
             otherwise
                error('unknown Format option');
            end
            
        end % convert.jd2date function
        
        % fraction to [H M S]
        function HMS=frac2hms(Frac)
            % Convert fraction to [H M S]
            % Package: @convert
            % Description: Convert fraction to [H M S].
            %              Note that if the seconds or minutes are printed
            %              as 60 it means that they are actually 60-eps.
            % Input  : - Vector of fractions.
            % Output : - A three column matrix of [H M S].
            % Example: frac2hms(rand(10,1))
            % Reliable: 2
            
            Hour = Frac(:).*24;
            H    = floor(Hour);
            Min  = (Hour - H).*60;
            M    = floor(Min);
            S    = (Min - M).*60;
            
%             Flag = S>59.99999999 || M>59.99999999;
%             Frac(Flag) = Frac(Flag) + 1e-15;
%             
%             Hour = Frac(:).*24;
%             H    = floor(Hour);
%             Min  = (Hour - H).*60;
%             M    = floor(Min);
%             S    = (Min - M).*60;
            
            HMS  = [H M S];
            
        end % convert.frac2hms
        
        % [H M S] to fraction
        function Frac=hms2frac(HMS)
            % Convert [H M S] to fraction
            % Package: @convert
            % Description: Convert [H M S] to fraction.
            % Input  : - A three column matrix of [H M S].
            % Output : - A fraction of day.
            % Example: convert.hms2frac([10 10 10]);
            % Reliable: 2
            
            Frac = (HMS(:,1).*3600 + HMS(:,2).*60 + HMS(:,3))./86400;
            
        end % convert.hms2frac function
        
        
    end % Static
    
    
    methods (Static)
        function Result = conversion(InUnits, OutUnits, Val)
            %
           
            % length: ang, nm, mm, cm, m, km, inch, mile, feet, yard, ly,
            % pc, kpc, Mpc, Gpc
            
            % mass: me, mp, g, kg, EM, JM, SM
            
            % time: s, min, hour, day, JY
            
            % Temp: 
            
        end
    end


    % Unit test
    methods(Static)   
        Result = unitTest()

    end    	

end % end class
            
