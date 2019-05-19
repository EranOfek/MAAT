% WeatherBoltwood class  
% Package: +obs/+sens/
% Description: A class for controlling the Boltwood meteorological sensor
%              (from http://diffractionlimited.com/product/boltwood-cloud-sensor-ii/).
% Tested : Matlab R2018a
%     By : David Polishook                    Dec 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: % Create a Boltwood object
%          WB=actxserver('ClarityII.CloudSensorII');
%          % update weather station values in WB object properties
%          update(WB) - update all parameters from theweather program
%          install ASCOM;
%          
% Reliable: 2
%--------------------------------------------------------------------------

classdef WeatherBoltwood < handle
    properties (SetAccess = public)
        % generic fields
        Status         = false;                         % false - readings are unreliable, true - ok
        Data           @ stack                          % Stack object containing data history
        DataCol        = {'JD','DayLightV','AmbientT','RelSkyT','SensorT','DewPointT','WindSpeed','Humidity','Rain'};   % Stack object columns
        DataUnits      = {'day','relative','deg C',   'deg C',  'deg C',  'deg C',    'km/h',     'percentage','mm'};          % Stack object column units

        % specific fields
        DayLightV      = NaN;                        % Day light value
        AmbientT       = NaN;                        % Ambient temperature
        AmbientTUnit   = 'deg C';                    % Ambient temperature units
        RelSkyT        = NaN;                        % Sky minus ambient temperature; an
                                                     % indicator for the clouds condition
        RelSkyTUnit    = 'deg C';                    % Sky minus ambient temperature units
        SensorT        = NaN;                        % Sensor temperature
        SensorTUnit    = 'deg C';                    % Sensor temperature unit
        DewPointT      = NaN;                        % Dew point temperature
        DewPointTUnits = 'deg C';                    % Dew point temperature units
        WindSpeed      = NaN;                        % Wind speed
        WindSpeedUnits = 'km/h';                     % wind speed units
        Humidity       = NaN;                        % Humiditiy
        HumidityUnits  = 'percentage';               % Humiditiy units
        Rain           = NaN;                        % Rain
        RainUnits      = 'mm';                       % Rain

        LastJD         = NaN;                        % JD of last sucessful reading

        % Text conditions from Boltwood
        DayCondition   = NaN;                        % Day condition
        CloudCondition = NaN;                        % Cloud condition
        WindCondition  = NaN;                        % Wind condition
        RainCondition  = NaN;                        % Rain condition        
    end
    
    properties (Constant = true)
        
    end
    
    properties (Hidden = true)
        
        ComObj 
        
    end
    
    
    % Constructor
    methods
        
        function WB=WeatherBoltwood()
            % WeatherBolton class constractor
            % Example: WB=obs.sens.WeatherBoltwood
            
            WB.Data     = stack(nan(100,numel(WB.DataCol)));
            WB.open              % Open Boltwood weather application.
                                 % Also open connection to application if
                                 % already opened.
            WB.update;           % update data from weather application
        end
        

    end
    
    % getters/setters
    methods

      
    end
        
    % static methods
    methods (Static)
        
    end
        
    methods
        function open(WB)
            % Open Boltwood weather application. Also open connection to
            % application if already opened.
            WB.ComObj = actxserver('ClarityII.CloudSensorII');
        end

        function update(WB)
            % Read weather parameters
            % Package: obs.sens.WeatherBoltwood
            % Description: Return observatories parameters
            % Input  : *
            % Output : - The return line.
            % Example: WB=obs.sens.WeatherBoltwood; update;

            % get current time [JD]
            WB.LastJD = celestial.time.julday;
            
            % If the weather application is down then reinitiate it
            try
                WB.ComObj.DataReady;
            catch
%                obs.sens.WeatherBoltwood(false); % call constructor without update
                pause(5)
                WB.Status = false;
                WB.open;
            end

            % send query to Boltwood and get answer
            WB.DayLightV      = WB.ComObj.DayLightV;          % Day light value
            WB.AmbientT       = WB.ComObj.AmbientT;           % Ambient temperature
            WB.RelSkyT        = WB.ComObj.RelSkyT;            % Sky minus ambient temperature; an
                                                    % indicator for the clouds condition
            WB.SensorT        = WB.ComObj.SensorT;            % Sensor temperature
            WB.DewPointT      = WB.ComObj.DewPointT;          % Dew point temperature
            WB.WindSpeed      = WB.ComObj.wind;               % Wind speed
            WB.Humidity       = WB.ComObj.HumidityPercent;    % Humiditiy percentage
            WB.Rain        = WB.ComObj.RainF;              % Rain

            % Text conditions from Boltwood
            WB.DayCondition   = WB.ComObj.DayCondition;       % Day condition
            WB.CloudCondition = WB.ComObj.CloudCondition;     % Cloud condition
            WB.WindCondition  = WB.ComObj.WindCondition;      % Wind condition
            WB.RainCondition  = WB.ComObj.RainCondition;      % Rain condition
            
            WB.Data.add([WB.LastJD, WB.DayLightV, WB.AmbientT, WB.RelSkyT, WB.SensorT, WB.DewPointT, ...
                         WB.WindSpeed, WB.Humidity, WB.Rain]);
            
    
            % if status is ok than set the time of the last query
            Dummy = WB.ComObj.DataReady;
            if (Dummy),
               % Update Status
               WB.Status = true;
            else
               WB.Status = false;
            end
             
        end
        
    end
    
    
    
end

            
