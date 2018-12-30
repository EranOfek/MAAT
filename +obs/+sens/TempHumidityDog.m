% TempHumidityDog class  
% Package: +obs/+sens/
% Description: A class for controlling the USB-TnH Type A (SHT10)
%              USB temperature + humidity sensor from DogRatian's.
%              See http://www.dogratian.com/products/index.php/menu-sensors/menu-usb-tnh-type-a-sht10
%              The class open a serial object to communicate with the
%              sensor.
% Tested : Matlab R2018a
%     By : Eran O. Ofek                    Dec 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: % Create a TempHumidityDog object
%          WC=obs.sens.TempHumidityDog('COM');
%          % update wind temp/humidity values in WC object properties
%          getWind(WC);
%          % close/open WC.ComObj tcpip object
%          open(WC); close(WC);
% Reliable: 2
%--------------------------------------------------------------------------

classdef TempHumidityDog < handle
    properties (SetAccess = public)
        % generic fields
        Status         = false;                         % false - readings are unreliable, true - ok
        Data           @ stack                          % Stack object containing data history
        DataCol        = {'JD','Temp','Humidity'};      % Stack object columns
        DataUnits      = {'day','C',''};                % Stack object column units
         
        
        % specific fields
        Temp             = NaN;                        % last Temp
        TempUnits        = 'C';                        % Temp units
        Humidity         = NaN;                        % last humidity
        HumidityUnits    = '';                         % humidity units
        LastJD           = NaN;                        % JD of last sucessful reading
        
    end
    
    properties (Constant = true)
        
    end
    
    properties (Hidden = true)
        
        Protocol
        Com 
        ComObj 
        RetLine
        TempStatus
        HumidityStatus
        
    end
    
    
    % Constructor
    methods
        
        function WC=TempHumidityDog(Com)
            % TempHumidityDog class constractor
            % Example: WC=obs.sens.TempHumidityDog('132.77.38.187',10001)
            BaudRate = 115200;
            
            T = serial(Com);
            T.BaudRate = BaudRate;
            
            WC.Data     = stack(nan(100,3));
            WC.Com      = Com;
            WC.Protocol = 'serial';
            WC.ComObj   = T;
            
            %WC.ComObj.Timeout = 10;
            WC.ComObj.Terminator = 'LF';  % set terminaotor to carige return
            
            % if close then open
            switch lower(T.Status)
                case 'closed'
                    fopen(T);
            end
            
            WC.update;
            
        end
        

    end
    
    % getters/setters
    methods
%         function WC=get.Data(WC)
%             % update the Data field (stack)
%             
%             getWind(WC);
%             
%         end
      
    end
    
    
    % static methods
    methods (Static)
     
    end
    
    methods 
        
        function open(WC)
            % Open tcp/ip connection to WindETH
            
            fopen(WC.ComObj);
            
        end
        
        function fopen(WC)
            % Open tcp/ip connection to WindETH
            
            fopen(WC.ComObj);
            
        end
        
        function close(WC)
            % Close tcp/ip connection to WindETH
            
            fclose(WC.ComObj);
            
        end
        
        function fclose(WC)
            % Close tcp/ip connection to WindETH
            
            fclose(WC.ComObj);
            
        end
        
        function [Line]=update(WC)
            % Read temperature and humidity
            % Package: obs.sense.TempHumidityDog
            % Description: Return observatories parameters
            % Input  : *
            % Output : - The return line.
            % Example: WC=obs.sens.TempHumidityDog('COM'); getWind(WC)
           
            % If tcp/ip connection is closed - than open
            switch lower(WC.ComObj.Status)
                case 'closed'
                    fopen(WC.ComObj);
            end
            
            % get current time [JD]
            JD = celestial.time.julday;
            
            % send query to WindETH and get answer
            %fwrite(WC.ComObj,unicode2native(sprintf('*B1MR0\r'), 'UTF-8'))
            %V=fread(WC.ComObj);
            % for query see Papouch s.r.o Spinel in AD4xxx manual: pp 10-11
            WC.RetLine{1} = query(WC.ComObj,sprintf('GT\n'));
            LineT = WC.RetLine{1};
            WC.RetLine{2} = query(WC.ComObj,sprintf('GH\n'));
            LineH = WC.RetLine{2};
            
            WC.Temp    = str2double(LineT);
            WC.Humidity = str2double(LineH);
            if (isnan(WC.Temp))
                WC.TempStatus = 'error';
            else
                WC.TempStatus = 'ok';
            end
            if (isnan(WC.Humidity))
                WC.HumidityStatus = 'error';
            else
                WC.HumidityStatus = 'ok';
            end
            
            % set Wind speed/dir status
            Status = [false false];
            switch WC.TempStatus
                case 'ok'
                    Status(1) = true;
            end
            switch WC.HumidityStatus
                case 'ok'
                    Status(2) = true;
            end
            
            % set status to true only if both speed and dir are ok
            if (all(Status))
                WC.Status = true;
            end
            
            % if status is ok than set the time of the last query
            if (WC.Status)
               % update time to that of last measurment
               WC.LastJD = JD;
            end
            
            
            WC.Data.add([WC.LastJD, WC.Temp, WC.Humidity]);
            
        end
        
    end
    
    
    
end

            
