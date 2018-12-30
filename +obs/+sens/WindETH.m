% WindETH class  
% Package: +obs/+sens/
% Description: A class for controlling the WindETH ethernet anemometer
%              (from papouch.com).
%              https://www.papouch.com/en/shop/product/windeth-ethernet-anemometer/
%              The class open a tcp/ip port to the instrument
% Tested : Matlab R2018a
%     By : Eran O. Ofek                    Nov 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: % Create a WindETH object
%          WC=obs.sens.WindETH('132.77.38.187',10001);
%          % update wind speed/az values in WC object properties
%          getWind(WC);
%          % close/open WC.ComObj tcpip object
%          open(WC); close(WC);
% Reliable: 2
%--------------------------------------------------------------------------

classdef WindETH < handle
    properties (SetAccess = public)
        % generic fields
        Status         = false;                         % false - readings are unreliable, true - ok
        Data           @ stack                          % Stack object containing data history
        DataCol        = {'JD','WindSpeed','WindAz'};   % Stack object columns
        DataUnits      = {'day','km/h','deg'};          % Stack object column units
        
        
        % specific fields
        WindSpeed      = NaN;                        % last wind speed
        WindSpeedUnits = 'km/h';                     % wind speed units
        WindAz         = NaN;                        % last wind Az
        WindAzUnits    = 'deg';                      % wind Az units
        LastJD             = NaN;                    % JD of last sucessful reading
        
    end
    
    properties (Constant = true)
        
    end
    
    properties (Hidden = true)
        
        Protocol
        IP       
        Port    
        ComObj 
        RetLine
        WindAzStatus
        WindSpeedStatus
        
    end
    
    
    % Constructor
    methods
        
        function WC=WindETH(IP,Port)
            % WindETH class constractor
            % Example: WC=obs.sens.WindETH('132.77.38.187',10001)
            
            T = tcpip(IP,Port);
            
            WC.Data     = stack(nan(100,numel(WC.DataCol)));
            WC.IP       = IP;
            WC.Port     = Port;
            WC.Protocol = 'tcpip';
            WC.ComObj   = T;
            
            %WC.ComObj.Timeout = 10;
            WC.ComObj.Terminator = 'CR';  % set terminaotor to carige return
            
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
        function Az=dir2az(DirStr)
            % Convert direction (e.g., NEE) to Azimuth
            % Input  : Direction (e.g., N, NE, NNW)
            % Output : Azimuth (deg)
            % Example: obs.sens.WindETH.dir2az('NE')
            
            N_N = numel(strfind(DirStr,'N'));
            N_E = numel(strfind(DirStr,'E'));
            N_S = numel(strfind(DirStr,'S'));
            N_W = numel(strfind(DirStr,'W'));
            
            N = [N_N N_E N_S N_W];
            
             
            
            switch sprintf('%d%d%d%d',N)
                case '1000'
                    Az = 0;
                case '0100'
                    Az = 90;
                case '0010'
                    Az = 180;
                case '0001'
                    Az = 270;
                case '2100'
                    Az = 22.5;
                case '1100'
                    Az = 45;
                case '1200'
                    Az = 67.5;
                case '0210'
                    Az = 112.5;
                case '0110'
                    Az = 135;
                case '0120'
                    Az = 157.5;
                case '0021'
                    Az = 202.5;
                case '0011'
                    Az = 225;
                case '0012'
                    Az = 247.5;
                case '1002'
                    Az = 292.5;
                case '1001'
                    Az = 315;
                case '2001'
                    Az = 337.5;
                otherwise
                    Az = NaN;
             end
            
        end
        
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
            % Read wind speed and direction
            % Package: obs.sense.WindETH
            % Description: Return observatories parameters
            % Input  : *
            % Output : - The return line.
            % Example: WC=obs.sens.WindETH('132.77.38.187',10001); getWind(WC)
           
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
            WC.RetLine = query(WC.ComObj,sprintf('*B1MR0\r'));
            Line = WC.RetLine;
            
            
            % single measurment command
%             HexCommand = ['2A';'61';'00';'06';'31';'02';'51';'00';'EA';'0D'];
%             DecCommand = hex2dec(HexCommand); 
%             fwrite(WC.ComObj,DecCommand,'uint8');
%             DecAns     = fread(WC.ComObj);
%             HexAns     = dec2hex(DecAns); 
%             
%             HexAns(7,:)     % command recieved OK
%             HexAns(8,:)     % status for ch1
%             HexAns(9:10,:)  % ch1 value
%             HexAns(12,:)    % status for ch2
%             HexAns(13:14,:) % ch2 value
%             
            % Extract data from output string
            %WC.RetLine
            Cell = regexp(WC.RetLine,'\s','split');
            Flag = cellfun(@isempty,Cell);
            Cell = Cell(~Flag);
            
            WindAzStatus    = Cell{3};   % WindAz status flag [80 | 88]
            WindAzString    = Cell{4};   % WindAz direction string
            WindSpeedStatus = Cell{6};   % WindSpeed status flag [80 | 88]
            WindSpeedString = Cell{7};   % WindSpeed 
            
            
            % set Wind speed/dir status
            Status = [false false];
            switch WindAzStatus
                case '80'
                    Status(1) = true;
                    WC.WindAzStatus = 'ok';
                case '88'
                    WC.WindAzStatus = 'out of range';
                otherwise
                    WC.WindAzStatus = 'error';
            end
            switch WindSpeedStatus
                case '80'
                    Status(2) = true;
                    WC.WindSpeedStatus = 'ok';
                case '88'
                    WC.WindSpeedStatus = 'out of range';
                otherwise
                    WC.WindSpeedStatus = 'error';
            end
            
            % set status to true only if both speed and dir are ok
            if (all(Status))
                WC.Status = true;
            end
            
            % convert speed/dir to numeric values
            if (WC.Status)
                % if status ok - convert wind Az/spped to numeric values
                WC.WindSpeed = str2double(WindSpeedString);
                
                Az = obs.sens.WindETH.dir2az(WindAzString);
                if (isnan(Az))
                    WC.Status = false;
                end
                WC.WindAz    = Az;
                
            end
            
            % if status is ok than set the time of the last query
            if (WC.Status)
               % update time to that of last measurment
               WC.LastJD = JD;
            end
            
            
            WC.Data.add([WC.LastJD, WC.WindSpeed, WC.WindAz]);
            
            
            %WC.WindSpeedHistory(end+1) = WC.WindSpeed;
            %WC.WindAzHistory(end+1)    = WC.WindAz;
            %WC.TimeHistory(end+1)      = WC.Time;
             
             
             
        end
        
    end
    
    
    
end

            
