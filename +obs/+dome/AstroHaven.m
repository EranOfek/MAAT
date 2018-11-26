% AstroHaven class  
% Package: +obs/+dome/
% Description: A class for controlling the AstroHaven dome.
%              The class open a serial port to the instrument
% Tested : Matlab R2018a
%     By : David Polishook                    Nov 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: %
% Reliable: ?
%--------------------------------------------------------------------------

classdef AstroHaven < handle
    properties
        % generic fields
        Status         = false;                         % false - readings are unreliable, true - ok
        Shutter1
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
            
            WC.Data     = stack(nan(100,3));
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
    
end
