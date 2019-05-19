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

classdef QHY367 < handle
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
        
        function QC=QHY367
            % QHY367 class constractor
            % Open link to the camera an initilaize the camera
            % Example: 
            
            
            
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
        function Flag=open(QHY)
            % open camera link
        end
        
        function Flag=close(QHY)
            
        end
        
        function Flag=set_temp(QHY,Temp)
            
        end
        
        function [Temp,Flag]=get_temp(QHY,Temp)
            
        end
        
        function Flag=set_exptime(QHY,ExpTime)
            % ExpTime in seconds
            
        end
        
        function Flag=set_gain(QHY,Gain)
            
        end
        
        function Flag=get_gain(QHY,Gain)
            
        end
        
        function Flag=set_binning(QHY,Binning)
            % default is 1x1
            
        end
        
        function Flag=set_bitdepth(QHY,BitDepth)
            % default is 16bit
            
        end
        
        function Flag=set_color(QHY,ColorMode)
            % default is bw
            
        end
        
        function [Image,Flag]=take_image(QHY)
            % take single image
            % no display
            
        end
        
        function [NImage,Flag]=take_live(QHY,Nimages)
            % take live image
            % currently set exptime to 2s.
            % Nimages default is 10.
            % Nimages of size [X, Y, N]
            
            
            
            
        end
        
        
        
    end
    
    
end

            
