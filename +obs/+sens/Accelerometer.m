classdef Accelerometer < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        hndl@Bluetooth; % bluetooth serial object
        
        Data@stack;
        
    end
    
    properties % inputs/outputs
        
        DataCol = {'X', 'Y', 'Z'};
        
        % acc_vec = acc_vec_data/gain + bias
        gain; % x,y,z
        bias; % x,y,z
        
        calibration_data;
        
    end
    
    properties % switches/controls
        
        bluetooth_name = 'HC-06';
        bluetooth_id = ''; 
        
        status = 0;
        time;
        jd;
        reply = '';
        acc_vec_raw;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        acc_vec;
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Accelerometer(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.sens.Accelerometer')
                if obj.debug_bit, fprintf('Accelerometer copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Accelerometer constructor v%4.2f\n', obj.version); end
                                
                obj.connect(varargin{:});
                
                obj.reset;
                
            end
            
        end
        
        function connect(obj, name, id)
            
            if nargin>1 && ~isempty(name)
                obj.bluetoot_name = name;
            end
            
            if nargin>2 && ~isempty(id)
                obj.bluetooth_id = id;
            end
            
            % first, make sure to close existing connections...
            try
                fclose(obj.hndl);
                delete(obj.hndl);
            end
            
            if isempty(obj.bluetooth_name)
                error('Must supply a name for a bluetooth device (e.g., HC-06)');
            end
            
            % must be paired to the bluetooth device! 
            obj.hndl = Bluetooth(obj.bluetooth_name, 1); % second argument is channel==1
            
            if isempty(obj.bluetooth_id) || isnumeric(obj.bluetooth_id)
                
                in = instrhwinfo('bluetooth', obj.bluetooth_name);
                
                if isnumeric(obj.bluetooth_id) && obj.bluetooth_id>0
                    idx = obj.bluetooth_id;
                else
                    idx = 1;
                end
                
                if idx>length(in)
                    error('Cannot open device number %d, there are only %d devices named "%s".', idx, length(in), obj.bluetooth_name);
                end
                
                if isempty(in(idx).RemoteID)
                    error('Device not found. Make sure to pair the device!');
                end
                
                if obj.debug_bit
                    fprintf('Found a bluetooth device with ID: %s\n', obj.bluetooth_id);
                end
                
                obj.bluetooth_id = in(idx).RemoteID(9:end);
                
            end
            
            obj.hndl.RemoteID = obj.bluetooth_id; 
            
            fopen(obj.hndl);
            
            pause(0.1);
            
            obj.update;
            
        end
        
    end
    
    methods % reset/clear
        
        function reset_calibration(obj)
            
            obj.calibration_data = [];
            obj.gain = [];
            obj.bias = [];
            
        end
        
        function reset(obj)
            
            obj.Data = stack(nan(10,3));
            
        end
        
    end
    
    methods % getters
        
        function val = get.acc_vec(obj)
            
            if isempty(obj.gain) || isempty(obj.bias)
                val = obj.acc_vec_raw;
            else
                val = (obj.acc_vec_raw - obj.bias)./obj.gain;
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % commands
        
        function update(obj)
            
            obj.status = 0;
            
            if ~strcmp(obj.hndl.Status, 'open')
                error('Device is closed, use fopen or connect function');
            end
            
            obj.hndl.BytesAvailableFcn = @obj.read_data;
            
            fprintf(obj.hndl, 'status;');
            
%             res = 0.001;
%             timeout = 5;
%             
%             if obj.hndl.BytesAvailable % flush any existing data from before
%                 fgetl(obj.hndl);
%             end
%             
%             for ii = 1:timeout/res
%                
%                 if obj.hndl.BytesAvailable
%                     
%                     obj.read_data;
%                     
%                     obj.status = 1;
%             
%                     return;
%                     
%                 end
%                 
%                 pause(res);
%                 
%             end
%             
%             error('Timeout reached after %f seconds while waiting for Bluetooth response...', timeout);
            
        end
        
        function read_data(obj, ~, ~)
            
            obj.reply = fgetl(obj.hndl);
            obj.acc_vec_raw = str2double(regexp(obj.reply,'-?\d*','Match'));
            obj.time = datetime('now', 'timezone', 'UTC');
            obj.jd = juliandate(obj.time);
            obj.status = 1;
            
            obj.Data.add(obj.acc_vec);
            
        end
        
        function measureCalibrationData(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('plotting', 0);
            input.input_var('time', 10, 'duration');
            input.input_var('N', 1000, 'number', 'measurements');
            input.scan_vars(varargin{:});
            
            obj.reset_calibration;
            
            for ii = 1:input.N
                
                obj.update;
                
                obj.calibration_data(ii,:) = obj.acc_vec;
                
                pause(input.time./input.N);
                
                if input.plotting
                    obj.showCalibration;
                    drawnow;
                end
                
            end
            
        end
        
        function sum_err_square = runCalibration(obj, varargin)
            
            data = obj.calibration_data;
            
            B_initial = mean(data);
            G_initial = mean(sqrt(sum((data-B_initial).^2,2)))*[1 1 1]; % assume equal gain at first

            b_initial = [B_initial, G_initial];

            func = @(b) obs.sens.Accelerometer.calcErrorsGainBias(data, b(1:3), b(4:6), B_initial, G_initial);

            b_final = fminsearch(func, b_initial);

            obj.bias = b_final(1:3);
            obj.gain = b_final(4:6);

            data_cal = (data-obj.bias)./obj.gain;
            
            S = sum( sum(data_cal.^2,2) - 1, 1);

            if obj.debug_bit
                
                fprintf('Calibration complete. GAIN= %f %f %f | BIAS= %f %f %f\n', ...
                    obj.gain(1), obj.gain(2), obj.gain(3), obj.bias(1), obj.bias(2), obj.bias(3));
                
                fprintf('---> Total summed errors= %f\n', S);
                
            end
            
            if nargout>0
                sum_err_square = S;
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function showCalibration(obj, varargin)
            
            X = obj.calibration_data(:,1);
            Y = obj.calibration_data(:,2);
            Z = obj.calibration_data(:,3);
            
            scatter3(X,Y,Z,'.b');
            
            if ~isempty(obj.gain) && ~isempty(obj.bias) % also show the calibrated results
                hold on;
                scatter3(X./obj.gain(1)+obj.bias(1), Y./obj.gain(2)+obj.bias(2), Z./obj.gain(3)+obj.bias(3), '.r');
                hold off;
            end
            
        end
        
        function liveView(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('N', 1000, 'number', 'measurements');
            input.input_var('interval', 0.1);
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            obj.reset;
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            for ii = 1:input.N
                
                obj.update; 
                
                pause(input.interval);
                
                color = linspace(0,1,size(obj.Data.St,1));                
                color = [color', zeros(size(color,2),2)];
                
                scatter3(input.ax, obj.Data.St(:,1), obj.Data.St(:,2), obj.Data.St(:,3), 3, color);
                
                input.ax.XLim = [-1 1]*1.2;
                input.ax.YLim = [-1 1]*1.2;
                input.ax.ZLim = [-1 1]*1.2;
                
                title(input.ax, num2str(obj.acc_vec));
                
            end
            
        end
        
    end    
   
    methods (Static=true)
       
        function sum_err_square = calcErrorsGainBias(data, Biases, Gains, B_initial, G_initial) % calculate the calibration fit

            new_data = (data-Biases)./Gains;

            errors = sqrt(sum((new_data).^2,2)) - 1;

            bayes_term = sum(cosh(Biases./B_initial).*cosh(Gains./G_initial)); % this keeps the parameters from getting too big

            errors = errors.*bayes_term.*1e-5;

            sum_err_square = sum(errors.^2,1);

        end 
        
    end
    
end

