% A class for Flux measurment objects
% Description: A class for flux measurments as a function of time
%              
% Tested : Matlab R2018b
%     By : Eran O. Ofek                    Jan 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 2
%--------------------------------------------------------------------------

classdef FluxClass
    properties (SetAccess = public)
        Time
        TimeUnits    % cell of str or str - default is JD
        Flux
        FluxErr      % 
        FluxUnits    % cell of str or str - default is AB mag
        Wave
        WaveUnits    % cell of str or str - default is Ang
        BandFamily   % cell of str or str, or ... 
        BandName     % cell of str or str
        Telescope
        Instrument
        Ref
    end
    
    properties (Hidden)
        
    end

    % constructor
    methods
        
        % constructor
        function S=FluxClass(varargin)
            % FluxClass constructor method
            % Package: @FluxClass
            % Input  : * Data,Keyword,value,...
            %            or keyword,value,... arguments.
            %            Data is a matrix with up to three coumns
            %            [Time, Flux, FluxErr], and keywords are any
            %            property of the FluxClass, while value is the
            %            value by which it will be populated.
            % Output : - A populated FluxClass object.
            
            if (nargin==0)
                Data = [];
            end
            
            Narg = numel(varargin);
            if (Narg./2)==floor(Narg./2)
                % assume key,val
                Data = [];
            else
                % assume Daya,key,val
                Data     = varargin{1};
                varargin = varargin(2:end);
                Narg     = Narg - 1;
            end
            
            for Iarg=1:2:Narg-1
                if ~isprop(S,varargin{Iarg})
                    error('Input keyword %s is not a property of FluxClass',varargin{Iarg});
                else
                    S.(varargin{Iarg}) = varargin{Iarg+1};
                end
            end
            
            Ncol = size(Data,2);
            DefColumns = {'Time','Flux','FluxErr'};
            if Ncol>numel(DefColumns)
                error('Number of input columns in first argument is larger than known default');
            end
            for Icol=1:1:Ncol
                S.(DefColumns{Icol}) = Data(:,Icol);
            end
            
        end

    end
    
    % get/set
    methods
        
        
%         function Array=get.St(S)
%             % get St property from stack (i.e., sorted stack)
%             
%             sorted_ind(S);
%             Array = S.Data(S.SI,:);
%         end
    end
    
    % conversions
    methods
        function F=convert(F,ColName,OutUnits,InUnits)
            % Convert the units of flux/time/wave values
            % Description: Conversion is done using convert.flux,
            %              convert.time and convert.energy, respectively.
            %              See these functions for options.
            % Input  : - A FluxClass object.
            %          - property/column name to convert (e.g., 'Flux').
            %          - Output units.
            %          - Input units. If not provided will attempt to use
            %            property units of object.
            % Output : - The FluxClass object with the new units.
            
            if (nargin<4)
                InUnits = [];
            end
                
            
            switch lower(ColName)
                case 'flux'
                    UnitsField = 'FluxUnits';
                    ValField   = 'Flux';
                    
                    % treat missing units
                    if (isempty(InUnits))
                        InUnits = F.(UnitsField);
                    end
                    if (isempty(InUnits))
                        if (isempty(F.(UnitsField)))
                            error('FluxUnits must be populated, or provided');
                        end
                    end
                    % convert values
                    if (isempty(F.Wave) || isempty(F.WaveUnits))
                        error('Wave and WaveUnits must be populated');
                    end
                    F.(ValField) = convert.flux(F.(ValField), InUnits, OutUnits,F.Wave,F.WaveUnits);

                case 'fluxerr'
                    UnitsField = 'FluxUnits';
                    ValField   = 'FluxErr';
                    
                     % treat missing units
                    if (isempty(InUnits))
                        InUnits = F.(UnitsField);
                    end
                    if (isempty(InUnits))
                        if (isempty(F.(UnitsField)))
                            error('FluxUnits must be populated, or provided');
                        end
                    end
                    % convert values
                    if (isempty(F.Wave) || isempty(F.WaveUnits))
                        error('Wave and WaveUnits must be populated');
                    end
                    F.(ValField) = convert.flux(F.(ValField), InUnits, OutUnits,F.Wave,F.WaveUnits);


                case 'time'
                    UnitsField = 'TimeUnits';
                    ValField   = 'Time';
                    
                     % treat missing units
                    if (isempty(InUnits))
                        InUnits = F.(UnitsField);
                    end
                    if (isempty(InUnits))
                        if (isempty(F.(UnitsField)))
                            error('FluxUnits must be populated, or provided');
                        end
                    end
                    % convert values
                    F.(ValField) = convert.time(F.(ValField), InUnits, OutUnits);

                case 'wave'
                    UnitsField = 'WaveUnits';
                    ValField   = 'Wave';
                    
                     % treat missing units
                    if (isempty(InUnits))
                        InUnits = F.(UnitsField);
                    end
                    if (isempty(InUnits))
                        if (isempty(F.(UnitsField)))
                            error('FluxUnits must be populated, or provided');
                        end
                    end
                    % convert values
                    F.(ValField) = convert.energy(InUnits, OutUnits, F.(ValField));

                    
                otherwise
                    error('Unknwon field name to convert');
            end
                  
            
            
        end
        
        
        function F=columnvec(F)
            % convert properties to column vectors
            
            Fields = fields(F);
            Nf     = numel(Fields);
            for If=1:1:Nf
                F.(Fields{If}) = F.(Fields{If})(:);
            end
            
        end
        
        function F=sort_time(F)
            % sort FluxClass object by time
            
            Fields = fields(F);
            Nf     = numel(Fields);
            Nobj   = numel(F);
            for Iobj=1:1:Nobj
                [~,SI] = sort(F(Iobj).Time);
                Nsi    = numel(SI);
                for If=1:1:Nf
                    if (numel(F(Iobj).(Fields{If}))==Nsi)
                        F(Iobj).(Fields{If}) = F(Iobj).(Fields{If})(SI);
                    end
                end
                
            end
            
        end
        
    end
    
    % fitting/interpolations
    methods
        function Fit=fit_pl(F,FlagT,DoConvert)
           
            
            if (nargin==2)
                DoConvert = false;
            end
            
            if (DoConvert)
                F = convert(F,'Time','JD');
                F = convert(F,'Wave','Hz');
                F = convert(F,'Flux','mJy');
                %F = convert(F,'FluxErr','mJy');
            end
            
            %FlagT = F.Time>=TimeRange(1) & F.Time<=TimeRange(2);
            
            
            % fit F = A*t^alpha * nu^beta
            %     log(F) = log(A) +alpha*log(t) + beta*log(nu)
            
            Fit.Nobs = sum(FlagT); % number of observations in time range
            % design matrix
            H    = [ones(Fit.Nobs,1), log10(F.Time(FlagT)), log10(F.Wave(FlagT))];
            Y    = log10(F.Flux(FlagT));
            ErrY = F.FluxErr(FlagT)./F.Flux(FlagT);
            [Fit.Par,Fit.ParErr] = lscov(H,Y,ErrY);
            Fit.Resid    = Y - H*Fit.Par;
            Fit.ResidStd = (Fit.Resid./ErrY).^2;
            Fit.Chi2     = sum(Fit.ResidStd);
            Fit.Npar     = size(H,2);
            Fit.Dof      = Fit.Nobs - Fit.Npar;
            Fit.RMS      = std(Fit.Resid);
            Fit.MeanTime = mean(F.Time(FlagT));
            
        end
       
        
    end
    
  
  
end
