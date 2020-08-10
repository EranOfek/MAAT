%--------------------------------------------------------------------------
% AstFilter class                                                    class
% Description: A class of structure array of astronomical transmission
%              filters.
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef AstFilter 
    properties (SetAccess = public)
        family
        band
        T
        Tunits
        nT
        min_wl
        max_wl
        eff_wl
        pivot_wl
        pivot_wl_photon
        half_width
        fwhm
        comments
        source
        UserData
    end
    
    % class constructor method
    methods
        function AF=AstFilter(N,M)
            % Description: AstFilter constructor method
            
            Field = 'family';
            
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
                    AF(I,J).(Field) = [];
                end
            end
        
            
        end
    end
    
    % Static functions
    methods (Static)
        
        function Ans=isAstFilter(Obj)
            % Return true if object is AstFilter
            % Description: Check if object is of AstFilter class.
            % Input  : - Object
            % Output : - {true|false}.
            %     By : Eran O. Ofek                    Oct 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: S=1; AstFilter.isAstFilter(S);
            % Reliable: 2
            Ans = isa(Obj,'AstFilter');
        end

        function AstF=get(varargin)
            % Get astronomical transmission filters by name
            % Package: @AstFilter
            % Description: Get astronomical filters in AstFilter class format using
            %              exact name search on family or band names.
            %              This function replaces get_filter.m.
            % Input  : - Family name or a cell array of family names.
            %            If empty then search only by band name.
            %            Note that this function have a different behavior if
            %            the first argument is of AstFilter class.
            %          - Band name or a cell array of band names.
            %            If empty then search only by family name.
            %          - Behaviour options:
            %            'i' - case insenstive search (default).
            %            'c' - case sensitive search.
            % Output : - Astronomical filter class structure array of selected
            %            filters.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Feb 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: AstF=AstFilter.get('SDSS','r');
            %          AstF=AstFilter.get([],'r');
            %          AstF=AstFilter.get;
            %          AstF=AstFilter.get({'PTF','sdss'});
            % Reliable: 2
            %--------------------------------------------------------------------------

            LoadCat = false;
            if (numel(varargin)>0)
                if (AstFilter.isAstFilter(varargin{1}))
                    AstFilterCat = varargin{1};
                    if (numel(varargin)>1)
                        varargin = varargin(2:end);
                    else
                        varargin = {};
                    end
                else
                    LoadCat = true;
                end
            else
                LoadCat = true;
            end

            if (LoadCat)
                %AstFilterCat = Util.IO.load_check('AstFilterCat.mat');
                AstFilterCat = cats.spec.filter.AstFilterCat;
            end

            Def.Family = [];
            Def.Band   = [];
            Def.Behav  = 'i'; 

            Narg = numel(varargin);
            if (Narg==0)
                Family = Def.Family;
                Band   = Def.Band;
                Behav  = Def.Behav;
            elseif (Narg==1)
                Family = varargin{1};
                Band   = Def.Band;
                Behav  = Def.Behav;
            elseif (Narg==2)
                Family = varargin{1};
                Band   = varargin{2};
                Behav  = Def.Behav;
            elseif (Narg==3)
                Family = varargin{1};
                Band   = varargin{2};
                Behav  = varargin{3};

            else
                error('Illegal number of input arguments');
            end

            switch lower(Behav)
                case 'i'
                    % case insensetive
                    BehavFun = @strcmpi;
                case 'c'
                    BehavFun = @strcmp;
                otherwise
                    error('Unknown Behav option');

            end

            Nf           = numel(AstFilterCat);

            if (isempty(Family) && isempty(Band))
                % get all filters
                IndF = (1:1:Nf);    
            else
                if (isempty(Band))
                    % Family is empty
                    % Band is given
                    if (iscell(Family))
                        Nfam = numel(Family);
                        IndF = [];
                        for Ifam=1:1:Nfam
                            IndF = [IndF, find(BehavFun({AstFilterCat.family},Family{Ifam}))];
                        end
                    else
                        IndF = find(BehavFun({AstFilterCat.family},Family));
                    end
                elseif (isempty(Family))
                    % Band is empty 
                    % Family is given
                    if (iscell(Band))
                        Nband = numel(Band);
                        IndF = [];
                        for Iband=1:1:Nband
                            IndF = [IndF, find(BehavFun({AstFilterCat.band},Band{Iband}))];
                        end
                    else
                        IndF = find(BehavFun({AstFilterCat.band},Band));
                    end
                else
                    % both band and family are provided
                    IndF = find(BehavFun({AstFilterCat.band},Band) & ...
                                BehavFun({AstFilterCat.family},Family));
                end
            end
            AstF = AstFilterCat(IndF);

        end    % AstFilter.get function
        
    end  % Static
    
    % queries
    methods
        
        % get all family names
        function Family=all_family(AstF)
            % Get all family names in an AstFilter object
            % Input  : - An AstFilter object
            % Output : - A cell array of family names
            % Example: all_family(AstFilter.get)
            % Reliable: 2
            
            Family=unique({AstF.family});
            
        end
            
        
        
    end
    
    methods

      
        
        function Obj=search(AstF,Family)
            % Search AstFilter class by family name substring
            % Package: @AstFilter
            % Description: Search AstFilter class by family name substring.
            % Input  : - AstFilter class.
            %          - Family substring to search.
            % Output : - AstFilter class containing all possible matches.
            
            Obj = AstF(~Util.cell.isempty_cell(strfind(lower({AstF.family}),lower(Family))));
            
        end
        
    end
    
    
    % construct filters
    methods (Static)
        
        function AF=top_hat(L1,L2)
            % construct a top hat filter between L1 and L2
            % Input  : - Lower wavelngth (L1; scalar or a vector). If
            %            vector than output is an array of AstFilter
            %            objects.
            %          - L2
            % Output : - An AstFilter object
            
            
            Nf = numel(L1);
            AF = AstFilter(Nf,1);
            for If=1:1:Nf
                W = [L1(If)-0.01, linspace(L1(If),L2(If),10), L2(If)+0.01]';
                T = zeros(size(W));
                T(2:end-1) = 1;
                AF(If).family = 'top_hat';
                AF(If).band   = sprintf('f%05d',If);
                AF(If).T      = [W,T];
                AF(If).nT     = AF(If).T;
                
            end
            AF        = AF.norm;
            AF        = AF.pop_wl;
        end
        
        
    end
    
    
    %----------------------
    %--- plot functions ---
    %----------------------
    methods
            
        function H=plot(AstF,varargin)
            % Plot the transmission of an AstFilter object
            % Package: @AstFilter
            % Description: Given an AstFilter class, plot the normalized
            %              transmissions of all the filters.
            % Input  : - AstFilter class.
            %          * Additional arguments to pass to plot.m
            % Output : - Vector of handles.
            % Example: F.plot
            
            Nf = numel(AstF);
            H  = zeros(Nf,1);
            for If=1:1:Nf
                H(If) = plot(AstF(If).nT(:,1),AstF(If).nT(:,2),varargin{:});
                hold on;
            end
            hold off;
            H = xlabel('Wavelength [\AA]');
            H.Interpreter = 'latex';
            H.FontSize    = 16;
            H = ylabel('Normalized transmission');
            H.Interpreter = 'latex';
            H.FontSize    = 16;
        end
            
        function H=stairs(AstF,varargin)
            % Plot (stairs) an AstFilter object
            % Description: Given an AstFilter class, plot the normalized
            %              transmissions of all the filters, using the
            %              stairs function.
            % Input  : - AstFilter class.
            %          * Additional arguments to pass to plot.m
            % Output : - Vector of handles.
            % Example: F.plot
            
            Nf = numel(AstF);
            H  = zeros(Nf,1);
            for If=1:1:Nf
                H(If) = plot(AstF(If).nT(:,1),AstF(If).nT(:,2),varargin{:});
                hold on;
            end
            hold off;
            H = xlabel('Wavelength [\AA]');
            H.Interpreter = 'latex';
            H.FontSize    = 16;
            H = ylabel('Normalized transmission');
            H.Interpreter = 'latex';
            H.FontSize    = 16;
        end
        
    end
    
    
    % write filters to file
    methods
        function write_filter(F,FileName)
            % Write normalized filter transmission curve to file
            % Description: Write normalized filter transmission curve to file
            % Input  : - A single element AstFilter object
            %          - File name to write
            % Output : null
            
            FID = fopen(FileName,'w');
            fprintf(FID,'%9.2f  %10.8f\n',F.nT');
            fclose(FID);
            
        end
    end
    
    %------------------
    %--- Statistics ---
    %------------------
    methods
        function [IntT,IntnT]=integral(AstF)
            % The integrals of an AstFilter object
            % Description: Calculate the integral of the .T and .nT fields
            %              in an AstFilter class.
            % Input  : - AstFilter class.
            % Output : - A vector of integrals of the .T field in each of
            %            the AstFilter class elements.
            %          - A vector of integrals of the .nT field in each of
            %            the AstFilter class elements.
            % Example: [IntT,IntnT]=AstF.integral;
            
            Nf = numel(AstF);
            IntT  = zeros(Nf,1);
            IntnT = zeros(Nf,1);
            for If=1:1:Nf
                if (isempty(AstF(If).T))
                    IntT(If)  = NaN;
                    IntnT(If) = NaN;
                else
                    IntT(If)  = trapz(AstF(If).T(:,1), AstF(If).T(:,2));
                    if (nargout>1)
                        IntnT(If) = trapz(AstF(If).nT(:,1),AstF(If).nT(:,2));
                    end
                end
            end
            
        end
        
        function [MaxT,MaxnT]=max_tran(AstF)
            % The maximum transmissions in AstFilter object
            % Description: Calculate the max of the .T and .nT fields
            %              in an AstFilter class.
            % Input  : - AstFilter class.
            % Output : - A vector of max values of the transmission in
            %            the .T field.
            %          - A vector of max values of the transmission in
            %            the .nT field.
            % Example: [MaxT,MaxnT]=max_tran(AstF);
            
            Nf = numel(AstF);
            MaxT  = zeros(Nf,1);
            MaxnT = zeros(Nf,1);
            for If=1:1:Nf
                MaxT(If)  = max(AstF(If).T(:,2));
                if (nargout>1)
                    MaxnT(If) = max(AstF(If).nT(:,2));
                end
            end
        end
        
        function [MaxWave,MaxT,MaxnT]=max_wave(AstF)
            % Wavelengths of maximum transmission in AstFilter object
            % Description: Get the wavelength of maximum transmission in
            %              of filter.
            % Input  : - AstFilter class.
            % Output : - A vector of wavelengths of maximum transmission.
            %          - Maximum transmission in .T
            %          - Maximum transmission in .nT.
            % Example: MaxWave = max_wave(AstF);
            
            Nf = numel(AstF);
            MaxWave = zeros(Nf,1);
            MaxT    = zeros(Nf,1);
            MaxnT   = zeros(Nf,1);
            for If=1:1:Nf
                [~,MaxInd] = max(AstF(If).nT(:,2));
                MaxWave(If)    = AstF(If).nT(MaxInd,1);
                if (nargout>1)
                    MaxT(If) = AstF(If).T(MaxInd,2);
                    if (nargout>2)
                        MaxnT(If) = AstF(If).nT(MaxInd,2);
                    end
                end
            end
            
        end
        
      
        function AstF=norm(AstF,Norm)
            % Normalize the integrals of the transmissions of AstFilter object
            % Description: normalize the integral of the .nT field to be Norm.
            % Input  : - AstFilter class.
            %          - Optional normalization. Default is 1.
            % Output : - AstFilter class in which the .nT field integral
            %            is normalized to Norm.
                
            if (nargin==1)
                Norm = 1;
            end
            Nf = numel(AstF);
            for If=1:1:Nf
                Integral = trapz(AstF(If).nT(:,1),AstF(If).nT(:,2));
                AstF(If).nT(:,2) = AstF(If).nT(:,2).*Norm./Integral;
            end
        end
             
    end
    
    % Shift operations
    methods
        function ShiftedAstF=shift(AstF,z)
            % Description: Shift the wavelength in AstFilter class by
            %              a multiplicative (1+z) redshift/blueshift factor.
            % Input  : - AstFilter class.
            %          - A blue shift (positive) or redshift (negative).
            % Output : - AstFilter class with shifted wavelength in the
            %            .T and .nT fields, and updated parameters.
            % Example: AstF.shift(0.1)
            Nf = numel(AstF);
            ShiftedAstF = AstFilter(Nf,1);
            for If=1:1:Nf
                ShiftedAstF(If).T  = [AstF(If).T(:,1).*(1+z),  AstF(If).T(:,2)];
                ShiftedAstF(If).nT = [AstF(If).nT(:,1).*(1+z), AstF(If).nT(:,2)];
            end
            ShiftedAstF = pop_wl(ShiftedAstF);
            
        end
        
    end
        
        
    %----------------------------------
    %--- Interpolation and sampling ---
    %----------------------------------
    methods
        function AstF=interp(AstF,W,varargin)
            % Description: Interpolate an astronomical filter class
            %               transmission cuve into a new wavelngth frid.
            % Input  : - AstFilter class.
            %          - Column vector of wavelength grid [Ang].
            % Output : - AstFilter class with a new wavelength grid in
            %            the .T and .nT fields.
            
            Nf = numel(AstF);
            for If=1:1:Nf
                AstF(If).nT = [W, interp1(AstF(If).nT(:,1),AstF(If).nT(:,2),W,varargin{:})];
                AstF(If). T = [W, interp1(AstF(If).T(:,1), AstF(If).T(:,2), W,varargin{:})];
            end
        end
            
        
        function AstF=equalize_sampling(AstF1,AstF2)
            % Description: Equalize the transmission sampling grid of
            %              two astronomicakl filter class objects.
            %              Change the sampling of the first object to
            %              be equal to that of the second object.
            % Input  : - AstFilter class 
            %          - AstFilter class
            % Output : - The first AstFilter class resampled.
                
            N1 = numel(AstF1);
            N2 = numel(AstF2);
            N  = max(N1,N2);
            AstF = AstFilter(N,1);
            for I=1:1:N
                If1 = min(I,N1);
                If2 = min(I,N2);
                
                AstF(I) = interp(AstF1(If1),AstF2(If2).nT(:,1));
                
            end
        end
        
        function [AstF1,AstF2]=common_sampling(AstF1,AstF2)
            % Description: Resample two AstFilter objects such they will
            %              have a grid which is common to two objects.
            %              I.e., the new grid will contain the total range
            %              of the two grids.
            % Input  : - AstFilt object.
            %          - AstFilt object.
            % Output : - The first AstFilt object resampled to a common
            %            grid.
            %          - The second AstFilt object resampled to a common
            %            grid.
            
            N1 = numel(AstF1);
            N2 = numel(AstF2);
            N  = max(N1,N2);
            AstF = AstFilter(N,1);
            for I=1:1:N
                If1 = min(I,N1);
                If2 = min(I,N2);
                
                % choose the minimal sampling rate
                Samp1 = min(diff(AstF1(If1).nT(:,1)));
                Samp2 = min(diff(AstF2(If2).nT(:,1)));
                Samp  = min(Samp1,Samp2);
       
                % choose the maximal range
                Min   = min([AstF1(If1).nT(:,1);AstF2(If2).nT(:,1)]);
                Max   = max([AstF1(If1).nT(:,1);AstF2(If2).nT(:,1)]);
                
                % common sampling grid
                Wave  = (Min:Samp:Max)';
                
                % interpolate to new grid
                AstF1(I) = interp(AstF1(If1),Wave);
                if (nargout>1)
                    AstF2(I) = interp(AstF2(If2),Wave);
                end
            end
        end
        
    end
        
                
           
    % populating the class
    methods
            
        % add Filter
        function OldF=add_filter(NewF,OldF)
            % Description: Given an AstFilter object with minimal fields
            %              populated (i.e., family, band, T), populate
            %              its additional fields and concat it to another
            %              AstFilter object (default is the AstFilterCat).
            % Input  : - AstFilter class with minimal fields populated.
            %          - An AstFilter object to which to concat the first
            %            argument (after the other fields populaNew.pted).
            %            Default is to load AstFilterCat.mat.
            %            Note that this function does not save the results
            %            to AstFilterCat.mat (see the save method).
            %            If empty then will create an empty AstFilter
            %            object.
            % Output : - Populated AstFilter object.
            % Example: AstFilterCat=add_filter(NewF);
                            
            if (nargin==1)
                OldF = Util.IO.load_check('AstFilterCat.mat');
            end
            
            Del1 = false;
            if (~AstFilter.isAstFilter(OldF))
                if (isempty(OldF))
                    OldF = AstFilter;
                    Del1 = true;
                else
                    error('Old AstFilter to which to add is not of AstFilter class');
                end
            end
            % check validity of New AstFilter
            Nf = numel(NewF);
            for If=1:1:Nf
                if (isempty(NewF(If).family))
                    error('New AstFilter number %d have empty family name',If);
                end
                if (isempty(NewF(If).band))
                    error('New AstFilter number %d have empty band name',If);
                end
                if (size(NewF(If).T,2)<2)
                    error('New AstFilter number %d have illegal T field',If);
                end
                % make sure T doesn't have negative transmission
                if any(NewF(If).T(:,2)<0)
                    warning('New AstFilter number %d have negative transmission',If);
                    Ans = input('Do you want to set negative transmission to zero [Y/n]?','s');
                    switch lower(Ans)
                        case 'n'
                            % skip
                        otherwise
                            NewF(If).T(NewF(If).T(:,2)<0,2) = 0;
                    end
                end
                
                % check that T is sorted
                if (~issorted(NewF(If).T(:,1)))
                    error('New AstFilter number %d have unsorted transmission table');
                end
                
                % look for lines with duplicated wavelength
                DupInd = find(diff(NewF(If).T(:,1))==0);
                if (~isempty(DupInd))
                    warning('Lines with duplicated wavelength deleted in AstFilter number %d',If);
                    NewF(If).T = Util.array.delete_ind(NewF(If).T,DupInd,1);
                end
                
                
                % populate the nT field
                NewF(If).nT = NewF(If).T;
            end
            % noramlize the nT field
            NewF = norm(NewF);
            
            % populate additional fields
            NewF = pop_wl(NewF);
            
            % concat to old AstFilter
            OldF = [OldF(:); NewF(:)];
            
            if (Del1)
                OldF = OldF(2:end);
            end
            
        end
           
        function AstF=pop_wl(AstF)
            % Description: Populate the min_wl, max_wl, eff_wl, and
            %              half_width fields of AstFilter class
            %              based on the .nT transmission field.
            % Input  : - AstFilter class.
            % Output : - AstFilter class with fields repopulated.
            
            Nf = numel(AstF);
            for If=1:1:Nf
                % calculate min_wl
                AstF(If).min_wl = min(AstF(If).nT(AstF(If).nT(:,2)>0,1));
                
                % calculate max_wl
                AstF(If).max_wl = max(AstF(If).nT(AstF(If).nT(:,2)>0,1));
                
                % calculate eff_wl
                AstF(If).eff_wl = nansum(AstF(If).nT(:,1).*AstF(If).nT(:,2))./nansum(AstF(If).nT(:,2));

                % calclulate pivol_wl
                ind = ~isnan(AstF(If).nT(:,2));
                AstF(If).pivot_wl = sqrt(trapz(AstF(If).nT(ind,1),AstF(If).nT(ind,2))./...
                                         trapz(AstF(If).nT(ind,1),AstF(If).nT(ind,2)./AstF(If).nT(ind,1).^2));
                AstF(If).pivot_wl_photon = sqrt(trapz(AstF(If).nT(ind,1),AstF(If).nT(ind,2).*AstF(If).nT(ind,1))./...
                                                trapz(AstF(If).nT(ind,1),AstF(If).nT(ind,2)./AstF(If).nT(ind,1)));
                
                % calculate filter half_width
                % the width tranmit half the flux
                Fnn = AstF(If).nT(~isnan(AstF(If).nT(:,2)),:);
                
%                 CumSum = cumsum((Fnn(1:end-1,2)+eps).*diff(Fnn(:,1)));
                CumTrapz = cumtrapz(Fnn(:,1),Fnn(:,2));
                if abs(CumTrapz(end)-1)>(1e-5)
                    warning('Filter #%d, %s %s is not normalized',If, AstF(If).family, AstF(If).band)
                    CumTrapz = CumTrapz./CumTrapz(end);
                end
%                 for i=2:length(CumSum)
%                     if CumSum(i)<=CumSum(i-1)
%                         CumSum(i)=CumSum(i)+eps+(CumSum(i-1)-CumSum(i));
%                     end
%                 end
                for i=2:length(CumTrapz)
                    if CumTrapz(i)<=CumTrapz(i-1)
                        CumTrapz(i)=CumTrapz(i-1)+10.*eps;
                    end
                end
%                 AstF(If).half_width = interp1(CumSum+eps,Fnn(1:end-1,1),0.75) - ...
%                                       interp1(CumSum+eps,Fnn(1:end-1,1),0.25)
%                 CumTrapz = CumTrapz + 10000.*eps.*(1:1:numel(CumTrapz)).';
                AstF(If).half_width = interp1(CumTrapz./max(CumTrapz),Fnn(:,1),0.75) - ...
                                      interp1(CumTrapz./max(CumTrapz),Fnn(:,1),0.25);


                                  
                % calculate fwhm
                % the width at which the transmission drops by 1/2 relative
                % to max
                MaxnT = max(AstF(If).nT(:,2));
                I1 = find(AstF(If).nT(:,2)>(0.5.*MaxnT),1,'first');
                I2 = find(AstF(If).nT(:,2)>(0.5.*MaxnT),1,'last');
                AstF(If).fwhm = AstF(If).nT(I2,1) - AstF(If).nT(I1,1);
                
                
                              
            end
        end
        
        
        function AstF=save(AstF,FileName)
            % Description: Save an AstFilter class array into a MAT file.
            % Input  : - AstFilter class.
            %          - Full path and file name in which to save the
            %            AstFilter object.
            %            Default is '~/matlab/data/+cats/+spec/+filter/AstFilterCat.mat'
            % Outout : - AstFilter object saved.
            
            if (nargin<2)
                if isunix
                    FileName = sprintf('~%smatlab%sdata%s+cats%s+spec%s+filter%sAstFilterCat.mat',filesep,filesep,filesep,filesep,filesep,filesep);
                elseif ispc
                    data_dir = get_data_folder;
                    FileName = sprintf('%s%s+cats%s+spec%s+filter%sAstFilterCat.mat',data_dir,filesep,filesep,filesep,filesep,filesep);
                end
            end
            
            if (exist(FileName,'file')>0)
                FileNameOld = sprintf('%s.%s',FileName,date);
                fprintf('old file name %s exist\n',FileName);
                fprintf('old file name will be copied into %s\n',FileNameOld);
            end
            AstFilterCat = AstF;
            save(FileName,'AstFilterCat');
            
            
        end
        
    end
    
    %-----------------
    %--- Operators ---
    %-----------------
    methods      
        function AstF=filter_fun2(AstF1,AstF2,Fun,varargin)
            % Description: Binary function on AstFilter class objects
            % Input  : - AstFilter array
            %          - AstFilter array
            %          - Function handle
            %          * Additional arguments to pass to the function
            % Output : - Function output
            % Example: filter_fun2(AstF1,AstF2,@plus)
            
            N1 = numel(AstF1);
            N2 = numel(AstF2);
            N  = max(N1,N2);
            AstF = AstFilter(N,1);
            for I=1:1:N
                If1 = min(I,N1);
                If2 = min(I,N2);
                
                [AstF1(If1),AstF2(If2)] = common_sampling(AstF1(If1),AstF2(If2));
                
                % add
                AstF(I).nT = [AstF1(If1).nT(:,1), Fun(AstF1(If1).nT(:,2),AstF2(If2).nT(:,2),varargin{:})];
                AstF(I).T  = [AstF1(If1).T(:,1),  Fun(AstF1(If1).T(:,2),AstF2(If2).T(:,2),varargin{:})];
                
                % populate other fields
                AstF(I) = pop_wl(AstF(I));
            end
            
        end
            
        function AstF=filter_fun1(AstF1,Fun,varargin)
            % Description: Unary function on AstFilter class objects
            % Input  : - AstFilter array
            %          - Function handle
            %          * Additional arguments to pass to the function
            % Output : - Function output
            % Example: filter_fun1(AstF1,@sin)
            
            N1 = numel(AstF1);
            AstF = AstFilter(N1,1);
            for If1=1:1:N1
                % function
                AstF(If1).nT = [AstF1(If1).nT(:,1), Fun(AstF1(If1).nT(:,2),varargin{:})];
                AstF(If1).T  = [AstF1(If1).T(:,1),  Fun(AstF1(If1).T(:,2),varargin{:})];
                
                % populate other fields
                AstF(If1) = pop_wl(AstF(If1));
            end
            
        end
        
        
        function AstF=plus(AstF1,AstF2)
            % Description: Add AstFilter arrays (+)
            % Input   : - AstFilter array
            %           - AstFilter array
            % Outoput : - Sum of the .T and .nT fields of the AstFilter arrays
            % Example : A+A
            
            AstF=filter_fun2(AstF1,AstF2,@plus);
            
        end

        function AstF=minus(AstF1,AstF2)
            % Description: Subtract AstFilter arrays (+)
            % Input   : - AstFilter array
            %           - AstFilter array
            % Outoput : - Subtraction of the .T and .nT fields of the AstFilter arrays
            % Example : A-A
            
            AstF=filter_fun2(AstF1,AstF2,@minus);
        end
        
        function AstF=times(AstF1,AstF2)
            % Description: Multiply AstFilter arrays (.*)
            % Input   : - AstFilter array
            %           - AstFilter array
            % Outoput : - Multiplication of the .T and .nT fields of the AstFilter arrays
            % Example : A.*A
            
            AstF=filter_fun2(AstF1,AstF2,@times);
        end
        
        function AstF=rdivide(AstF1,AstF2)
            % Description: Multiply AstFilter arrays (./)
            % Input   : - AstFilter array
            %           - AstFilter array
            % Outoput : - Right division of the .T and .nT fields of the AstFilter arrays
            % Example : A./A
            
            AstF=filter_fun2(AstF1,AstF2,@rdivide);
        end
        
        function AstF=lt(AstF1,AstF2)
            % Description: little than operator on AstFilter arrays (<)
            % Input   : - AstFilter array
            %           - AstFilter array
            % Outoput : - Larger than operator on .T and .nT fields of the AstFilter arrays
            % Example : A<A
            
            AstF=filter_fun2(AstF1,AstF2,@lt);
        end
        
        function AstF=gt(AstF1,AstF2)
            % Description: greater than operator on AstFilter arrays (>)
            % Input   : - AstFilter array
            %           - AstFilter array
            % Outoput : - Larger than operator on .T and .nT fields of the AstFilter arrays
            % Example : A>A
            
            AstF=filter_fun2(AstF1,AstF2,@gt);
        end
        
        function AstF=le(AstF1,AstF2)
            % Description: lesser-equal than operator on AstFilter arrays (<=)
            % Input   : - AstFilter array
            %           - AstFilter array
            % Outoput : - Larger than operator on .T and .nT fields of the AstFilter arrays
            % Example : A<=A
            
            AstF=filter_fun2(AstF1,AstF2,@le);
        end
        
        function AstF=ge(AstF1,AstF2)
            % Description: greater-equal than operator on AstFilter arrays (>=)
            % Input   : - AstFilter array
            %           - AstFilter array
            % Outoput : - Larger than operator on .T and .nT fields of the AstFilter arrays
            % Example : A>=A
            
            AstF=filter_fun2(AstF1,AstF2,@ge);
        end
        
        function AstF=ne(AstF1,AstF2)
            % Description: not-equal than operator on AstFilter arrays (~=)
            % Input   : - AstFilter array
            %           - AstFilter array
            % Outoput : - Larger than operator on .T and .nT fields of the AstFilter arrays
            % Example : A~=A
            
            AstF=filter_fun2(AstF1,AstF2,@ge);
        end
        
        function AstF=eq(AstF1,AstF2)
            % Description: equal than operator on AstFilter arrays (==)
            % Input   : - AstFilter array
            %           - AstFilter array
            % Outoput : - Larger than operator on .T and .nT fields of the AstFilter arrays
            % Example : A==A
            
            AstF=filter_fun2(AstF1,AstF2,@ge);
        end
        
    end
    
    %--------------------------
    %--- Structre functions ---
    %--------------------------
    methods
        function obj=isfield(Sim,Field)
            % Description: isfield function for AstFilter
            % Input  : - AstFilter class object
            %          - Field name
            % Output : - Logical flag indicating if field exist.
            
            obj = any(strcmp(fieldnames(Sim),Field));
        end

        function obj=isstruct(Sim)
            % Description: isstruct function for AstFilter
            % Input  : - AstFilter class object
            % Output : - true.
            
            obj = true;  %isstruct(Sim) || isa(Sim,'SIM');
        end

    end
        
end

            
