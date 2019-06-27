%--------------------------------------------------------------------------
% ClassPSF class                                                     class
% Description: A class of structure (single elemnt array) of point spread
%              function (PSF).
%              PSF are stored either as a single image, a function
%              that return the PSF (in a specific image position),
%              or a cube pf PSF at differenr positions in the image.
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef ClassPSF
    properties (SetAccess = public)
        PSF       % either a matrix, cell array of matrices, or structure of .ImPSF, or a function handle [ImPSF,ErrPSF]=F(Par,X,Y)
        
    end
    properties (Hidden = true)
        ErrPSF    
        ParPSF = {};  % Default parameters for PSF function
        CooPSF = [];  % corrdinates for PSFs [X,Y] - empty for all image...

    end
    

    %-------------------
    %--- Constructor ---
    %-------------------
    methods
        
        function Psf=ClassPSF(N,M)
            % Description: ClassPSF class constructor.
            % Input  : - Number of rows, or [row, columns] Default is 1.
            %          - Number of columns. Default is 1.
            % Output : - A ClassPSF object of the requested size.

            PSF_Field = 'PSF';

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
                    Psf(I,J).(PSF_Field) = [];
                end
            end
            
        end
        
                
    end
    
    % static methods
    methods (Static)
        function Ans=isClassPSF(Obj)
            % Return true if object is ClassPSF
            % Description: Check if object is of ClassPSF class.
            % Input  : - Object
            % Output : - {true|false}.
            %     By : Eran O. Ofek                    Oct 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: S=1; ClassPSF.isClassPSF(S);
            % Reliable: 2
            Ans = isa(Obj,'ClassPSF');
        end
    end % static methods
    
    
    %------------------------------------
    %--- Get PSF from ClassPSF object ---
    %------------------------------------
    methods

        function Psf=clear_psf(Psf)
            % Description: Clear the PSF field in a ClassPSF object.
            % Input  : - A ClassPSF object.
            % Output : - A ClassPSF object.
            PSF_Field    = 'PSF';
            ErrPSF_Field = 'ErrPSF';

            Npsf = numel(Psf);
            for Ipsf=1:1:Npsf
                Psf(Ipsf).(PSF_Field)    = [];
                Psf(Ipsf).(ErrPSF_Field) = [];
            end
        end
        
        function Psf=setpsf(Psf,MatP)
            % Description: set the PSF
            % Input  : - A ClassPSF object
            %          - PSF to populate in the ClassPSF object.
            %            If a cell array of matrices, than each cell will
            %            populate a different element in the ClassPSF
            %            object.
            % Output : - A ClassPSF object with the populated PSF
           
            PSF_Field    = 'PSF';

            Npsf = numel(Psf);
            for Ipsf=1:1:Npsf
                if (iscell(MatP))
                    Psf(Ipsf).(PSF_Field) = MatP{Ipsf};
                else
                    Psf(Ipsf).(PSF_Field) = MatP;
                end
            end
                
            
        end

        function [P,Pe]=getpsf(Psf,Coo,Par)
            % Get the PSF from a single element ClassPSF object.
            % Description: Get the PSF from a single element ClassPSF
            %              object.
            % Input  : - A single element ClassPSF object.
            %          - A two elements vector of [X,Y] coordinates at
            %            which to return the PSF.
            %            If empty then attempt to ignore coordinates.
            %            Default is empty.
            %          - Cell array of parameters to pass to the PSF
            %            function handle. Default is {}.
            % Output : - A matrix of PSF.
            %          - A matrix of error in PSF estimate.
            % See also: getmpsf.m
            % Example: [P,Pe]=getpsf(Psf(1));
            
            
            PSF_Field    = 'PSF';
            ErrPSF_Field = 'ErrPSF';
            
            Def.Coo = [];
            Def.Par = {};
            if (nargin==1)
                Coo = Def.Coo;
                Par = Def.Par;
            elseif (nargin==2)
                Par = Def.Par;
            else
                % do nothing
            end
            
            if (numel(Psf)>1)
                error('Input must be a single element ClassPSF object');
            end
                
           
           if (isnumeric(Psf.(PSF_Field))),
               % PSF is numeric - a single PSF for the entire image
               P  = Psf.(PSF_Field);
               Pe = Psf.(ErrPSF_Field);
           elseif isa(Psf.(PSF_Field),'function_handle')
               % PSF is a function
               if (isempty(Coo)),
                   % Coo is not provided will attempt to ask for PSF
                   % without Coo
                   [P,Pe] = Psf.(PSF_Field)(Par{:});
               else
                   if (nume(Coo)~=2),
                       error('Coo must have two elements [X,Y]');
                   end
                   [P,Pe] = Psf.(PSF_Field)(Par{:},Coo(1),Coo(2));
               end
           elseif (iscell(Psf.(PSF_Field))),
               % A cell array of PSFs each one corresponds to a different
               % location.
               % The PSF locations are specified by Par parameter.
               
               error('Cell array PSF is not supported yet');
           elseif (isstruct(Psf.(PSF_Field))),
               error('Structure array PSF is not supported yet');
           else
               error('Unknown Psf.PSF type');
           end
        end
        
        function [P,Pe]=getmpsf(Psf,Coo,Par)
            % Get the PSF from a multiple elements ClassPSF object.
            % Description: Get the PSF from a multiple elements ClassPSF
            %              object.
            % Input  : - A ClassPSF object.
            %          - A two column matrix of [X,Y] coordinates at
            %            which to return the PSF.
            %            If ClassPSF object contain a single element and
            %            Coo multiple elements, then return the PSF for
            %            each location.
            %            If empty then attempt to ignore coordinates.
            %            Default is empty.
            %          - Cell array of parameters to pass to the PSF
            %            function handle. Default is {}.
            % Output : - A cell array of PSF matrices.
            %          - A cell array of PSF matrices of error in
            %            PSF estimate.
            % See also: getpsf.m
            % Example: [P,Pe]=getmpsf(Psf(1));
        
            Def.Coo = [];
            Def.Par = {};
            if (nargin==1)
                Coo = Def.Coo;
                Par = Def.Par;
            elseif (nargin==2)
                Par = Def.Par;
            else
                % do nothing
            end
            
            Npsf = numel(Psf);
            Ncoo = size(Coo,1);
            N    = max(Npsf,Ncoo);
            P    = cell(N,1);
            Pe   = cell(N,1);
            for I=1:1:N,
                Ipsf = min(I,Npsf);
                Icoo = min(I,Ncoo);
               
                if (isempty(Coo)),
                    Coo1 = [];
                else
                    Coo1 = Coo(Icoo,:);
                end
                [P{I},Pe{I}] = getpsf(Psf(Ipsf),Coo1,Par);
            end
            
        end
        
    end
    
    % display PSF
    methods
        function showallpsf(P)
            % Display PSFs in surface plot
            
            Np = numel(P);
            for Ip=1:1:Np
                % for each PSF
                R = input(sprintf('Enter to show next PSF (numer %d)',Ip),'s');
                clf;
                surface(P(Ip).PSF);
            end
        end
        
    end
    
    %------------------------------
    %--- Some useful PSF shapes ---
    %------------------------------
    methods
        function Psf=insert_psf(Psf,P,varargin)
            % Description: Insert a PSF into a ClassPSF object.
            % Input  : - A ClassPSF object.
            %          - A PSF to insert into the ClassPSF object.
            %            This can be one of the followings:
            %            A numeric matrix to insert to each one of the
            %            ClassPSF object elements;
            %            A cell array of matrices. Each cell element will
            %            be inserted to the corresponding ClassPSF object
            %            element; A function handle that return a PSF
            %            matrix.
            %          * Arbitrary number of arguments to pass to the
            %            function handle that generate the PSF
            % Output : - A ClassPSF object with the PSF field populated.
            % Example: P=insert_psf(ClassPSF,ones(3,3);
            %          P=insert_psf(CalssPSF,@Kernel2.gauss);
            %          P=insert_psf(CalssPSF,@Kernel2.aper,3);
            %          P=insert_psf(CalssPSF,@Kernel2.annulus);
            %          P=insert_psf(CalssPSF,@Kernel2.exp,2.5,19,19);
            
            PSF_Field = 'PSF';
            
            Np   = numel(P);
            Npsf = numel(Psf);
            for Ipsf=1:1:Npsf
                if (iscell(P))
                    % P is a cell array of matrices
                    Icell = min(Ipsf,Np);
                    Psf(Ipsf).(PSF_Field) = P{Icell};
                elseif (isnumeric(P))
                    Psf(Ipsf).(PSF_Field) = P;
                elseif (isa(P,'function_handle'))
                    Psf(Ipsf).(PSF_Field) = P(varargin{:});
                    
                else
                    error('Unknown PSF type to insert');
                end
                
            end
            
        
        end
        
    end
    
    
    %-------------------------------------
    %--- Statistics of ClassPSF object ---
    %-------------------------------------
    methods
        
        
         function [SumPsf]=sum_psf(Psf,varargin)
            % Description: Calculate the sum of each PSF in a ClassPSF
            %              object.
            % Input  : - A ClassPSF object.
            %          - A two column matrix of [X,Y] coordinates at
            %            which to return the PSF.
            %            If ClassPSF object contain a single element and
            %            Coo multiple elements, then return the PSF for
            %            each location.
            %            If empty then attempt to ignore coordinates.
            %            Default is empty.
            %          - Cell array of parameters to pass to the PSF
            %            function handle. Default is {}.
            % Output : - An array in which each element is the sum of each
            %            PSF.
            % See also: rms_psf.m
            % Example: [S]=sum_psf(Psf(1));


            Npsf   = numel(Psf);
            P      = getmpsf(Psf,varargin{:});
            
            SumPsf = zeros(size(Psf));
            for Ipsf=1:1:Npsf
                SumPsf(Ipsf) = sum(P{Ipsf}(:));
            end
            
         end
        
         function [RmsPsf]=rms_psf(Psf,varargin)
            % Description: Calculate the sqrt of sum of squares of each
            %              PSF in a ClassPSF object.
            % Input  : - A ClassPSF object.
            %          - A two column matrix of [X,Y] coordinates at
            %            which to return the PSF.
            %            If ClassPSF object contain a single element and
            %            Coo multiple elements, then return the PSF for
            %            each location.
            %            If empty then attempt to ignore coordinates.
            %            Default is empty.
            %          - Cell array of parameters to pass to the PSF
            %            function handle. Default is {}.
            % Output : - An array in which each element is the sqrt of sum
            %            of squares of each PSF.
            % See also: sum_psf.m
            % Example: [S]=rms_psf(Psf(1));


            Npsf   = numel(Psf);
            P      = getmpsf(Psf,varargin{:});
            
            RmsPsf = zeros(size(Psf));
            for Ipsf=1:1:Npsf
                RmsPsf(Ipsf) = sqrt(sum(P{Ipsf}(:).^2));
            end
            
         end
        
         function [X,Y,X2,Y2,XY] = moment_psf(Psf,MomentSigma,varargin)
                % Description: Calculate 1st and 2nd moments of PSFs.
                %              Assumes all the PSF have the same size.
                % Input  : - A ClassPSF object.
                %          - Sigma of Gaussian by which to weight the pixels.
                %            Default is 1.5.
                %          - A two column matrix of [X,Y] coordinates at
                %            which to return the PSF.
                %            If ClassPSF object contain a single element and
                %            Coo multiple elements, then return the PSF for
                %            each location.
                %            If empty then attempt to ignore coordinates.
                %            Default is empty.
                %          - Cell array of parameters to pass to the PSF
                %            function handle. Default is {}.
                % Output : - X first moment.
                %          - Y first moment.
                %          - X^2 second moment.
                %          - Y^2 second moment.
                %          - X*Y second moment.
             
                Def.MomentSigma = 1.5;
                if (nargin==1)
                    MomentSigma = Def.MomentSigma;
                end
                if (isempty(MomentSigma))
                    MomentSigma = Def.MomentSigma;
                end
                
                SizeP = size(Psf);
                Npsf  = numel(Psf);
                P     = getmpsf(Psf,varargin{:});
                % allocate output
                X     = zeros(SizeP);
                Y     = zeros(SizeP);
                if (nargout>1)
                    X2     = zeros(SizeP);
                    Y2     = zeros(SizeP);
                    XY     = zeros(SizeP);
                end
                for Ipsf=1:1:Npsf
                    % calc moments for each PSF:
                    Size = size(P{Ipsf});
                    if (nargout>2)
                        [Mom,Mom2] = ImUtil.Im.im_moments(P{Ipsf}, floor(Size(2).*0.5),...
                                                         floor(Size(1).*0.5),...
                                                         floor(min(Size).*0.5), MomentSigma);
                        X(Ipsf)  = Mom.X;
                        Y(Ipsf)  = Mom.Y;
                        X2(Ipsf) = Mom2.X2;
                        Y2(Ipsf) = Mom2.Y2;
                        XY(Ipsf) = Mom2.XY;
                    else
                        [Mom]      = ImUtil.Im.im_moments(P{Ipsf}, floor(Size(2).*0.5),...
                                                         floor(Size(1).*0.5),...
                                                         floor(min(Size).*0.5), MomentSigma);
                        X(Ipsf)  = Mom.X;
                        Y(Ipsf)  = Mom.Y;
                    end
                end
                
                
                
         end
         
         function [CurveGrowth,RadHalf]=curve_growth_psf(Psf,varargin)
                % Description: Calculate the curve of growth for the PSFs
                %              in ClassPSF object.
                % Input  : - A ClassPSF object.
                %          - A two column matrix of [X,Y] coordinates at
                %            which to return the PSF.
                %            If ClassPSF object contain a single element and
                %            Coo multiple elements, then return the PSF for
                %            each location.
                %            If empty then attempt to ignore coordinates.
                %            Default is empty.
                %          - Cell array of parameters to pass to the PSF
                %            function handle. Default is {}.
                % Output : - A matrix of curve of growth.
                %            Each column represent a PSF. The first line is
                %            for the inner radius (1) and the last line is
                %            for the outer radius. Radii steps is 1.
                %            The curve of growth is normalize such that the
                %            sum of each PSF is 1.
                %          - Vector of radii indicating for each PSF, the
                %            radius that contains half of the light.
                
                RadStep = 1;
                
                Npsf  = numel(Psf);
                P     = getmpsf(Psf,varargin{:});
                % get the size of the 1st PSF
                % assumes all PSFs have the same size
                Size     = size(P{1});
                MaxRad   = floor(min(Size).*0.5);
                RadVec   = (1:RadStep:MaxRad);
                % allocate output
                CurveGrowth = zeros(MaxRad,Npsf);
                RadHalf     = zeros(Npsf,1);
                % calculate moments
                [X,Y] = moment_psf(Psf,[],varargin{:});
                [MatX,MatY] = meshgrid((1:1:Size(2)),(1:1:Size(1)));
                
                for Ipsf=1:1:Npsf
                    MatR = sqrt((MatX-X(Ipsf)).^2 + (MatY-Y(Ipsf)).^2);
                    % make sure PSF is normalized to unity
                    P{Ipsf} = P{Ipsf}./sum(P{Ipsf}(:));
                    for Rad=1:RadStep:MaxRad
                        CurveGrowth(Rad,Ipsf) = Util.stat.sumnd(P{Ipsf}(MatR<=Rad));
                    end
                    
                    Epsilon       = (1:1:numel(CurveGrowth)).*eps*1000;
                    RadHalf(Ipsf) = interp1(CurveGrowth.'+Epsilon,RadVec,0.5);
                end
                
         end
         
         % function fit_psf
        
    end
    
    
   
    %------------------------------
    %--- Pad, Shift and FFT PSF ---
    %------------------------------
    methods
          
         function S=psf_pad2mat(P,ImageSize,FftShift,varargin)
            % Description: shift and pad PSF and copy it to a matrix.
            % Input  : - ClassPSF object
            %          - Image size to pad the PSF.
            %            This is a two column matrix of [Y X] size.
            %          - A flag indicating if to fftshift the PSF.
            %            Default is false.
            % Output : - Matrix of PSF
            
            if (nargin<3)
                FftShift = false;
            end
            
            Nsize = size(ImageSize,1);
            
            ImageField   = SIM.ImageField;
            PSF_Field    = 'PSF';
            ErrPSF_Field = 'ErrPSF';
            
            PadVal = 0;
            [CellP,~] = getmpsf(P,varargin{:});
            
            Npsf = numel(P);
            if (Npsf>1)
                error('psf_pad2mat can work on a single PSF');
            end
            
            for Ipsf=1:1:Npsf
                % for each PSF
                % shift PSF
                
                Size = ImageSize(min(Nsize,Ipsf),:);
                % Pad the matrix format PSF by zeros
                SizeP = size(CellP{Ipsf});
                % number of pixels in each side
                PadHalfX = (Size(2)-SizeP(2)).*0.5;
                PadHalfY = (Size(1)-SizeP(1)).*0.5;

                % If need to pad by a non-integer amount (0.5)
                % than first shift the PSF by 0.5 pix and than pad
                % by a whole number of pixels.
                ShiftX = PadHalfX - floor(PadHalfX);
                ShiftY = PadHalfY - floor(PadHalfY);
                if (ShiftX==0 && ShiftY==0),
                    % no need to shift - do nothing
                    ShiftedP  = CellP{Ipsf};
                else
                    ShiftedP = ImUtil.Im.image_shift_fft(CellP{Ipsf},ShiftX,ShiftY);
                    % if PSF shifted in opne of the axes need to remove
                    % a corresponding column/row:
                    Xstart   = 1 + ceil(ShiftX);
                    Ystart   = 1 + ceil(ShiftY);
                    ShiftedP = ShiftedP(Ystart:end,Xstart:end);
                    % Fix the pad value accordingly
                    PadHalfX = ceil(PadHalfX);
                    PadHalfY = ceil(PadHalfY);
                end
                
                
                % pad PSF
                S    = padarray(ShiftedP,[PadHalfY,PadHalfX],PadVal,'both');

                if (FftShift),
                    S    = fftshift(S);
                end
                                    
                
            end
            
         end
        
         
        function S=psf_pad2sim(P,ImageSize,FftShift,varargin)
            % Description: shift and pad PSF and copy it to a SIM image
            % Input  : - ClassPSF object
            %          - Image size to pad the PSF.
            %            This is a two column matrix of [Y X] size.
            %          - A flag indicating if to fftshift the PSF.
            %            Default is false.
            % Output : - SIM array with the PSF at the image field.
            
            if (nargin<3)
                FftShift = false;
            end
            
            
            S = SIM(size(P));
            S = P;
            
            Nsize = size(ImageSize,1);
            
            ImageField   = SIM.ImageField;
            PSF_Field    = 'PSF';
            ErrPSF_Field = 'ErrPSF';
            
            PadVal = 0;
            [CellP,~] = getmpsf(P,varargin{:});
            
            Npsf = numel(P);
            for Ipsf=1:1:Npsf,
                % for each PSF
                % shift PSF
                
                Size = ImageSize(min(Nsize,Ipsf),:);
                % Pad the matrix format PSF by zeros
                SizeP = size(CellP{Ipsf});
                % number of pixels in each side
                PadHalfX = (Size(2)-SizeP(2)).*0.5;
                PadHalfY = (Size(1)-SizeP(1)).*0.5;

                % If need to pad by a non-integer amount (0.5)
                % than first shift the PSF by 0.5 pix and than pad
                % by a whole number of pixels.
                ShiftX = PadHalfX - floor(PadHalfX);
                ShiftY = PadHalfY - floor(PadHalfY);
                if (ShiftX==0 && ShiftY==0),
                    % no need to shift - do nothing
                    ShiftedP  = CellP{Ipsf};
                else
                    ShiftedP = ImUtil.Im.image_shift_fft(CellP{Ipsf},ShiftX,ShiftY);
                    % if PSF shifted in opne of the axes need to remove
                    % a corresponding column/row:
                    Xstart   = 1 + ceil(ShiftX);
                    Ystart   = 1 + ceil(ShiftY);
                    ShiftedP = ShiftedP(Ystart:end,Xstart:end);
                    % Fix the pad value accordingly
                    PadHalfX = ceil(PadHalfX);
                    PadHalfY = ceil(PadHalfY);
                end
                
                
                % pad PSF
                S(Ipsf).(ImageField)    = padarray(ShiftedP,[PadHalfY,PadHalfX],PadVal,'both');

                if (FftShift),
                    S(Ipsf).(ImageField)    = fftshift(S(Ipsf).(ImageField));
                end
                                    
                
            end
            
        end
            
        
        
        
        function P=psf_shift_fft(P,ShiftX,ShiftY,varargin)
            % Description: Shift (using the fft shift theorem) the PSF.
            %              Note that this may be useful for sub pixel
            %              shifts of the PSF. Only the PSF field will be
            %              shifted, and the shifted PSF will be stored in
            %              an array format.
            % Input  : - A ClassPSF object.
            %          - ShiftX [pix]. A scalar or vector of the length of
            %            ClassPSF. Note that the shift is applied to the
            %            source position.
            %          - ShiftY [pix]. A scalar or vector of the length of
            %            ClassPSF.
            %          * Additional arguments to pass to getmpsf.m.
            % Output : - A ClassPSF object in which the PSF field is the
            %            shifted PSF in an array format.
            % See also: ImUtil.Im.image_shift_fft.m
            % Example: P=psf_shift_fft(P,0.5,0.5);
            
            PSF_Field = 'PSF';
            
            CellP  = getmpsf(P,varargin{:});
            Npsf   = numel(CellP);
            ShiftX = ShiftX(:).*ones(Npsf,1);
            ShiftY = ShiftY(:).*ones(Npsf,1);
            for Ipsf=1:1:Npsf
                P(Ipsf).(PSF_Field)  = ImUtil.Im.image_shift_fft(CellP{Ipsf},ShiftX(Ipsf),ShiftY(Ipsf));
            end
            
        end
        
        function P=psf_pad_fftshift(P,SizeIm,FftShift,varargin)
            % Description: Given a PSF stamp of size NxM, which center is
            %              located in the middle of a the central pixel,
            %              where the central pixel index is given by
            %              ceil([N M].*0.5+0.1), pad the PSF (to some requested
            %              size) and ifftshift it (such that the PSF center
            %              will be located at the origin).
            % Input  : - A ClassPSF object in which the PSF centers are
            %            located in the middle of the central pixel, where
            %            the central pixel index is given by
            %            ceil([N M].*0.5+0.1)
            %          - Final image size of the PSF. The PSF will be pad
            %            to this size. If input is SIM, and this parameter
            %            is empty or not given then set the image size to
            %            the SIM image size.
            %          - A flag indicating if to ifftshift the PSF such
            %            that its center will be at the origin.
            %            Default is true.
            %          * Additional arguments to pass to getmpsf.m
            % Output : - A ClassPSF object in whic the PSF and ErrPSF
            %            fields are padded and shifted.
            % Example: P=psf_pad_fftshift(P,size(Image));
            
                
            PSF_Field    = 'PSF';
            ErrPSF_Field = 'ErrPSF';
            
            if (nargin==1),
                SizeIm   = [];
                FftShift = true;
            elseif (nargin==2),
                FftShift = true;
            else
                % do nothing
            end
            if (isempty(SizeIm))
                if (SIM.issim(P))
                    SizeIm = imagesize(P);
                else
                    error('SizeIm parameter must provided or PSF should be a SIM class');
                end
            end
            
            
            
            PadVal = 0;
            
            [CellP,CellPe] = getmpsf(P,varargin{:});
            Npsf  = numel(CellP);
            Nsize = size(SizeIm,1);
            
            for Ipsf=1:1:Npsf
                SizePSF   = size(CellP{Ipsf});
                SizeFinal = SizeIm(min(Nsize,Ipsf),:);
                CenterPSF = ceil(SizePSF.*0.5+0.1);
                CenterIm  = ceil(SizeFinal.*0.5+0.1);
                
                PadBefore = CenterIm - CenterPSF;
                PadAfter  = (SizeFinal-CenterIm) - (SizePSF-CenterPSF);
                CellP{Ipsf} = padarray(CellP{Ipsf},PadBefore,PadVal,'pre');
                CellP{Ipsf} = padarray(CellP{Ipsf},PadAfter,PadVal,'post');
                if (~isempty(CellPe{Ipsf}))
                    CellPe{Ipsf} = padarray(CellPe{Ipsf},PadBefore,PadVal,'pre');
                    CellPe{Ipsf} = padarray(CellPe{Ipsf},PadAfter,PadVal,'post');
                end
                
                if (FftShift)
                    P(Ipsf).(PSF_Field) = ifftshift(CellP{Ipsf});
                    if (~isempty(CellPe{Ipsf}))
                        P(Ipsf).(ErrPSF_Field) = ifftshift(CellPe{Ipsf});
                    end
                else
                    P(Ipsf).(PSF_Field) = CellP{Ipsf};
                    if (~isempty(CellPe{Ipsf}))
                        P(Ipsf).(ErrPSF_Field) = CellP{Ipsf};
                    end
                end
            end
            
        end
        
        function P=pad_psf(P,Size,Shift,varargin)
            % Description: Pad PSFs with zeros or return PSF in matrix
            %              format. The function treating shifts of PSF by
            %              half a pixel in order to keep the PSF centered.
            %              Optionaly also shift the PSF using fftshift.
            % Input  : - A ClassPSF object.
            %          - Size [Y,X] of final output PSF. The output PSF
            %            will be always in matrix format.
            %            If empty, then the output size is the same as the
            %            input PSF, but always in matrix format.
            %          - A flag indicating if to shift the PSF using
            %            fftshift post padding. Default is false.
            %          - A two column matrix of [X,Y] coordinates at
            %            which to return the PSF.
            %            If ClassPSF object contain a single element and
            %            Coo multiple elements, then return the PSF for
            %            each location.
            %            If empty then attempt to ignore coordinates.
            %            Default is empty.
            %          - Cell array of parameters to pass to the PSF
            %            function handle. Default is {}.
            % Output : - A ClassPSF object in which the PSF and ErrPSF
            %            fields are pad with zeros to the requested size
            %            and are in matrix format.
            %            Note that while the PSF field is shifted by 0.5pix
            %            if necessery. The ErrPSF is not shifted and should
            %            be regarded as an estimate based on nearest pixel.
            % See also: shift_psf.m
            % Example: P1=pad_psf(P,[]);
            %          P1=pad_psf(P,[1024 1024]);
            
            
            error('do not use!')
            
            % Shift parameter default
            if (nargin==2),
                Shift = false;
            end
            if (isempty(Shift)),
                Shift = false;
            end
            
            PSF_Field    = 'PSF';
            ErrPSF_Field = 'ErrPSF';
            
            PadVal = 0;
            [CellP,CellPe] = getmpsf(P,varargin{:});
            
            Npsf = numel(P);
            for Ipsf=1:1:Npsf
                % for each PSF
                
                if (isempty(Size))
                    % if input Size is empty then return the PSF in matrix
                    % format
                    P(Ipsf).(PSF_Field)    = CellP{Ipsf};
                    P(Ipsf).(ErrPSF_Field) = CellPe{Ipsf};
                else
                
                    % Pad the matrix format PSF by zeros
                    SizeP = size(CellP{Ipsf});
                    % number of pixels in each side
                    PadHalfX = (Size(2)-SizeP(2)).*0.5;
                    PadHalfY = (Size(1)-SizeP(1)).*0.5;

                    % If need to pad by a non-integer amount (0.5)
                    % than first shift the PSF by 0.5 pix and than pad
                    % by a whole number of pixels.
                    ShiftX = PadHalfX - floor(PadHalfX);
                    ShiftY = PadHalfY - floor(PadHalfY);
                    if (ShiftX==0 && ShiftY==0)
                        % no need to shift - do nothing
                        ShiftedP  = CellP{Ipsf};
                        ShiftedPe = CellPe{Ipsf};
                    else
                        ShiftedP = ImUtil.Im.image_shift_fft(CellP{Ipsf},ShiftX,ShiftY);
                        % if PSF shifted in opne of the axes need to remove
                        % a corresponding column/row:
                        Xstart   = 1 + ceil(ShiftX);
                        Ystart   = 1 + ceil(ShiftY);
                        ShiftedP = ShiftedP(Ystart:end,Xstart:end);
                        % Fix the pad value accordingly
                        PadHalfX = ceil(PadHalfX);
                        PadHalfY = ceil(PadHalfY);

                        % treat the PSF error
                        % using interpolation
                        if (isempty(CellPe{Ipsf}))
                            ShiftedPe = [];
                        else
                            ShiftedPe = CellPe{Ipsf}(Ystart:end,Xstart:end);
                        end
                        
                    end

                    % pad the PSF and PSF error
                    P(Ipsf).(PSF_Field)    = padarray(ShiftedP,[PadHalfY,PadHalfX],PadVal,'both');
                    if (isempty(ShiftedPe))
                        P(Ipsf).(ErrPSF_Field) = [];
                    else
                        P(Ipsf).(ErrPSF_Field) = padarray(ShiftedPe,[PadHalfY,PadHalfX],PadVal,'both');
                    end
                    
                    % Shift PSF using fftshift
                    if (Shift)
                        P(Ipsf).(PSF_Field)    = ifftshift(P(Ipsf).(PSF_Field));
                        P(Ipsf).(ErrPSF_Field) = ifftshift(P(Ipsf).(ErrPSF_Field));
                    end
                end
            end
        
            
        end
        
        function P=shift_psf(P,varargin)
            % Description: Use fft shift to generate a shifted version of
            %              the ClassPSF object in matrix format.
            %              Its recomended that the shifted image will have
            %              even size. Note that if you pad the PSF you
            %              first need to pad it and only than shift it.
            %              If you need to pad and shift the PSF it is
            %              recomended to use: pad_psf.m.
            % Input  : - A ClassPSF object.
            %          - A two column matrix of [X,Y] coordinates at
            %            which to return the PSF.
            %            If ClassPSF object contain a single element and
            %            Coo multiple elements, then return the PSF for
            %            each location.
            %            If empty then attempt to ignore coordinates.
            %            Default is empty.
            %          - Cell array of parameters to pass to the PSF
            %            function handle. Default is {}.
            % Output : - A ClassPSF object in which the PSF is presented in
            %            matrix form and is shifted using fftshift.
            %            If the input PSF is centered then the output PSF
            %            is shifted to the origin.
            % See also: pad_psf.m
            % Example: P = shift_psf(P);
            
            
            PSF_Field    = 'PSF';
            ErrPSF_Field = 'ErrPSF';
            
            % get THE PSF in matrix format
            [CellP,CellPe] = getmpsf(P,varargin{:});
            
            Npsf = numel(P);
            for Ipsf=1:1:Npsf
                % for each PSF
                P(Ipsf).(PSF_Field)    = fftshift(CellP{Ipsf});
                P(Ipsf).(ErrPSF_Field) = fftshift(CellPe{Ipsf});
            end
            
            
        end
        
        function P=psf_fft2(P,varargin)
            % Description: run fft2 on PSF
            % Input  : - A ClassPSF object.
            %          * Additional arguments to pass to getmpsf.m
            % Output : - A ClassPSF object with ffted PSF
            
            PSF_Field    = 'PSF';
            ErrPSF_Field = 'ErrPSF';
            
            % get THE PSF in matrix format
            [CellP,CellPe] = getmpsf(P,varargin{:});
            
            Npsf = numel(P);
            for Ipsf=1:1:Npsf
                % for each PSF
                P(Ipsf).(PSF_Field)    = fft2(CellP{Ipsf});
            end
        end
        
    end

        % See list of overload functions
        % http://www.mathworks.com/help/matlab/matlab_oop/implementing-operators-for-your-class.html

        
    methods
        function Psf=psf_extrapolate(Psf,Method,varargin)
            % 
           
            if (nargin==1)
                Method = 'poly2';
            end
            
            [CellP,~] = getmpsf(Psf,varargin{:});

            Npsf = numel(Psf);
            for Ipsf=1:1:Npsf
                % for each PSF
                
                switch lower(Method)
                    case 'poly2'
                        SizePSF = size(CellP{Ipsf});
                        Npix    = numel(CellP{Ipsf});
                        [MatX,MatY] = meshgrid((1:1:SizePSF(2)),(1:1:SizePSF(1)));
                        P           = CellP{Ipsf};

                        H   = [ones(Npix,1), MatX(:), MatY(:), MatX(:).^2, MatY(:).^2, MatX(:).*MatY(:)];
                        Ind0 = P>0;
                        Par = H(Ind0,:)\log10(P(Ind0));
                        Ppred = 10.^(H*Par);
                        P(P==0) = Ppred(P==0);
                        Psf(Ipsf).PSF = P./sum(P(:));
                    otherwise
                        error('unknown Method option');
                end
            
                
                
            end
            
            
        end
  
    end
        
end

            
