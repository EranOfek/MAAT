%--------------------------------------------------------------------------
% AstTran class                                                      class
% Description: An astronomical transformation class.
%              Define a 2-D transformation for astronomical images.
%              This is a structure array in which each element refers to
%              additive part of the transformation (e.g., the rotation).
%              Fields are:
%              'Type' - A char array of the transformation name
%                       'x_shift', 'y_shift'
%                       'x_scale', 'y_scale'
%                       'x_orirntation', 'y_orientation'
%                       'x_rot', 'y_rot'
%                       'xy_rot'
%                       'x_pixphase', 'y_pixphase'
%                       'x_platetilt', 'y_platetilt'
%              'Fun' - Function
%              'Par' - Parameters of transformation
%              'UserData' - General user data.
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Dec 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef TranClass
    properties (SetAccess = public)
        Fun     % A cell array (element per dimension) of cell array function handles
                % 
                % e.g., {{@FunOne, @FunX, @FunY},{@FunOne}} represent a transformation of
                % the form: X' = alpha + beta*X + gamma*Y
                %           Y' = delta
        FunArg  % A cell (element per dimension) of cell array of function arguments
        NPar    % A cell arry (element per dimension) of
                % optional vector in which each element is the number of free
                % parameters corresponding to the function in Fun
        Par     % A cell array (element per dimension) of cell array of vectors.
                % Each vector is the free parameters
                % corresponding to the function in Fun
        ParErr  % Error in the parameters (like Par)
        %Poly    % Symbolic polinomial representation
        %ChPoly  % Seprated symbolic polinomial representation
        
    end
    
    % constructor
    methods (Static)
        % TranClass constructor method
        function TranC=TranClass(varargin)
            % TranClass object constructor
            % Description: Constructor method for a TranClass object.
            % Input  : * Arbitarary number of elements. Element per
            %            dimension.
            %            Each lement is a pairs of
            %            ...,function, fun-argument,...
            % Outout : - A TranClass object.
            % Example: TranC=TranClass({@FunOne,[],@FunPolyXY,[1 2;2 3]}, {@FunOne,[],@FunPolyXY,[1 3; 1 1]})
            FunField     = TranClass.FunField;
            FunArgField  = TranClass.FunArgField;
            FunNParField = TranClass.NParField;
            
            IC = 1;
            
            Ndim = numel(varargin);
            % Each input element represent a dimension
            % for each dimension
            for Idim=1:1:Ndim
                Argument = varargin{Idim};
                
                Narg = numel(Argument);
                K = 0;
                for I=1:2:Narg
                    K = K + 1;
                    TranC(IC).(FunField){Idim}{K}     = Argument{I};
                    TranC(IC).(FunArgField){Idim}{K}  = Argument{I+1};
                    [~,NPar] = TranClass.evalfun(Argument{I},Argument{I+1});
                    TranC(IC).(FunNParField){Idim}{K} = NPar; 
                end
                
            end
        end
        
    end
    
    % static methods for field names
    methods (Static)
        function [Field]=FunField
            % Return the 'Fun' field name
            % Package: @TranClass
            Field = 'Fun';
        end
        function [Field]=FunArgField
            % Return the 'FunArg' field name
            % Package: @TranClass
            Field = 'FunArg';
        end
        function [Field]=NParField
            % Return the 'NPar' field name
            % Package: @TranClass
            Field = 'NPar';
        end
        function [Field]=ParField
            % Return the 'Par' field name
            % Package: @TranClass
            Field = 'Par';
        end
        function [Field]=ParErrField
            % Return the 'ParErr' field name
            % Package: @TranClass
            Field = 'ParErr';
        end
        function [Field]=PolyField
            % Return the 'Poly' field name
            % Package: @TranClass
            Field = 'Poly';
        end
        function [Field]=ChPolyField
            % Return the 'ChPoly' field name
            % Package: @TranClass
            Field = 'ChPoly';
        end
    end
    
    
    
    
    
    % static methods for transformation functions
    methods (Static)
        
        function [Out,NPar,SymFunT]=evalfun(FunT,Arg,varargin)
            % A static method to evaluate the transformation functions
            % Description: A static method to evaluate the transformation
            %              functions in the TranClass object.
            % Input  : - A function handle or a function name (string).
            %          - Matrix of arguments to pass to the function
            %            The Transformation functions:
            %            FunPolyXY, FunPolyChebyshev1XY, and FunPolyChebyshev2XY
            %            have the following arguments:
            %            A two rows matrix in which each
            %            column provide the orders of the X and Y
            %            polynomials.
            %            E.g., [1 1 2 2; 1 2 1 3] corresponds to:
            %            X.*Y, X.*Y.^2, X.^2.*Y, X.^2.*Y.^3
            %          * Arbitrary number of pairs of arguments:
            %            ...,keyword,value,...
            %            where keyword are parameters that needed for the
            %            evaluation of the function - obvious option may
            %            include:
            %            'X'    - Vector of X coordinates
            %            'Y'    - Vector of Y coordinates
            %            'Mag'  - Vector of magnitudes
            %            'Color'- Vector of colors
            % Output : - Vector or matrix of evaluated columns
            %          - Number of free parameters in the function
            % Example: [Out,NPar]=TranClass.evalfun(@FunOne,[],'X',[1;1]);
            %          [Out,NPar]=TranClass.evalfun(@FunPolyXY,[1 2 3;2 2 2],'X',[1;3],'Y',[1;2]);
            
            if (isa(FunT,'function_handle'))
                FunStr = func2str(FunT);
            else
                FunStr = FunT;
            end
            
            switch FunStr
                case 'FunOne'
                    DefV.X                    = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = ones(numel(InPar.X),1);
                    NPar  = 1;
                    if (nargout>2)
                        syms SymFunT;
                        SymFunT = 1;
                    end
                case 'FunX'
                    DefV.X                   = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = InPar.X;
                    NPar  = 1;
                    if (nargout>2)
                        syms X;
                        SymFunT = X;
                    end
                case 'FunMX'
                    DefV.X                   = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = -InPar.X;
                    NPar  = 1;
                    if (nargout>2)
                        syms X;
                        SymFunT = -X;
                    end
                case 'FunY'
                    DefV.Y                   = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = InPar.Y;
                    NPar  = 1;
                    if (nargout>2)
                        syms Y;
                        SymFunT = Y;
                    end
                case 'FunMY'
                    DefV.Y                   = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = -InPar.Y;
                    NPar  = 1;
                    if (nargout>2)
                        syms Y;
                        SymFunT = -Y;
                    end
                case 'FunXY'
                    DefV.X                   = [];
                    DefV.Y                   = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = InPar.X.*InPar.Y;
                    NPar  = 1;
                    if (nargout>2)
                        syms X Y;
                        SymFunT = X*Y;
                    end
                case 'FunTiltXp'
                    DefV.X                    = [];
                    DefV.Y                    = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = InPar.X.*(InPar.X+InPar.Y);
                    NPar  = 1;
                    if (nargout>2)
                        syms X Y;
                        SymFunT = [X*Y, X*X];
                    end
                case 'FunTiltYp'
                    DefV.X                    = [];
                    DefV.Y                    = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = InPar.Y.*(InPar.X+InPar.Y);
                    NPar  = 1;
                    if (nargout>2)
                        syms X Y;
                        SymFunT = [X*Y, Y*Y];
                    end
                case 'FunTiltXn'
                    DefV.X                    = [];
                    DefV.Y                    = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = InPar.X.*(InPar.X-InPar.Y);
                    NPar  = 1;
                    if (nargout>2)
                        syms X Y;
                        SymFunT = [-X*Y, X*X];
                    end
                case 'FunTiltYn'
                    DefV.X                    = [];
                    DefV.Y                    = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = InPar.Y.*(InPar.X-InPar.Y);
                    NPar  = 1;
                    if (nargout>2)
                        syms X Y;
                        SymFunT = [X*Y, -Y*Y];
                    end
                    
                case 'FunTiltX2p'
                    DefV.X                    = [];
                    DefV.Y                    = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = InPar.X.^2.*(InPar.X+InPar.Y);
                    NPar  = 1;
                    if (nargout>2)
                        syms X Y;
                        SymFunT = [X^3, X^2*Y];
                    end
                case 'FunTiltX2n'
                    DefV.X                    = [];
                    DefV.Y                    = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);
                    Out   = InPar.X.^2.*(InPar.X-InPar.Y);
                    NPar  = 1;
                    if (nargout>2)
                        syms X Y;
                        SymFunT = [X^3, -X^2*Y];
                    end    
                    
                    
                    
                case 'FunPolyXY'
                    DefV.X                    = [];
                    DefV.Y                    = [];
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);

                    N     = numel(InPar.X);
                    NPar  = size(Arg,2);
                    Out = zeros(N,NPar);
                    for I=1:1:NPar
                        Out(:,I) = InPar.X.^Arg(1,I) .* InPar.Y.^Arg(2,I);
                    end
                    if (nargout>2)
                        syms X Y;
                        SymFunT = Util.symbolic.symbolic_poly(Arg(1,:),X) .* Util.symbolic.symbolic_poly(Arg(2,:),Y);
                    end
                                        
                case 'FunPolyChebyshev1XY'
                    DefV.X                    = [];
                    DefV.Y                    = [];                    
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);

                    Kind = 1;

                    N    = numel(InPar.X);
                    NPar = size(Arg,2);

                    Out = zeros(N,NPar);
                    for I=1:1:NPar
                        Out(:,I) = Util.fun.chebyshev_poly(InPar.X,Arg(1,I),Kind).*Util.fun.chebyshev_poly(InPar.Y,Arg(2,I),Kind);
                    end
                    if (nargout>2)
                        syms X Y CX CY;
                        SymFunT = chebyshevT(Arg(1,:),X) .* chebyshevT(Arg(2,:),Y);
                    end
                case 'FunPolyChebyshev2XY'
                    DefV.X                    = [];
                    DefV.Y                    = [];                    
                    InPar = InArg.populate_keyval(DefV,varargin,mfilename,false);

                    Kind = 2;

                    N    = numel(InPar.X);
                    NPar = size(Arg,2);

                    Out = zeros(N,NPar);
                    for I=1:1:NPar
                        Out(:,I) = Util.fun.chebyshev_poly(InPar.X,Arg(1,I),Kind).*Util.fun.chebyshev_poly(InPar.Y,Arg(2,I),Kind);
                    end
                    if (nargout>2)
                        syms X Y CX CY;
                        SymFunT = chebyshevU(Arg(1,:),X) .* chebyshevU(Arg(2,:),Y);
                    end
                otherwise
                    error('Unknown Function option');
            end
        end
    end
    
    
    methods (Static)
        function [OutX,OutY]=apply_tran(X,Y,TranC,RotationUnits)
            % Apply a TranClass object or simple affine transformation
            % Package: @TranClass
            % Description: Apply a TranClass object or simple affine transformation
            % Input  : - Vector of X coordinates.
            %          - Vector of Y coordinates.
            %          - A TranClass object or a vectotr of
            %            [ShiftX, ShiftY]  (shift)
            %            or [ShiftX, ShiftY, Rot(rad)]  (rotation)
            %            or [ShiftX, ShiftY, RotX(rad), RotY(rad)]  (affine)
            %          - Rotation units (if TranC is a vector).
            %            Default is 'rad'.
            % Output : - Output X coordinates.
            %          - Output Y coordinates.
            
            if (nargin<4)
                RotationUnits    = 'rad';
            end

            if (isa(TranC,'TranClass'))
                % TranC is a TranClass object
                [OutX,OutY] = apply_tranclass(TranC, X, Y);
            else
                % Assume TranC is
                % [ShiftX, ShiftY]  (shift)
                % or [ShiftX, ShiftY, Rot(rad)]  (rotation)
                % or [ShiftX, ShiftY, RotX(rad), RotY(rad)]  (affine)
                switch numel(TranC)
                    case 2
                        % shift transformation
                        OutX = TranC(1) + X;
                        OutY = TranC(2) + Y;
                    case 3
                        % shift+rotation transformation
                        Rot  = convert.angular(RotationUnits,'rad',TranC(3));
                        OutX = TranC(1) + X.*cos(Rot) - Y.*sin(Rot);
                        OutY = TranC(2) + X.*sin(Rot) + Y.*cos(Rot);
                    case 4
                        % affine transformation
                        RotX  = convert.angular(RotationUnits,'rad',TranC(3));
                        RotY  = convert.angular(RotationUnits,'rad',TranC(4));
                        OutX = TranC(1) + X.*cos(RotX) - Y.*sin(RotX);
                        OutY = TranC(2) + X.*sin(RotY) + Y.*cos(RotY);
                    otherwise
                        error('Unknown TranC option');
                end
            end
        end
    end
    
    methods
        function TranC=poly_representation()
            %
            
            Ix = 1;
            Iy = 2;
            
            for Idim=1:1:2
                Nfun = numel(TranC.Fun{Idim});
                for Ifun=1:1:Nfun
                    Fun = TranC.Fun{Idim}{Ifun};
                    
                    switch func2str(Fun)
                        case 'FunPolyXY'
                            Arg = TranC.FunArg{Idim}{Ifun};
                            Poly = Arg;
                            Coef = ones(size(Arg));
                        case 'FunPolyChebyshev2XY'
                            [Cheb]=Util.fun.chebyshev_polyrep(2);
                            Arg = TranC.FunArg{Idim}{Ifun};
                            Narg = size(Arg,2);
                            if (isempty(TranC.Par))
                                Par = ones(1,Narg);
                            else
                                Par = TranC.Par{Idim}{Ifun};
                            end
                            
                            % Order is poly deg
                            
                            OrderVec = cell(size(Arg));
                            CoefVec  = cell(size(Arg));
                            
                            for Iarg=1:1:Narg
                                OrderX = Arg(1,Iarg);
                                OrderY = Arg(2,Iarg);
                                
                                %OrderVecX{Iarg}  = (OrderX:-1:0);
                                CoefVecX{Iarg}   = Cheb{OrderX+1};
                                CoefVecY{Iarg}   = Cheb{OrderY+1};
                                
                                CoefVecX{Iarg} = CoefVecX{Iarg}.*Par(Iarg);
                                CoefVecY{Iarg} = CoefVecY{Iarg}.*Par(Iarg);
                            end
                            
                        otherwise
                            error('Unknown function name');
                    end
                        
                end
            end
            
            
        end
        
    end
    
    
    methods
        % Construct the design matrix (per axis)
        function H=design_matrix(TranC,varargin)
            % Construct a design matrix for the TransClass object
            % Package: @TranClass
            % Description: Construct a design matrix for the TransClass
            %              object for a single axis.
            % Input  : - A TranClass object containing a single element.
            %            The TranClass object may contains multiple
            %            functions.
            %          * Arbitrary number of pairs of ...,key,val,...
            %            to transfer to the TranClass functions.
            %            e.g., ...'X',X,'Y','Mag',Mag,...
            % Output : - A cell array of design matrix. A design matrix per
            %            dimension.
            % Example: TranC = TranClass({@FunOne, [],@FunX,[],@FunY,[],@FunPolyXY,[2 3;1 2]}, {@FunOne, [],@FunX,[],@FunY,[]})
            %          X=rand(10,1); Y=rand(10,1);
            %          H=design_matrix(TranC,'X',X,'Y',Y)
            FunField     = TranClass.FunField;
            FunArgField  = TranClass.FunArgField;
            NParField    = TranClass.NParField;
            
            if (numel(TranC)>1)
                error('TranClass must contain a single element');
            end
            
            IC = 1;
            
            % for each dimension
            Ndim = numel(TranC(IC).(FunField));
            H    = cell(1,Ndim);
            for Idim=1:1:Ndim
                
            
                if (~iscell(TranC(IC).(FunField)))
                    TranC(IC).(FunField) = {TranC(IC).(FunField)};
                end

                % number of transformations
                Ntran = numel(TranC(IC).(FunField){Idim});

                for Itran=1:1:Ntran
                    FunName = TranC(IC).(FunField){Idim}{Itran};
                    Arg     = TranC(IC).(FunArgField){Idim}{Itran};
                    
%                     if (ischar(FunName)),
%                         FunStr  = FunName;
%                     elseif isa(FunName,'function_handle'),
%                         FunStr  = func2str(FunName);
%                     else
%                         error('TranClass object Fun field contains non function handle or string element');
%                     end

                    % Note that TranClass.(FunStr) is the static function name
                    Hadd = TranClass.evalfun(FunName,Arg,varargin{:});
                    %Hadd = TranClass.(FunStr)(varargin{:});
                    Npar = size(Hadd,2);  % number of free parameters per function
                    if (Itran==1)
                        H{Idim} = Hadd;
                    else
                        H{Idim} = [H{Idim}, Hadd];
                    end
                    TranC(IC).(NParField){Idim}{Itran} = Npar;

                end
            end
            
        end
        
        function [OutX,OutY] = apply_tranclass(TranC,X,Y)
            % Apply forward transformation to X,Y coordinates
            % Package: @TranClass
            % Description: Apply forward transformation to X,Y coordinates
            % Input  : - A TranClass object
            %          - A vector of X coordinates
            %          - A vector of Y coordinates
            % Output : - A vector of forward transformed X coordinates
            %          - A vector of forward transformed Y coordinates
            
            IC = 1;
            H = design_matrix(TranC,'X',X,'Y',Y);

            Nh = numel(H);
            if (Nh~=2)
                error('Number of dimensions must be 2');
            end
            % X par
            Idim = 1;
            %ParX = cell2mat(TranC(IC).Par{Idim}).';
            ParX = par2vector(TranC(IC),Idim).';
            % Y par
            Idim = 2;
            %ParY = cell2mat(TranC(IC).Par{Idim}).';
            ParY = par2vector(TranC(IC),Idim).';

            H = design_matrix(TranC(IC),'X',X,'Y',Y);

            Idim = 1;
            OutX = H{Idim}*ParX;
            Idim = 2;
            OutY = H{Idim}*ParY;
                    
        end
        
        function [OutX,OutY] = apply_tranclass_inv(TranC,X,Y,Thresh)
            % Apply inverse transformation to X,Y coordinates
            % Package: @TranClass
            % Description: Apply inverse transformation to X,Y coordinates
            % Input  : - A TranClass object
            %          - A vector of X coordinates
            %          - A vector of Y coordinates
            %          - Accuracy for convergence. Default is 1e-3.
            % Output : - A vector of inverse transformed X coordinates
            %          - A vector of inverse transformed Y coordinates
            
            IC = 1;

            if (nargin<4)
                Thresh               = 1e-3;
            end
            
            % 
            GuessX = X;
            GuessY = Y;

            % X par
            Idim = 1;
            %ParX = cell2mat(TranC(IC).Par{Idim}).';
            ParX = par2vector(TranC(IC),Idim).';
            % Y par
            Idim = 2;
            %ParY = cell2mat(TranC(IC).Par{Idim}).';
            ParY = par2vector(TranC(IC),Idim).';


            NotConverge = true;
            Counter = 0;
            while NotConverge
                Counter = Counter + 1;

                H = design_matrix(TranC(IC),'X',GuessX,'Y',GuessY);

                Idim = 1;
                X1 = H{Idim}*ParX;
                Idim = 2;
                Y1 = H{Idim}*ParY;

                DX = X1 - X;
                DY = Y1 - Y;

                %[DX, X, X1, GuessX]

                if (max(abs(DX))<Thresh && max(abs(DY))<Thresh)
                    % converged
                    NotConverge = false;
                end

                GuessX = GuessX - DX;
                GuessY = GuessY - DY;

            end

            OutX = GuessX;
            OutY = GuessY;
        end
                
        function TranC=populate_par(TranC,Idim,Par,ParErr)
            % Populate the parameters and errors in a TranClass object
            % Package: @TranClass
            % Description: Populate the parameters and errors in a TranClass object
            % Input  : - A single element TranClass object.
            %          - Index of dimension.
            %          - Vector of parameters.
            %          - Vector of errors in parameters
            % Output : - A TranClass object with the Par and ParErr field
            %            populated.
            % Example: TranC = TranClass({@FunOne,[],@FunX,[],@FunY,[],@FunTiltXp,[],@FunPolyXY,[1 2;2 2]});
            %          TranC=populate_par(TranC,1,[0.1 0.2 1 0 0.01 0.02]);
            
            FunField     = TranClass.FunField;
            NParField    = TranClass.NParField;
            ParField     = TranClass.ParField;
            ParErrField  = TranClass.ParErrField;
         
            if (nargin<4)
                ParErr = [];
            end
            
            if (numel(TranC)>1)
                error('TranClass object must have a single element');
            end
            
            IC = 1;
            Nfun = numel(TranC(IC).(FunField){Idim});
            K = 1;
            for Ifun=1:1:Nfun
                Npar = TranC(IC).(NParField){Idim}{Ifun};
                I1 = K;
                I2 = K+Npar-1;
                TranC(IC).(ParField){Idim}{Ifun} = Par(I1:I2);
                if (~isempty(ParErr))
                    TranC(IC).(ParErrField){Idim}{Ifun} = ParErr(I1:I2);
                end
                K  = K + Npar;
            end
            
            
        end
        
        function Par=par2vector(TranC,Idim)
            % Construct a vector of parameters from TranClass object
            % Description: Construct a vector of parameters from TranClass
            %            object in a specified dimension.
            % Input  : - A single element TranClass object.
            %          - Index of dimension.
            % Output : - A vector of parameters.
            % Example: Par=par2vector(TranC,1)
            % Reliable: 2
            
            
            ParField    = TranClass.ParField;
            NParField    = TranClass.NParField;
            
            if (numel(TranC)>1)
                error('TranClass object must have a single element');
            end
            IC = 1;
            if (isempty(TranC(IC).Par))
                % set all Par to 1
                TotalNPar = sum(cell2mat(TranC(IC).(NParField){Idim}));
                Par       = ones(1,TotalNPar);
            else
                % read parameters into a vector
                Par = cell2mat(TranC(IC).(ParField){Idim});
            end
                
        end
        
        function Par=getpar(TranC,Idim,Fun)
            % Get parameter for specific function in TranClass object
            % Package: @TranClass
            % Input  : - A TranClass object
            %          - Dimension index
            %          - Function name
            % Output : - Parameters for the specific function.
            %            If function name not found then empty.
            
            
            
            if (isa(Fun,'function_handle'))
                FunName = func2str(Fun);
            else
                FunName = Fun;
            end
            
            CellFunStr = cellfun(@func2str,TranC.Fun{Idim},'UniformOutput',false);
            Flag       = (strcmp(CellFunStr,FunName));
            if (any(Flag))
                Par    = TranC.Par{Idim}{Flag};
            else
                Par    = [];
            end
            
        end
        
        function [PolyS,OutOrder,OutCoef,ChPoly]=tranclass2poly(TranC,Idim,MultByPar,NormXY)
            % Convert a TranClass object into a polynomial representation
            % Description: Construct a polynomial representation for
            %              an TranClass object.
            % Input  : - A single TranClass object.
            %          - The dimension in the TranClass object for which
            %            to find the polynomial representation.
            %          - A flag (true|false) indicating if to multiply the 
            %            polynomials by their parameters (coeficiants).
            % Output : - A symbolic polynomial that represent the
            %            transformation.
            %          - Two rows matrix of the orders identified in the
            %            polynomials. Each column represent the order of
            %            the X and Y polynomials.
            %            For example: [1 2;0 2] represent:
            %            X*Y^2 + Y^2.
            %          - Vector of output coeficients that multiply each
            %            one of the polynomial parts.
            %          - A vector of the polynomial parts.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: TranC=TranClass({@FunOne,[],@FunPolyXY,[1 1 2;0 1 2]});
            %          TranC.Par{1}{1} = 30; TranC.Par{1}{2}=[0.01 0.02 0.03];
            %          [Poly,OutOrder,OutCoef,ChPoly]=tranclass2poly(TranC,1,true)
            % Reliable: 2
            
            
            if (nargin<4)
                NormXY = 1;
            end
            
            FunField    = TranClass.FunField;
            FunArgField = TranClass.FunArgField;
            ParField    = TranClass.ParField;
            %PolyField   = TranClass.PolyField;
            %ChPolyField = TranClass.ChPolyField;
            
            if (numel(TranC)>1)
                error('TranClass object must have a single element');
            end
            IC = 1;
            
            %IsPolyEmpty   = isempty(TranC(IC).(PolyField){Idim});
            %IsChPolyEmpty = isempty(TranC(IC).(ChPolyField){Idim});
            
            if isempty(TranC(IC).(ParField))
                if (MultByPar)
                    error('ParField must be populated for MultByPar true');
                end
            end
            TFun    = TranC(IC).(FunField){Idim};
            TFunArg = TranC(IC).(FunArgField){Idim};
            
            % for each function in the TranClass object:
            Nfun = numel(TFun);
            syms PolyS;
            for Ifun=1:1:Nfun
                [~,~,SymFun]=TranClass.evalfun(TFun{Ifun},TFunArg{Ifun});
                
                if (MultByPar)
                    % multiply each function by its coef.
                    SymFun = SymFun.*TranC(IC).(ParField){Idim}{Ifun}(:)';
                end
                if (Ifun==1)
                    PolyS = sum(SymFun);
                else
                    PolyS = PolyS + sum(SymFun);
                end
                
            end
            if (nargout>1)
                [OutOrder,OutCoef,ChPoly]=Util.symbolic.sympoly2d_2orders(PolyS);
                % renormalize coef.
                Nc = size(OutOrder,2);
                for Ic=1:1:Nc
                    ReNorm = NormXY.^OutOrder(1,Ic) .* NormXY.^OutOrder(2,Ic);
                    OutCoef(Ic) = OutCoef(Ic)./ReNorm;
                end 
                
            end
        end
        
        function TformPoly2D=tranclass2PolynomialTransformation2D(TranC,MultByPar)
            % Convert a TranClass object to a PolynomialTransformation2D object
            % Description: Convert a TranClass object into a 
            %              images.geotrans.PolynomialTransformation2D
            %              object.
            % Input  : - A TranClass object.
            %          - A flag (true|false) indicating if to multiply the 
            %            polynomials by their parameters (coeficiants).
            % Output : - A images.geotrans.PolynomialTransformation2D with
            %            the polynomial transformations.
            % Example: TranC=TranClass({@FunOne,[],@FunPolyXY,[1 1 2;0 1 2]},{@FunOne,[],@FunPolyXY,[1 1 2;0 1 2]});
            %          TranC.Par{1}{1} = 30; TranC.Par{1}{2}=[0.01 0.02 0.03];
            %          TranC.Par{2}{1} = 30; TranC.Par{2}{2}=[0.01 0.02 0.03];
            %          TformPoly2D=tranclass2PolynomialTransformation2D(TranC,true)
            % Reliable: 2
            
            FunField    = TranClass.FunField;
            
            Ntran = numel(TranC);
            for Itran=1:1:Ntran
                Ndim  = numel(TranC(Itran).(FunField));
                Coef  = cell(1,Ndim);
                for Idim=1:1:Ndim
                    
                    [~,OutOrder,OutCoef,~]=tranclass2poly(TranC(Itran),Idim,MultByPar);
                    % constrct an images.geotrans.PolynomialTransformation2D
                    % object
                    % Note that this object supports polynomials only to the
                    % 4th degree.
                
                    %U = A(1) + A(2).*X + A(3).*Y + A(4).*X.*Y + A(5).*X.^2 + A(6).*Y.^2 +...
                    TemplatePoly2DX = [0  1 0  1 2 0  2 1 3 0  2 3 1 4 0];
                    TemplatePoly2DY = [0  0 1  1 0 2  1 2 0 3  2 1 3 0 4];
                    TemplatePoly2D  = [TemplatePoly2DX', TemplatePoly2DY'];
                    
                    [~,J] = ismember(TemplatePoly2D,OutOrder','rows');
                    Coef{Idim} = Util.array.nangetind(OutCoef',J).';
                    Coef{Idim}(isnan(Coef{Idim})) = 0;
                end
                
                TformPoly2D(Itran) = images.geotrans.PolynomialTransformation2D(Coef{:});
            end
            
        end
        
        
        function Res=fit_transform(TranC,MatchedCat,varargin)
            % Fit a 2-D transformation to a list of matched coordinates
            % Description: Given a transformation object (TranClass) and
            %              a list of matched coordinates, find the best fit
            %              transformation between the lists of matched
            %              coordinates.
            %              The fit does not use user supplied  positional
            %              errors, but rather estimate the errors as a
            %              function of magnitude.
            %              OBSOLETE - use ImUtil.pattern.fit_transform
            % Input  : - An TranClass object.
            %            The TranClass object defines the transformation
            %            that will be fitted.
            %          - A MatchedCat structure generated by
            %            astcat2matched_array.m.
            %            This is a structure in which each field represent
            %            a column in the table of sources (e.g.,
            %            'XWIN_IMAGE'). The size of the matrix is number of
            %            sources (rows) by the number of epochs (columns).
            %            The matrix represent matched sources (e.g., the
            %            first line represent the first source as detected
            %            over all ephocs).
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'CooField' - Cell array of field names in the
            %                         MatchedCat structure containing the
            %                         X and Y coordinates, respectively.
            %                         The X and Y coordinates are used in
            %                         order to fit the transformation.
            %                         Default is:
            %                         {'XWIN_IMAGE','YWIN_IMAGE'}.
            %            'MagField' - A string indicating the magnitude
            %                         field in the MatchedCat structure.
            %                         Default is 'MAG_PSF'.
            %                         If empty then will not estimate the
            %                         astrometric error as a function of
            %                         magnitude and the fit will be
            %                         un-weighted. In this case the
            %                         AssymptoticRMS is not calculated.
            %            'ColorField'- A string indicating the color
            %                         field in the MatchedCat structure.
            %                         Default is empty.
            %            'UserFlagMat'- A matrix of true|false flags of
            %                         size number of sources by number of
            %                         epochs, indicating which sources to
            %                         use in the fit. This can be used to
            %                         exclude from the fit problematic
            %                         sources (e.g., saturated, extended).
            %                         If empty, then all true.
            %                         Default is empty.
            %            'SelectRefMethod' - The epoch index to use as the
            %                         reference epoch (i.e., the
            %                         transformation will be fitted
            %                         relative to this epoch).
            %                         Alternatively one of the following
            %                         methods:
            %                         'max' - use the epoch with the
            %                                 largest number of sources.
            %                         Default is 'max'.
            %            'MethodLSQ' - Least square fitting method.
            %                         Options are:
            %                         '\'|'lscov'|
            %                         {'pcg','cgs','bicg','bicgstab','bicgstabl','qmr','minres','gmres'}
            %                         Default is 'lscov'.
            %            'lscovPar' - Cell array of additional parameters
            %                         to pass to lscov.m.
            %                         Default is {}.
            %            'ls_conjgradPar' - Cell array of additional parameters
            %                         to pass to ls_conggrad.m.
            %                         Default is {}.
            %            'MagBinSize' - The magnitude bin size for
            %                         estimating the astrometric error as a
            %                         function of magnitude.
            %                         Default is 0.5 mag.
            %            'MagBinMaxErr' - The maximum astrometric error
            %                         allowed for bins in the rms. vs.
            %                         magnitude histogram. Bins with larger
            %                         error will be excluded.
            %                         Default is 0.5 mag.
            %            'MagBinMinNpt' - Minimum number of points allowed
            %                         for bins in the rms. vs. magnitude
            %                         histogram. Bins with smaller
            %                         number of points will be excluded.
            %                         Default is 2.
            %            'ClipNsigma' - Sigma clipping in sigmas.
            %                         Default is 2.
            %                         This will be applied only if number
            %                         of iteration is >1.
            %            'Niter' -    Number of sigma clipping iterations.
            %            'Verbose' -  Verbose. Default is true.
            %            'VerbosePlot' - Plot rms vs. mag.
            %                         Default is false.
            % Output : - A structure array of results. One element per
            %            epoch. The following fields are available:
            %            'TranC' - A TranClass object with the fitted
            %                      transformation.
            %            'Dim'   - Structure array with element per
            %                      dimension.
            %            'Resid' - Vector of all residuals.
            %            'Flag'  - Vector of flags indicating if source was
            %                      used in the fit.
            %            'wmed_rms'-weighted median rms of best fit [pix].
            %            'rms'     - rms of best fit [pix].
            %            'rrms'    - robust rms of best fit [pix].
            %            'rmsAll'  - rms of all sources [pix].
            %            'AssymptoticRMS' - minimum of the rms fitted to
            %                        the rms vs. mag plot [pix].
            % % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Jan 2017
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: S=FITS.read2sim('test_PTF_*f02_p004162_c08.fits');
            %          S=mextractor(S);
            %          [AstM]=match(S);
            %          [MatchedCat]=astcat2matched_array(AstM,{'XWIN_IMAGE','YWIN_IMAGE','MAG_PSF'});
            %          TranC = TranClass({@FunOne, [],@FunX,[],@FunY,[],@FunTiltXp,[],@FunTiltXn,[]}, {@FunOne, [],@FunX,[],@FunY,[],@FunTiltYp,[],@FunTiltYn,[] })
            %          Res=fit_transform(TranC,MatchedCat)
            %          %plot residuals:
            %          %X vs ResX
            %          plot(Res(2).Dim(1).Coo(Res(2).Flag), Res(2).Dim(1).Resid(Res(2).Flag),'.')
            %          %Y vs. ResX
            %          plot(Res(2).Dim(2).Coo(Res(2).Flag), Res(2).Dim(1).Resid(Res(2).Flag),'.')
            %          %mark stars that were used for the best fit transformation
            %          %ds9(S(2)); ds9.plot([Res(2).Dim(1).Coo( Res(2).Flag ), Res(2).Dim(2).Coo( Res(2).Flag )]);
            % Reliable: 2
            
            
            DefV.CooField             = {'XWIN_IMAGE','YWIN_IMAGE'};
            DefV.MagField             = 'MAG_PSF';
            DefV.ColorField           = [];
            DefV.UserFlagMat          = [];
            DefV.SelectRefMethod      = 'max';   % index, 'max'
            DefV.MethodLSQ            = 'lscov'; % '\'|'lscov','pcg','cgs','bicg','bicgstab','bicgstabl','qmr','minres','gmres'
            DefV.lscovPar             = {};
            DefV.ls_conjgradPar       = {};
            DefV.MagBinSize           = 0.5;
            DefV.MagBinMaxErr         = 0.5;
            DefV.MagBinMinNpt         = 2;
            DefV.ClipNsigma           = 2;
            DefV.Niter                = 2;
            DefV.Verbose              = true;
            DefV.VerbosePlot          = false;
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

          
            Ndim = numel(InPar.CooField);
            for Idim=1:1:Ndim
                Mat(Idim).Coo = MatchedCat.(InPar.CooField{Idim});
            end
            [Nsrc,Nep] = size(MatchedCat.(InPar.CooField{1}));

            if (~isempty(InPar.MagField))
                MatMag = MatchedCat.(InPar.MagField);
            else
                MatMag = nan(Nsrc,Nep);
            end
            if (~isempty(InPar.ColorField))
                MatColor = MatchedCat.(InPar.ColorField);
            else
                MatColor = nan(Nsrc,Nep);
            end

            % define the user flag matrix
            % this matrix indicating which sources to use in the fit
            if (isempty(InPar.UserFlagMat))
                InPar.UserFlagMat = true(Nsrc,Nep);
            end

            % select the column index for the reference image
            if (ischar(InPar.SelectRefMethod))
                switch lower(InPar.SelectRefMethod)
                    case 'max'
                        [~,RefInd] = max(sum(~isnan(MatchedCat.(InPar.CooField{1})),1)./Nsrc);
                    otherwise
                        error('Unknown SelectRefMethod option');
                end
            else
                % assume InPar.SelectRefMethod is and index:
                RefInd = InPar.SelectRefMethod;
            end


            % The ref catalog properties
            for Idim=1:1:Ndim
                Ref(Idim).Coo = Mat(Idim).Coo(:,RefInd);
            end
            
            Mref = MatMag(:,RefInd);
            Cref = MatColor(:,RefInd);

            
            
            
            % contruct the design matrix
            H        = design_matrix(TranC,'X',Ref(1).Coo,'Y',Ref(2).Coo,'Mag',Mref,'Color',Cref);
            FlagH    = ~any(isnan(H{1}),2);
                            
            % find transformation for each epoch separately
            Res = Util.struct.struct_def({'TranC','Dim','Resid','Flag','wmed_rms','rms','rrms','rmsAll','AssymptoticRMS'},Nep,1);
            for Iep=1:1:Nep
                %Iep
                % for each epoch
                Res(Iep).TranC = TranC;
                New = Util.struct.struct_def({'Coo','CooErr'},Ndim,1);
                for Idim=1:1:Ndim
                    New(Idim).Coo    = Mat(Idim).Coo(:,Iep);
                    New(Idim).CooErr = ones(size(New(Idim).Coo));
                    Res(Iep).Dim(Idim).Coo = Mat(Idim).Coo(:,Iep);
                end
                
                M = MatMag(:,Iep);
                C = MatColor(:,Iep);
                
                % Initiate FlagClip (all true)
                FlagClip = true(size(FlagH));

                % sigma clipping and error estimation - iterations
                Iiter = 0;
                Cont  = true;
                while Iiter<InPar.Niter && Cont
                    Iiter = Iiter + 1;
                
                %for Iiter=1:1:InPar.Niter
                    
                    % for each dimension
                    Res(Iep).Resid = zeros(Nsrc,1);
                    for Idim=1:1:Ndim

                        Flag     = FlagH & ~isnan(New(Idim).Coo) & FlagClip & InPar.UserFlagMat(:,Iep);
                        
                        Res(Iep).Flag = Flag;

                        switch lower(InPar.MethodLSQ)
                            case 'lscov'
                                [Res(Iep).Dim(Idim).Par,Res(Iep).Dim(Idim).ParErr] = lscov(H{Idim}(Flag,:), New(Idim).Coo(Flag), 1./(New(Idim).CooErr(Flag).^2),InPar.lscovPar{:});
                            case '\'
                                Res(Iep).Dim(Idim).Par = H{Idim}\New(Idim).Coo;
                            otherwise
                                % use ls_conjgrad
                                 Out = ls_conjgrad(H{Idim}(Flag,:), New(Idim).Coo(Flag), New(Idim).CooErr(Flag), InPar.MethodLSQ,InPar.ls_conjgradPar{:});
                                 Res(Iep).Dim(Idim).Par    = Out.Par;
                                 Res(Iep).Dim(Idim).ParErr = Out.ParErr;
                        end
                        % Best fit residuals
                        Res(Iep).Dim(Idim).Resid = New(Idim).Coo - H{Idim}*Res(Iep).Dim(Idim).Par;

                        Res(Iep).Resid = Res(Iep).Resid + Res(Iep).Dim(Idim).Resid.^2;
                    end  % end Idim
                    % estimate the errors as a function of magnitude
                    % combine residuals
                    Res(Iep).Resid = sqrt(Res(Iep).Resid);
                    %semilogy(Mref,Res(Iep).Resid,'.')
                    
                    % populate the RMS information
                    Res(Iep).wmed_rms = Util.stat.wmedian(Res(Iep).Resid(Flag),New(Idim).CooErr(Flag)); %Err(Flag));
                    Res(Iep).rms      = std(Res(Iep).Resid(Flag));
                    Res(Iep).rrms     = Util.stat.rstd(Res(Iep).Resid(Flag));
                    Res(Iep).rmsAll   = nanstd(Res(Iep).Resid);
                    
                    if (Res(Iep).wmed_rms<1e-10)
                        % residual is tiny - no need to go for 2nd
                        % iteration
                        Cont = false;
                    end
                    
                    %--- calculate Errors and Flags for next iteration ---
                    if (isempty(InPar.MagField))
                        % user did not supply the MagField in the
                        % MatchedCat structure - do not apply errors as a
                        % function of magnitude
                        % Use simple sigma clipping instead
                        
                        SigmaClipLimit = nanmedian(Res(Iep).Resid) + Util.stat.rstd(Res(Iep).Resid).*InPar.ClipNsigma;
                        FlagClip       = Res(Iep).Resid<SigmaClipLimit;
                        
                        Err = ones(Nsrc,1);
                        for Idim=1:1:Ndim
                            New(Idim).CooErr = Err;
                            Res(Iep).Dim(Idim).CooErr = Err;
                        end
                    else
                        % use Magnitude in order to estimate the
                        % astrometric errors as a function of magnitude
                        
                        B = timeseries.binning([Mref,Res(Iep).Resid],InPar.MagBinSize,[NaN NaN],{'MidBin',@numel,@median,@Util.stat.rstd});

                        FlagBin = B(:,2)>InPar.MagBinMinNpt & B(:,3)<InPar.MagBinMaxErr;
                        B = B(FlagBin,:);
                        if (isempty(B))
                            % do something about it...
                        end
                            
                        BclipInterp = interp1(B(:,1), B(:,3)+B(:,4).*InPar.ClipNsigma, Mref,'linear');


                        %hold on;
                        %semilogy(Mref(FlagClip),Res(Iep).Resid(FlagClip),'.')
                        Err = interp1(B(:,1),B(:,3),Mref,'linear');
                        %New(Idim).CooErr = interp1(B(:,1),B(:,3),Mref,'linear');
                        for Idim=1:1:Ndim
                            New(Idim).CooErr = Err;
                            Res(Iep).Dim(Idim).CooErr = Err;
                        end
                        FlagClip = Mref>min(B(:,1)) & Mref<max(B(:,1)) & Res(Iep).Resid<BclipInterp & ~isnan(New(Idim).CooErr);
                    end

                end % end Iiter
                
                if (InPar.VerbosePlot)
                    semilogy(Mref(FlagClip),Res(Iep).Resid(FlagClip),'.')
                    hold on;
                end
                
                % Calculate asymptotic rms
                % the asymtiotitc rms is defined as the minimum of the
                % function describes the rms vs. mag.
                if (Iep==RefInd || isempty(InPar.MagField))
                    % If current image is the ref image than skip
                    Res(Iep).AssymptoticRMS = NaN;
                else
                    
                    P=polyfit(Mref(FlagClip),log10(Res(Iep).Resid(FlagClip)),2);
                    ExtramumPX = -P(2)./(2.*P(1));
                    if (ExtramumPX>min(B(:,1)) && ExtramumPX<max(B(:,1)))
                        ExtramumPY = 10.^polyval(P,ExtramumPX);
                        Res(Iep).AssymptoticRMS = ExtramumPY;
                    else
                        Res(Iep).AssymptoticRMS = NaN;
                    end

                end
                
                % Populate the output TranClass object
                Res(Iep).TranC = TranC;
                %Res(Iep).TranC.FunArg = 
                for Idim=1:1:Ndim
                    Res(Iep).TranC = populate_par(Res(Iep).TranC,Idim, Res(Iep).Dim(Idim).Par, Res(Iep).Dim(Idim).ParErr);
                end
            end % end Iep
            






            
            
        end
        
        function WCS=tranclass2wcs(TranC,NormXY,varargin)
            % OBSOLETE
            
            [~,OutOrderX,OutCoefX,~]=tranclass2poly(TranC,1,true);
            [~,OutOrderY,OutCoefY,~]=tranclass2poly(TranC,2,true);
            
            
        end
        
    end
    
end

            
