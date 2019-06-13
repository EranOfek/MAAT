%--------------------------------------------------------------------------
% SIM class                                                          class
% Description: A class of structure array of images (SIM).
%              This is an in-memory image container and calculator.
%              Note that this class is a superset of the AstCat class.
%              Fields are:
%              .Im 
%              .Header
%              .ImageFileName
%              .PSF
%              .Mask
%              .BackIm
%              .ErrIm
%              .WeightIm
%              + the fields of AstCat
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------


classdef SIM < AstCat & MASK & ClassPSF % & IM
    properties (SetAccess = public)
        Im            
        ImageFileName = '';
        BackIm        
        ErrIm
        WeightIm
        %PSF = ClassPSF;
    end

    %-----------------------
    %--- SIM constructor ---
    %-----------------------
    methods
        
        function Sim=SIM(N,M)
            % SIM constructor method
            % Package: @SIM
            % Package Category: class
            
            ImageField = 'Im';
            
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
                    Sim(I,J).(ImageField) = [];
                end
            end
        end

    end
    
    %--------------------
    %--- Clear fields ---
    %--------------------
    methods
        function Sim=clear_back(Sim)
            % clear background image in SIM
            % Package: @SIM
            % Package Category: class
            % Input  : - A SIM object
            % Output : - A SIM object
            
            BackField = SIM.BackField;
            Nsim = numel(Sim);
            for Isim=1:1:Nsim
                Sim(Isim).(BackField) = [];
            end

        end

        function Sim=clear_err(Sim)
            % clear background image in SIM
            % Package: @SIM
            % Package Category: class
            % Input  : - A SIM object
            % Output : - A SIM object
            
            ErrField = SIM.ErrField;
            Nsim = numel(Sim);
            for Isim=1:1:Nsim
                Sim(Isim).(ErrField) = [];
            end

        end
    end
    
  
    %----------------------------------------
    %--- Static methods for the Class SIM ---
    %----------------------------------------
    % Static methods for field names
    methods (Static)
 
        function Ans=issim(Obj)
            % Return true if object is SIM
            % Package: @SIM
            % Package Category: class
            % Description: Check if object is of SIM class.
            % Input  : - Object
            % Output : - {true|false}.
            %     By : Eran O. Ofek                    Oct 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: S=1; SIM.issim(S);
            % Reliable: 2
            Ans = isa(Obj,'SIM');

         end
        
        function Sim=scalar2sim(Array,varargin)
            % Array of scalars to SIM of scalars
            % Package: @SIM
            % Package Category: convert
            % Description: Convert an array of scalars to a SIM of scalars.
            %              Each element in the array will be stored in a SIM
            %              image consist of one scalar.
            % Input  : - An array of numerics.
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'TargetField' - A string or a cell array of string containing
            %                            the field names to which to copy the array
            %                            values.
            % Output : - A SIM object in which each elemnt contains an image
            %            with the corresponding scalar from the array.
            % License: GNU general public license version 3
            %     By : Eran O. Ofek                    Apr 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Sim=SIM.scalar2sim(rand(10,3));
            % Reliable: 2


            DefV.TargetField       = SIM.ImageField;
            if (isempty(varargin))
                InPar = DefV;
            else
                InPar = InArg.populate_keyval(DefV,varargin,mfilename);
            end

            if (~iscell(InPar.TargetField))
                InPar.TargetField = {InPar.TargetField};
            end
            Nf = numel(InPar.TargetField);
            Sim  = simdef(size(Array));
            Nsim = numel(Array);
            for Isim=1:1:Nsim
                for If=1:1:Nf
                    Sim(Isim).(InPar.TargetField{If}) = Array(Isim);
                end
            end

        end
        
        function OutSim=struct2sim(InSim)
            % Convert Structure to SIM
            % Input  : - A structure array with the SIM fields.
            % Output : - A SIM object
            
            if (SIM.issim(InSim))
                OutSim = InSim;
            else
                FN     = fieldnames(InSim);
                Nfn    = numel(FN);

                OutSim = SIM;
                Nsim   = numel(InSim);

                for Isim=1:1:Nsim
                    for Ifn=1:1:Nfn
                        OutSim(Isim).(FN{Ifn}) = InSim(Isim).(FN{Ifn});
                    end
                end
            end

        end
            
        function Name=ImageField
            % Return the image field name in the SIM class.
            % Package: @SIM
            % Package Category: class
            % Input  : null
            % Output : Image field name.
            Name = 'Im';
        end
        
        function Name=BackField
            % Return the image background field name in the SIM class.
            % Package: @SIM
            % Package Category: class
            % Input  : null
            % Output : Image background field name.
            Name = 'BackIm';
        end
        
        function Name=ErrField
            % Return the image error field name in the SIM class.
            % Package: @SIM
            % Package Category: class
            % Input  : null
            % Output : Image error field name.
            Name = 'ErrIm';
        end
       
        function Name=WeightField
            % Return the image weight field name in the SIM class.
            % Package: @SIM
            % Package Category: class
            % Input  : null
            % Output : Image weight field name.
            Name = 'WeightIm';
        end
        
        function Name=FileNameField
            % Return the image file name field name in the SIM class.
            % Package: @SIM
            % Package Category: class
            % Input  : null
            % Output : Image file name field name.
            Name = 'ImageFileName';
        end
    end  % methods static
    
    
    % methods for image functions
    methods 
        
        % vector cut
        function Res=vector_prof(Sim,Start,End,varargin)
            % Image values along a line defined by two points.
            % Package: @SIM
            % Description: Return the SIM image values along a line defined
            %              by two points.
            % Input  : - A SIM object.
            %          - Two column matrix of start points [X,Y].
            %            Each line corresponds to one image.
            %          - Two column matrix of end points [X,Y].
            %          * Arbitrary number of ...,key,val,... pairs.
            %            Available keywords are:
            %            'Plot' - Plot the lines. Default is false.
            %            'PlotPar' - Cell array of additional parameters to
            %                     pass to the plot function. Default is {}.
            %            'InterpMethod' - Interpolation method.
            %                     See interp1.m for options.
            %                     Default is 'linear'.
            %            'Step'    - Steps along the line. Default is 1.
            %            'AxisType'- X axis in the plot:
            %                     'x' - plot intensity as a function of x.
            %                     'y' - plot intensity as a function of y.
            %                     'd' - plot intensity as a function of
            %                           distance from line edge.
            %                     Default is 'd'.
            %            'Field'   - SIM field from which to get the data.
            %                        Default is 'Im'.
            % Output : - A structure array (element per image or line) with
            %            the following fields:
            %            .LineV  - Values along the vector.
            %            .LineX  - X coordinate along the vector.
            %            .LineY  - Y coordinate along the vector.
            % Example: Res=vector_prof(Sim,[1 1],[6 10]);
            % Reliable: 2
            
            
            DefV.Plot                = false;
            DefV.PlotPar             = {};
            DefV.InterpMethod        = 'linear';
            DefV.Step                = 1;
            DefV.AxisType            = 'd';
            DefV.Field               = SIM.ImageField;
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            Npt = size(Start,1);
            Nim = numel(Sim);
            if ~(Npt==Nim || Npt==1 || Nim==1)
                error('Number of points and number of images is inconsistent');
            end
            N   = max(Npt,Nim);
            
            % allocate output structure
            Res = Util.struct.struct_def({'LineV','LineX','LineY'},N,1);
            for I=1:1:N
                Ipt = min(I,Npt);
                Iim = min(I,Nim);
                
                Lengths = abs(End-Start);   % vector [X,Y] lengths
                if (Lengths(1)>Lengths(2))
                    LineX = (Start(Ipt,1):InPar.Step:End(Ipt,1))';
                    LineY = interp1( [Start(Ipt,1),End(Ipt,1)], [Start(Ipt,2),End(Ipt,2)], LineX,'linear');
                else
                    LineY = (Start(Ipt,2):InPar.Step:End(Ipt,2))';
                    LineX = interp1( [Start(Ipt,2),End(Ipt,2)], [Start(Ipt,1),End(Ipt,1)], LineY,'linear');
                end
                
                Res(I).LineV = interp2(Sim(Iim).(InPar.Field),LineX,LineY,InPar.InterpMethod);
                
                Res(I).LineX = LineX;
                Res(I).LineY = LineY;
                
                if (InPar.Plot)
                    switch lower(InPar.AxisType)
                        case 'x'
                            plot(Res(I).LineX,Res(I).LineV,InPar.PlotPar{:});
                            H = xlabel('X [pix]');
                            H.FontSize = 18;
                        case 'y'
                            plot(Res(I).LineY,Res(I).LineV,InPar.PlotPar{:});
                            H = xlabel('Y [pix]');
                            H.FontSize = 18;
                        case 'd'
                            Dist = sqrt((LineX-LineX(1)).^2 + (LineY-LineY(1)).^2);
                            plot(Dist,Res(I).LineV,InPar.PlotPar{:});
                            H = xlabel('Distance [pix]');
                            H.FontSize = 18;
                        otherwise
                            error('Unknown AxisType option');
                    end
                    hold on;
                    grid on;
                    H = ylabel('Value');
                    H.FontSize = 18;
                    hold off;
                    drawnow;
                end
            end
            
        end
        
        % radial profile
        function Res=rad_prof(Sim,Coo,Radius,varargin)
            % Mean radial profile around points in SIM images
            % Package: @SIM
            % Description: Calculate and plot the mean radial profile
            %              around points in SIM images.
            % Input  : - A SIM object.
            %          - A two column matrix of [X,Y] coordinates around
            %            which to calculate the radial profile.
            %            Each line corresponds to one image.
            %          - Radius. Default is 15.
            %          * Arbitrary number of ...,key,val,... pairs.
            %            Available keywords are:
            %            'Plot' - Plot the lines. Default is false.
            %            'PlotPar' - Cell array of additional parameters to
            %                     pass to the plot function. Default is {}.
            %            'InterpMethod' - Interpolation method.
            %                     See interp1.m for options.
            %                     Default is 'linear'.
            %            'Bin'     - Radial bin size. Default is 1.
            %            'MeanFun' - Function ti use in calculating the
            %                        radial profile. Default is @mean.
            %            'Field'   - SIM field from which to get the data.
            %                        Default is 'Im'.
            % Output : - A structure array (element per image or line) with
            %            the following fields:
            %            .Rad  - Radius.
            %            .Val  - Value.
            % Example: Res=rad_prof(Sim,[60 67.2]);
            % Reliable: 
            
            if (nargin==2)
                Radius = 15;
            end
            
            DefV.Plot                = false;
            DefV.PlotPar             = {};
            DefV.InterpMethod        = 'linear';
            DefV.Bin                 = 1;
            DefV.MeanFun             = @nanmean;
            DefV.Field               = SIM.ImageField;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            Npt = size(Coo,1);
            Nim = numel(Sim);
            if ~(Npt==Nim || Npt==1 || Nim==1)
                error('Number of points and number of images is inconsistent');
            end
            N   = max(Npt,Nim);
            
            RadVec      = (0:InPar.Bin:Radius);
            Nrad        = numel(RadVec);
            % allocate output structure
            Res = Util.struct.struct_def({'Rad','Val'},N,1);
            for I=1:1:N
                Ipt = min(I,Npt);
                Iim = min(I,Nim);
                
                % trim image around coordinate
                TrimIm      = ImUtil.Im.trim_image(Sim(Iim).(InPar.Field),[Coo(Ipt,:), Radius, Radius],'center');
                Size        = size(TrimIm);
                [MatX,MatY] = meshgrid((1:1:Size(2)),(1:1:Size(1)));
                Xcen = Radius + 1;
                Ycen = Radius + 1;
                MatR        = sqrt((MatX-Xcen).^2 + (MatY-Ycen).^2);
                
                % radial profile
                Res(I).Rad = RadVec(1:end-1)+0.5.*InPar.Bin;
                Res(I).Val = zeros(1,Nrad-1);
                for Irad=1:1:Nrad-1
                    Flag = (MatR>=RadVec(Irad) & MatR<RadVec(Irad+1));
                    Res(I).Val(Irad) = InPar.MeanFun(TrimIm(Flag));
                end
                    
              if (InPar.Plot)
                    plot(Res(I).Rad,Res(I).Val,InPar.PlotPar{:});
                    hold on;
                    grid on;
                    H = xlabel('R [pix]');
                    H.FontSize = 18;
                    H = ylabel('Value');
                    H.FontSize = 18;
                    hold off;
                    drawnow;
                end
            end
            
        end
        
        % surface plot
        function TrimIm=local_surface(Sim,Coo,HalfSize,varargin)
            % Local surface around a point in a single SIM image
            % Package: @SIM
            % Description: Present the local surface around a point in a
            %              single SIM image.
            % Input  : - A single element SIM object.
            %          - A two column matrix of [X,Y] coordinates around
            %            which to calculate the radial profile.
            %            Each line corresponds to one image.
            %          - Half size [X,Y] for surface region.
            %            Default is [15 15].
            %          * Arbitrary number of ...,key,val,... pairs.
            %            Available keywords are:
            %            'Plot' - Plot the surface. Default is true.
            %            'PlotPar' - Cell array of additional parameters to
            %                     pass to the plot function. Default is {}.
            %            'PlotFun' - plot function (e.g., @surface,
            %                     @contour). Default is @surface.
            %            'Field'   - SIM field from which to get the data.
            %                        Default is 'Im'.
            % Output : - A SIM object with the trimmed image around the
            %            requested coordinates.
            % Example: Res=local_surface(Sim,[60 67.2]);
            % Reliable: 
            
            if (nargin==2)
                HalfSize = [15 15];
            end
            
            if (numel(HalfSize)==1)
                HalfSize = [HalfSize, HalfSize];
            end
            
            DefV.Plot                = true;
            DefV.PlotPar             = {};
            DefV.PlotFun             = @surface;
            DefV.Field               = SIM.ImageField;
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            Npt = size(Coo,1);
            Nim = numel(Sim);
            if ~(Npt==1 || Nim==1)
                error('Number of points and number of images is inconsistent');
            end
            
            % works on single image/coordinate
            TrimIm = trim_image(Sim,[Coo HalfSize],'SectionMethod','center');
            
            if (InPar.Plot)
                InPar.PlotFun(TrimIm.(InPar.Field));
                H = colorbar;
                ylabel(H,'Intensity');
                
                H = xlabel('X [pix]');
                H.FontSize = 18;
                H = ylabel('Y [pix]');
                H.FontSize = 18;
                H = zlabel('Value');
                H.FontSize = 18;
                drawnow;
            end
        end
        
        % central moments around a single point (see ImUtil.Im.im_moments)
        function [AstC,Col]=moments(Sim,Coo,varargin)
            % 1st and 2nd moments around specific positions in SIM images
            % Package: @SIM
            % Description: Calculate 1st and 2nd moments around specific
            %              positions in SIM images.
            % Input  : - A SIM object.
            %          - Two column matrix of [X,Y] coordinates in which
            %            to measure moments.
            %          * Arbitrary number of ...,key,val,... pairs.
            %            Available keywords are:
            %            'Radius'  - Radius in which to calculate the moments
            %                        and properties. Default is 10.
            %            'Sigma'   - Sigma of Gaussian by which to weight the
            %                        pixels. Default is 1.5.
            %            'MaxIter' - Maximum number of position iterations.
            %                        If 1 then use initial coordinates.
            %                        Default is 3.
            %            'Field'   - SIM field from which to get the data.
            %                        Default is 'Im'.
            %            'OutType' - 'AstCat'|'mat'. Default is 'AstCat'.
            %                        Note that 'mat' will work only for a
            %                        single element SIM object.
            % Output : - An AstCat object containing the following columns:
            %            'XWIN_IMAGE' - 1st central weighted moment in X.
            %            'YWIN_IMAGE' - 1st central weighted moment in Y.
            %            'X2WIN_IMAGE'- 2nd central weighted moment in X^2.
            %            'Y2WIN_IMAGE'- 2nd central weighted moment in Y^2.
            %            'XYWIN_IMAGE'- 2nd central weighted moment in X*Y.
            %          - Structure of column names and indices.
            % Example: [AstC]=moments(Sim,[1;2],[101;200]);
            % Reliable: 2
            
            
            CatField     = AstCat.CatField;
            ColField     = AstCat.ColField;
            ColCellField = AstCat.ColCellField;
            
            DefV.Radius              = 10;
            DefV.Sigma               = 1.5;
            DefV.MaxIter             = 3;
            DefV.Field               = SIM.ImageField;
            DefV.OutType             = 'AstCat';
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            Nsim = numel(Sim);
            AstC = AstCat(size(Sim));
            ColCell = {'XWIN_IMAGE','YWIN_IMAGE','X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE'};
            for Isim=1:1:Nsim
                % for each image
                % calculate moments
                [Mom,Mom2,~] = ImUtil.Im.im_moments(Sim(Isim).(InPar.Field),Coo(:,1),Coo(:,2),InPar.Radius,InPar.Sigma,InPar.MaxIter);
                % populate catalog
                AstC(Isim).(CatField) = [Mom.X, Mom.Y, Mom2.X2, Mom2.Y2, Mom2.XY];
                AstC(Isim).(ColCellField) = ColCell;
            end
            AstC = colcell2col(AstC);
            if (nargout>1)
                Col  = AstC(1).(ColField);
            end
            
            switch lower(InPar.OutType)
                case 'mat'
                    AstC = AstC.(CatField);
                otherwise
                    % do nothing
            end
            
        end
        
        
    end 
    
    
    methods
        function [N,X]=hist(Sim,varargin)
            % hist function for pixels in a SIM image
            % Package: @SIM
            
            if (numel(Sim)>1)
                error('hist is defined for a single SIM class image');
            end
            [N,X] = hist(Sim.Im(:),varargin);
        end
        
        function [N,X]=histc(Sim,varargin)
            % histc function for pixels in a SIM image
            % Package: @SIM
            if (numel(Sim)>1)
                error('histc is defined for a single SIM class image');
            end
            [N,X] = histc(Sim.Im(:),varargin);
        end
        
        
    
        

        %--- Structre functions ---
        function obj=isfield(Sim,Field)
            % isfield function for SIM object
            % Package: @SIM
            obj = any(strcmp(fieldnames(Sim),Field));
        end

        function obj=isstruct(Sim)
            % return true for isstruct of SIM
            % Package: @SIM
            obj = true;  %isstruct(Sim) || isa(Sim,'SIM');
        end

        

        
        % additional functions to add:
        % flip ???
        % fliplr
        % rot90
        % imagesize
        % resize
        % get_stamp
        % crdetect
        % saturated
        
        % near_coo
        % nearest_xy
        % nearest_coo
        % addextcat
        % ds9
        % ds9_regions
        % ds9_cat
        % ds9_extcat
        % back_std
        % imshow
        % filter/conv
        % getpsf (need a PSF class)
        % align
        % matchcat
        % relphot
        % relastrom
        % addastrom
        % twflat
        % domeflat
        % skyflat
        % coadd
        % combine
        % subtract
        % 
        
    end
end

            
