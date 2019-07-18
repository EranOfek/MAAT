%--------------------------------------------------------------------------
% Kernel2 class                                                      class
% Description: A static class for 2-D kernels (functions)
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Aug 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef Kernel2 
             
    % 2-D functions
    methods (Static)
        
        % gaussian
        function K=gauss(varargin)
            % Generate a Gaussian kernel with zero background.
            % Package: @Kernel2
            % Description: Generate a Gaussian kernel with zero background and
            %              normalized to unity sum.
            % Input  : - SigmaX of Gaussian. Default is 1.5 pix.
            %          - SigmaY of Gaussian. Default is 1.5 pix.
            %          - Rho of Gaussian. Default is 0.
            %          - Stamp size in X direction. Default is 11.
            %          - Stamp Size in Y direction. Default is 11.
            % Output : - A Gaussian kernel located at the stamp center.
            %            The stamp center is defined to be the middle of the central
            %            pixel (i.e., delta function span over one pixel), where the
            %            central pixel is give by ceil(Size./2 + 0.1).
            % See also: kernel_aper.m, kernel_annulus.m, kernel_exp.m
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Jun 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: K=Kernel2.gauss;
            % Reliable: 2
            %--------------------------------------------------------------------------

            %         SigmaX SigmaY Rho SizeX SizeY
            DefPar = [1.5, 1.5, 0, 11, 11];
            Nvar   = numel(varargin);
            Par    = [varargin{:}, DefPar(Nvar+1:end)];
            SigmaX = Par(1);
            SigmaY = Par(2);
            Rho    = Par(3);
            SizeX  = Par(4);
            SizeY  = Par(5);

            [MatX,MatY] = meshgrid((1:1:SizeX),(1:1:SizeY));
            X0  = ceil(SizeX.*0.5+0.1);
            Y0  = ceil(SizeY.*0.5+0.1);
            MatX = MatX - X0;
            MatY = MatY - Y0;

            K = 1./(2.*pi.*SigmaX.*SigmaY.*sqrt(1-Rho.^2)) .* ...
                                   exp(-1./(2.*(1-Rho.^2)) .* ...
                                       (MatX.^2./SigmaX.^2 + ...
                                        MatY.^2./SigmaY.^2 - ...
                            2.*Rho.*MatX.*MatY./(SigmaX.*SigmaY)));


        end % function gauss
        
        % exponential
        function K=exp(varargin)
            % Generate an exponential radial kernel with zero background.
            % Package: @Kernel2
            % Description: Generate an exponential radial kernel with zero background 
            %              and normalized to unity sum.
            % Input  : - Exponential scale radius. Default is 3 pix.
            %          - Stamp size in X direction. Default is 21.
            %          - Stamp Size in Y direction. Default is 21.
            % Output : - An exponential kernel located at the stamp center.
            %            The stamp center is defined to be the middle of the central
            %            pixel (i.e., delta function span over one pixel), where the
            %            central pixel is give by ceil(Size./2 + 0.1).
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Jun 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: K=Kernel2.exp;
            % Reliable: 2
            %--------------------------------------------------------------------------

            %         Radius SizeX SizeY
            DefPar = [3, 21, 21];
            Nvar   = numel(varargin);
            Par    = [varargin{:}, DefPar(Nvar+1:end)];
            Radius = Par(1);
            SizeX  = Par(2);
            SizeY  = Par(3);

            [MatX,MatY] = meshgrid((1:1:SizeX),(1:1:SizeY));
            %X0   = SizeX.*0.5;
            %Y0   = SizeY.*0.5;
            X0  = ceil(SizeX.*0.5+0.1);
            Y0  = ceil(SizeY.*0.5+0.1);

            MatX = MatX - X0;
            MatY = MatY - Y0;
            MatR2=MatX.^2 + MatY.^2;

            K = exp(-sqrt(MatR2)./Radius);
            K = K./sum(K(:));

        end % function exp

        % sersic
        function K=sersic(varargin)
            % Generate a Sersic radial kernel with zero background.
            % Package: @Kernel2
            % Description: Generate a Sersic radial kernel with zero background 
            %              and normalized to unity sum.
            %              Functional form: e^(-k*r^(1/n))
            % Input  : - Sersic k parameter. Default is 1.
            %          - sersic n parameter. Default is 4 (de Vaucouleurs)
            %          - Stamp size in X direction. Default is 21.
            %          - Stamp Size in Y direction. Default is 21.
            % Output : - A sersic kernel located at the stamp center.
            %            The stamp center is defined to be the middle of the central
            %            pixel (i.e., delta function span over one pixel), where the
            %            central pixel is give by ceil(Size./2 + 0.1).
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Sep 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: K=Kernel2.sersic;
            % Reliable: 2
            %--------------------------------------------------------------------------

            %         Radius SizeX SizeY
            DefPar = [1, 4, 21, 21];
            Nvar   = numel(varargin);
            Par    = [varargin{:}, DefPar(Nvar+1:end)];
            SerK   = Par(1);
            SerN   = Par(2);
            SizeX  = Par(3);
            SizeY  = Par(4);

            [MatX,MatY] = meshgrid((1:1:SizeX),(1:1:SizeY));
            %X0   = SizeX.*0.5;
            %Y0   = SizeY.*0.5;
            X0  = ceil(SizeX.*0.5+0.1);
            Y0  = ceil(SizeY.*0.5+0.1);

            MatX = MatX - X0;
            MatY = MatY - Y0;
            MatR2=MatX.^2 + MatY.^2;

            K = exp(-SerK.*sqrt(MatR2).^(1./SerN));
            K = K./sum(K(:));

        end % function sersic

        % aperture (top hat)
        function K=aper(varargin)
            % Generate an aperture (top-hat) circular kernel.
            % Package: @Kernel2
            % Description: Generate an aperture (top-hat) circular kernel with zero
            %              background and normalized to unity sum.
            % Input  : - Aperture radius. Default is 5 pix.
            %          - Stamp size in X direction. Default is 11.
            %          - Stamp Size in Y direction. Default is 11.
            % Output : - A Gaussian kernel located at the stamp center.
            %            The stamp center is defined to be the middle of the central
            %            pixel (i.e., delta function span over one pixel), where the
            %            central pixel is give by ceil(Size./2 + 0.1).
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Jun 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: K=Kernel2.aper;
            % Reliable: 2
            %--------------------------------------------------------------------------

            %         Radius SizeX SizeY
            DefPar = [5, 11, 11];
            Nvar   = numel(varargin);
            Par    = [varargin{:}, DefPar(Nvar+1:end)];
            Radius = Par(1);
            SizeX  = Par(2);
            SizeY  = Par(3);

            [MatX,MatY] = meshgrid((1:1:SizeX),(1:1:SizeY));
            %X0   = SizeX.*0.5;
            %Y0   = SizeY.*0.5;
            X0  = ceil(SizeX.*0.5+0.1);
            Y0  = ceil(SizeY.*0.5+0.1);


            MatX = MatX - X0;
            MatY = MatY - Y0;
            MatR2=MatX.^2 + MatY.^2;

            K = zeros(SizeY,SizeX);
            K(MatR2<=Radius.^2) = 1;
            K = K./sum(K(:));

        end % function aper

        % annulus
        function K=annulus(varargin)
            % Generate an annulus circular kernel.
            % Package: @Kernel2
            % Description: Generate an annulus circular kernel with zero
            %              background and normalized to unity sum.
            % Input  : - Annulus inner radius. Default is 5 pix.
            %          - Annulus outer radius. Default is 7 pix.
            %          - Stamp size in X direction. Default is 15.
            %          - Stamp Size in Y direction. Default is 15.
            % Output : - A Gaussian kernel located at the stamp center.
            %            The stamp center is defined to be the middle of the central
            %            pixel (i.e., delta function span over one pixel), where the
            %            central pixel is give by ceil(Size./2 + 0.1).
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Jun 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: K = Kernel2.annulus;
            % Reliable: 2
            %--------------------------------------------------------------------------


            %         InnerRadius OuterRadius SizeX SizeY
            DefPar = [5, 7 15, 15];
            Nvar   = numel(varargin);
            Par    = [varargin{:}, DefPar(Nvar+1:end)];
            InRad  = Par(1);
            OutRad = Par(2);
            SizeX  = Par(3);
            SizeY  = Par(4);

            [MatX,MatY] = meshgrid((1:1:SizeX),(1:1:SizeY));
            %X0   = SizeX.*0.5;
            %Y0   = SizeY.*0.5;
            X0  = ceil(SizeX.*0.5+0.1);
            Y0  = ceil(SizeY.*0.5+0.1);

            MatX = MatX - X0;
            MatY = MatY - Y0;
            MatR2=MatX.^2 + MatY.^2;

            K = zeros(SizeY,SizeX);
            K(MatR2>InRad.^2 & MatR2<=OutRad.^2) = 1;
            K = K./sum(K(:));

        end % function annulus
        
        % line
        function K=line(varargin)
            % Line kernel
            % Description: Generate a line kernel.
            % Input  : - Line length. Default is 11.
            %          - Line Width. Default is 1.
            %          - Line angle [deg]. Default is 0.
            %            Measured relative to the X axis in the
            %            trigonometric direction.
            %          - Size of gap in line center. Default is 0.
            %          - Stamp size in X direction. Default is line length.
            %          - Stamp size in Y direction. Default is line length.
            % Output : The kernel
            % Example: K=Kernel2.line
            % Reliable : 2
            
            
            DefPar = [11, 1, 0, 0, NaN, NaN];
            Nvar   = numel(varargin);
            Par    = [varargin{:}, DefPar(Nvar+1:end)];
            Length = Par(1);
            Width  = Par(2);
            Angle  = Par(3);
            Gap    = Par(4);
            SizeX  = Par(5);
            SizeY  = Par(6);
            
            % set angle to the [-90..+90] range
            Angle = mod(Angle,180);
            if (Angle>90)
                Angle = Angle - 180;
            end
            
            
            if (isnan(SizeX))
                SizeX = Length;
            end
            if (isnan(SizeY))
                SizeY = Length;
            end
            
            K = zeros(SizeY,SizeX);
            % Center
            CX = SizeX.*0.5 + 0.5;
            CY = SizeY.*0.5 + 0.5;
            
            X  = (1:1:SizeX);
            if (abs(Angle)>45)
                Angle = Angle - sign(Angle).*90;
                Rot   = true;
            else
                Rot = false;
            end
            Y  = CY + tand(Angle).*(X-CX);
            for k=-round(Width/2/cosd(Angle)):1:round(Width/2/cosd(Angle))
            Ind = sub2ind([SizeY, SizeX],round(Y)+k,round(X));
            K(Ind) = 1;
            end
            
            [MatX,MatY] = meshgrid((1:1:SizeX),(1:1:SizeY));
           
            MatX = MatX - CX;
            MatY = MatY - CY;
            MatR2=MatX.^2 + MatY.^2;
            K(MatR2>(Length.*0.5).^2)=0;
            if (Gap>0)
                K(MatR2<(Gap.*0.5).^2) = 0;
            end
            if (Rot)
                K = rot90(K);
            end
        
            K = K./sum(K(:));
            
        end
        
        % Poisson Matched filter
        function K=pmf(varargin)
            % Generate a Poisson-noise matched filter kernel
            % Package: @Kernel2
            % Description: Generate a Poisson-noise matched filter kernel
            %              from a PSF. The PSF is normalized, and the
            %              filter includes the background.
            % Input  : - Background [counts].
            %            Default is 0.005.
            %          - A PSF. If Scalar than a Gaussian sigma [pix].
            %            Default is 1.5.
            %          - Flux threshold (F_th). Default is 3 counts.
            %          - Stamp size in X direction. Default is 11.
            %          - Stamp Size in Y direction. Default is 11.
            % Output : - A Poisson=-noise kernel located at the stamp center.
            %            The stamp center is defined to be the middle of the central
            %            pixel (i.e., delta function span over one pixel), where the
            %            central pixel is give by ceil(Size./2 + 0.1).
            % See also: kernel_aper.m, kernel_annulus.m, kernel_exp.m
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    May 2017
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: K=Kernel2.pmf;
            %          K=Kernel2.pmf(0.005,Kernel2.gauss)
            % Reliable: 2
            %--------------------------------------------------------------------------

            %         Back   PSF  Fth SizeX SizeY
            DefPar = {0.005, 1.5, 3, 11, 11};
            Nvar   = numel(varargin);
            Par    = {varargin{:}, DefPar{Nvar+1:end}};
            Back   = Par{1};
            PSF    = Par{2};
            Fth    = Par{3};
            SizeX  = Par{4};
            SizeY  = Par{5};

            if (numel(PSF)==1)
                % PSF is a scalar
                [MatX,MatY] = meshgrid((1:1:SizeX),(1:1:SizeY));
                X0  = ceil(SizeX.*0.5+0.1);
                Y0  = ceil(SizeY.*0.5+0.1);
                MatX = MatX - X0;
                MatY = MatY - Y0;
                
                Rho = 0;
                P = 1./(2.*pi.*PSF.*PSF.*sqrt(1-Rho.^2)) .* ...
                                       exp(-1./(2.*(1-Rho.^2)) .* ...
                                           (MatX.^2./PSF.^2 + ...
                                            MatY.^2./PSF.^2 - ...
                                2.*Rho.*MatX.*MatY./(PSF.*PSF)));
            else
                % PSF is a matrix
                P = PSF;
                % normalize PSF
                P = P./sum(P(:));
            end
            
            K = log(1+ Fth.*P./Back);

            
        end
        
        
        % addition kernel
        
        % subtraction kernel
        function fftPd=psf_diff(Pr,Pn,Fr,Fn,SigmaR,SigmaN,FlagFFT)
            % Calculate the image subtraction PSF (P_D)
            % Description: Calculate the fft of the image subtraction PSF (P_D)
            %              of the ZOGY algorithm, given the ref and new
            %              PSFs.
            % Input  : - Pr: A matrix of the ref image PSF (or its fft).
            %          - Pn: A matrix of the new image PSF (or its fft).
            %          - F_r (flux ZP of ref image).
            %          - F_n (flux ZP of new image).
            %          - SigmaR (noise of ref image).
            %          - SigmaN (noise of new image).
            %          - A flag indicating if the input Pr, Pn PSFs are
            %            ffted (default is false).
            % Output : - The fft of the P_D.

            if (nargin==6)
                FlagFFT = false;
            end
            
            if (FlagFFT)
                fftPr = Pr;
                fftPn = Pn;
            else
                fftPr = fft2(Pr);
                fftPn = fft2(Pn);
            end
            SnFr2 = (SigmaN.*Fr).^2;
            SrFn2 = (SigmaR.*Fn).^2;
            Fd = Fr.*Fn./sqrt(SnFr2 + SrFn2);
            fftPd = Fr.*Fn.*fftPr.*fftPn./(Fd.*sqrt(SnFr2.*fftPr.*conj(fftPr) + SrFn2.*fftPn.*conj(fftPn)));
            
            
        end
        
        % motion/pwdvariability subtraction kernel
        function fftPdShift=psf_diff_shift(fftPd,Fr,Fn,SigmaR,SigmaN,Shift,AlphaR,AlphaN,varargin)
            % Calculate the image subtraction shifted/variable PSF
            % Description: Calculate the image subtraction shifted/variable PSF
            %              of the ZOGY algorithm
            %              fft(P_D)*(alpha_R + alpha_N fft(Shift),
            %              where Shift is the shift operator
            % Input  : - Pd: A matrix of the fft of thesubtraction image
            %            PSF.
            %          - F_r (flux ZP of ref image).
            %          - F_n (flux ZP of new image).
            %          - SigmaR (noise of ref image).
            %          - SigmaN (noise of new image).
            %          - Shift [X,Y] of the new relative to the ref.
            %          - AlphaR (flux zp of ref).
            %          - AlphaN (flux zp of new).
            %          * Additional parameters to pass to
            %            ImUtil.Im.imagefft_shift_fft.m ...,NY,NX,Nr,Nc.
            %            By default they will be calculated.
            % Output : - The fft of the shifted subtraction PSF.
            % Example: Pn = Kernel2.gauss(1.5,1.5); Pr = Kernel2.gauss(2,2);
            %       fftPd = Kernel2.psf_diff(Pr,Pn,1,1,0.01,0.01,false);
            %       fftPdShift=Kernel2.psf_diff_shift(fftPd,1,1,0.01,0.01,[1.4 0.3],1,1);  
            %       surface(fftshift(ifft2(fftPdShift)))
            
            warning('INCORRECT?!')
            
            %fftPd = Kernel2.psf_diff(Pr,Pn,Fr,Fn,SigmaR,SigmaN,FlagFFT);
            
            %[ShiftedPd,NY,NX,Nr,Nc]=AstroIm.imagefft_shift_fft(fftPd,Shift(1),Shift(2),varargin{:}); %NY,NX,Nr,Nc)
            [ShiftedPd]=ImUtil.Im.imagefft_shift_fft(fftPd,Shift(1),Shift(2),varargin{:}); %NY,NX,Nr,Nc)
        
            fftPdShift = AlphaR.*fftPd + AlphaN.*fft2(ShiftedPd);
            
        end    % Kernel2.psf_diff_shift
    end % static
    
    % Kernel utilities
    methods (Static)
        % kernel to cell of matrices
        function KernelCell=cellmat(Kernel,varargin)
            % Convert a set of kernels into a cell array of matrices.
            % Package: @Kernel2
            % Description: Convert a set of kernels into a cell array of matrices.
            % Input  : - A set of kernels. This can be:
            %            A ClassPSF object.
            %            A matrix.
            %            A function.
            %            A cell array of matrices.
            %            A cell array of functions.
            %          * Arbitrary number of arguments.
            %            If the first input argument is a ClassPSF object these
            %            arguments will be passed to getmpsf.m. If the first input
            %            argument is a function or cell of functions, then these
            %            arguments will be passed to the function.
            %            If this is a single cell array than it will be assumed that
            %            each element of the cell array is a cell array that will be
            %            passed to the corresponding function (see example).
            % Output : - A cell array of matrices. Each matrix represent the kenel.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Jun 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: KernelCell=Kernel2.cellmat(rand(11,11));
            %          KernelCell=Kernel2.cellmat({rand(11,11)});
            %          KernelCell=Kernel2.cellmat(@Kernel2.gauss);
            %          KernelCell=Kernel2.cellmat({@Kernel2.exp},3);
            %          KernelCell=Kernel2.cellmat({@Kernel2.exp,@Kernel2.gauss},{{3},{2,2}});
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (ClassPSF.isClassPSF(Kernel))
                % kernel is a ClassPSF object
                KernelCell = getmpsf(Kernel,varargin{:});
            elseif (isnumeric(Kernel))
                % Kernel is a matrix
                KernelCell = {Kernel};
            elseif (isa(Kernel,'function_handle'))
                % Kernel is a function
                KernelCell = {Kernel(varargin{:})};
            elseif (iscell(Kernel))
                Ncell      = numel(Kernel);
                KernelCell = cell(Ncell,1);
                for Icell=1:1:Ncell
                    if (isnumeric(Kernel{Icell}))
                        % Kernel cell element is a matrix
                        KernelCell{Icell} = Kernel{Icell};
                    elseif (isa(Kernel{Icell},'function_handle'))
                        % Kernel cell element is a function 
                        if (numel(varargin)>0)
                            if (iscell(varargin{1}))
                                Par = varargin{1}; %{Icell};
                            else
                                Par = varargin;
                            end
                        else
                            Par = {};
                        end
                        KernelCell{Icell} = Kernel{Icell}(Par{Icell}{:});
                    else
                        error('Unknown Kernel input option');
                    end
                end
            else
                error('Unknown Kernel input option');
            end

        
        end % function cellmat
        
    end % static
    
end % end class
            
