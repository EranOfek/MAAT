function [Sim,Psf,CorrFactor]=filter(Sim,Filter,varargin)
% Apply a 2-D filter to all the images in a SIM object.
% Package: @SIM
% Description: Apply a 2-D filter to all the images in a SIM object.
%              The filtering process is done by either:
%              filter2.m, filter2_fft.m, or filter2_fftc.m.
%              Also update the noise in the 'ErrIm' field to account for
%              the filtering process.
% Input  : - A SIM object.
%          - An optional kernel by which to filter the image.
%            This can be a matrix, a cell array of matrices, a function
%            handle (parameters via 'FilterFunPar'), or a ClassPSF object
%            that will be read using getmpsf (parameters via 'GetPsfPar').
%            Kernel is expected to be centered.
%            If empty then will attempt to read the kernel from the
%            ClassPSF object within the SIM object. Default is empty.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FilterFun' - One of the following filtering function to use:
%                          'auto' - select automatically.
%                          'filter2' - Use filter2.m.
%                          'filter2_fft' - Use filter2_fft.m.
%                          'filter2_fftc' - Use filter2_fftc.m.
%            'ExecField' - Cell array of SIM fields on which to execute the
%                          filter. Default is {'Im'}.
%            'GetPsfPar' - A cell array of arguments to pass to the
%                          ClassPSF/getmpsf.m function.
%                          Default is {}.
%            'FilterFunPar' - A cell array of arguments to pass to the
%                          to the Filter argument function handle.
%                          Default is {1.5,1.5,0,15,15}.
%                          The default is sutiable for the @Kernel2.gauss
%                          function.
%            'UpdateErrIm' - A flag indicating if to correct the 'ErrIm'
%                          field noise to represent the noise of the
%                          filtered image. The correction factor is
%                          sqrt(sum(Psf.^2)). Default is true.
%            'NormErr'   - Normalize the filtered image by the ErrIm field,
%                          so that the filtered image will be in units of
%                          sigmas. Default is false.
%                          Will be applied only if UpdateErrIm is true.
% Output : - A SIM object with the filtered images.
%          - A cell array of the PSFs (kernel filters) used.
%          - A vector of the noise correction factors due to the filtering
%            (if applicable).
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=filter(Sim,@Kernel2.gauss);
%          Sim=filter(Sim,Kernel2.gauss(2,2));
% Reliable: 2
%--------------------------------------------------------------------------

ThresholdFilter = 0.05;  % The threshold by which to choose which algorithm to use

ImageField     = 'Im';
ErrField       = 'ErrIm';

if (nargin==1)
    % Default Filter is to use the PSF object in the SIM
    Filter = [];
end


DefV.FilterFun          = 'auto';   % 'auto' | 'filter2' | 'filter2_fft' | 'filter2_fftc'
DefV.ExecField          = {ImageField};
DefV.GetPsfPar          = {};  % additional parameters to pass to getmpsf.m
DefV.FilterFunPar       = {1.5,1.5,0,15,15};   % appropriate for Filter = @Kernel2.gauss
DefV.UpdateErrIm        = true;
DefV.NormErr            = false;
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end
Nf = numel(InPar.ExecField);

% get the PSF object from the SIM object
% The Psf is returned a cell array of PSFs
if (isempty(Filter))
    % no user supplied Filter
    % attempt to read the PSF object from SIM
    Psf = getmpsf(Sim,InPar.GetPsfPar{:});
    if any(Util.cell.isempty_cell(Psf))
        error('Some of the PSF objects doesnt exist');
    end
else
    % Use the user supplied Filter (in the second input argument)
    if (isnumeric(Filter))
        % Filter is numeric
        Psf{1} = Filter;
    elseif (isa(Filter,'function_handle'))
        % Filter is a function handle
        Psf{1} = Filter(InPar.FilterFunPar{:});
    elseif (iscell(Filter))
        % assume Filter is a cell array of numerical PSFs
        Psf = Filter;
    elseif (ClassPSF.isClassPSF(Filter))
        Psf = getmpsf(Filter,InPar.GetPsfPar{:});
    else
        error('Illegal Filter input');
    end
end

% check if all the images have the same size:
% going to assume that all the PSFs have the same size
SizeIm = imagesize(Sim);
if (all(SizeIm(1,1)==SizeIm(:,1)) && all(SizeIm(1,1)==SizeIm(:,1)))
    SameSize = true;
else
    SameSize = false;
end


Npsf = numel(Psf);
Nsim = numel(Sim);
CorrFactor = zeros(Nsim,1);
for Isim=1:1:Nsim
    % for each SIM object element
    Ipsf    = min(Isim,Npsf);
    
    % Get the filtering function
    % nominally, this will be executed once - unless requested or input
    % images have different sizes.
    % Note: PSF with variable stamp size are not taken into account when
    % estimating the faster algorithm (i.e., the 'auto' option).
    if (~SameSize || Isim==1)
        SizePsf = size(Psf{Ipsf});
        SizeIm  = size(Sim(Isim).(ImageField));

        % Select the filtering function
        % either: filter2, filter2_fft, filter2_fftc
        switch lower(InPar.FilterFun)
            case 'auto'
                if (all(SizeIm==SizePsf))
                    % Image and PSF have the same size
                    FilterFun = @filter2_fftc;
                else
                    if (max(SizePsf./SizeIm)<ThresholdFilter)
                        FilterFun = @filter2;
                    else
                        FilterFun = @Util.filter.filter2_fft;
                    end
                end
            case 'filter2'
                FilterFun = @filter2;
            case 'filter2_fft'
                FilterFun = @Util.filter.filter2_fft;
            case 'filter2_fftc'
                FilterFun = @filter2_fftc;
            otherwise
                error('Unknown FilterFun option');
        end
    end
    
    % Filter the images
    for If=1:1:Nf
        % for each SIM field
        switch func2str(FilterFun)
            case 'filter2'
                Sim(Isim).(InPar.ExecField{If}) = FilterFun(Psf{Ipsf},Sim(Isim).(InPar.ExecField{If}),'same');
            otherwise
                Sim(Isim).(InPar.ExecField{If}) = FilterFun(Sim(Isim).(InPar.ExecField{If}),Psf{Ipsf});
        end
    end
    
    % Update the ErrIm field that it will represent the noise of the
    % filtered image
    if (InPar.UpdateErrIm && ~isempty(Sim(Isim).(ErrField)))
        % The noise correction factor is given by:
        CorrFactor(Isim) = sqrt(sum(Psf{Ipsf}(:).^2));
        Sim(Isim).(ErrField) = Sim(Isim).(ErrField).*CorrFactor(Isim);
        
        if (InPar.NormErr)
            Sim(Isim).(InPar.ExecField{If}) = Sim(Isim).(InPar.ExecField{If})./Sim(Isim).(ErrField);
        end
    end
    
    
        
    
end

        
                        
                        
           
