function [Sim,PSF]=image_art(varargin)
%--------------------------------------------------------------------------
% image_art function                                               ImBasic
% Description: Generate artificial images or add artificial sources
%              to real images.
% Input  : * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Image'- Image matrix or a single SIM image to which
%                     to add sources. If empty a new image will be created.
%                     Default is empty.
%            'ImSize' - Size of new image to create. Default is [256 256].
%            'Bias'   - Bias level to add. Default is 0.
%            'RN'     - Read noise (electron) to add.
%            'Gain'   - Gain. Default is 1.0.
%            'AddNoise'- Add noise to image. Default is true.
%            'Back'   - Background level [electrons]. Default is 300.
%            'PSF'    - PSF type. Options are {'gauss'}.
%                       Default is 'gauss'.
%            'ParPSF' - PSF parameters. Default is [1.5 1.5 0]
%                       i.e., sigmaX, sigmaY Rho for 'gauss'.
%            'MaxRad' - PSF maximum radius. Default is 10.
%            'StarList' - Table of [X, Y, Counts] per source.
%                       Default is [128 128 2000].
%            'OutFileName' - FITS output file name to write.
%                       Default is empty (do not write file).
%            'Nnoisereal' - Number of noise realizations.
%                        Default is 1.
%                        If >1 multiple images will be returned/written.
%                        If 0 than no noise will be added.
%            'Header' - FITS header to write. Default is [].
%            'DataType' - FITS datatype to write. Default is 'single'.
%            'OutSIM' - Return SIM images (true) or a matrix (false).
%                       Default is true.
% Output : - Output images in SIM or matrix format.
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Aug 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: OutImage=image_art;
%          [MatX,MatY] = meshgrid((30:30:1000)',(30:30:1000)'); N=numel(MatX);
%          List = [MatX(:)+rand(N,1), MatY(:)+rand(N,1), 100000.*ones(N,1)];
%          List = rand(200,3); List(:,1:2)=List(:,1:2).*1024;
%          List(:,3)=5000;
%          Out=image_art('ImSize',[1024 1024],'StarList',List,'OutFileName','N.fits');
% Reliable: 2
%--------------------------------------------------------------------------

ImageField       = 'Im';
FileNameField    = 'ImageFileName';


DefV.Image          = [];
DefV.ImSize         = [256 256]; %[1024 1024];
DefV.Bias           = 0;
DefV.RN             = 0;
DefV.Gain           = 1.0;
DefV.AddNoiseSource = true;
DefV.AddNoiseBack   = true;
DefV.AddNoiseRN     = true;
DefV.Back           = 300;       % e-
DefV.PSF            = 'gauss';
DefV.ParPSF         = [1.5 1.5 0];   % sigmaX, sigmaY Rho for 'gauss'
DefV.MaxRad         = 10;
DefV.StarList       = [128 128 2000];
DefV.OutFileName    = [];
DefV.Nnoisereal     = 1;         % number of noise realizations - if 0 no noise added
%DefV.PSFhalfSize    = 20;
DefV.Header         = [];
DefV.DataType       = 'single';
DefV.OutSIM         = true;
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% create image
if (SIM.issim(InPar.Image)),
    OutImage = InPar.Image.(ImageField);
else
    OutImage = InPar.Image;
end
if (isempty(OutImage)),
    OutImage = zeros(InPar.ImSize);
end
% override ImSize if input image is provided
InPar.ImSize = size(OutImage);

% add bias
OutImage = OutImage + InPar.Bias;


% coordinates
VecX = (1:1:InPar.ImSize(2)).';
VecY = (1:1:InPar.ImSize(1)).';
[MatX,MatY] = meshgrid(VecX,VecY);

% stars
Nstar = size(InPar.StarList,1);
switch lower(InPar.PSF)
    case 'gauss'
        for Istar=1:1:Nstar,
            %Flux = poissrnd(InPar.StarList(Istar,3));
            %   [X0, Y0, SigmaX, SigmaY, Rho, Norm, Const],
            OutImage = OutImage + Util.fun.bivar_gauss(MatX,MatY,...
                                              [InPar.StarList(Istar,1:2), InPar.ParPSF(1:3), InPar.StarList(Istar,3), 0]);
                                            
        end
    otherwise
        error('Unknown PSF option');
end

% if (~InPar.AddNoise),
%     OutImageN = OutImage;
% end

if (InPar.OutSIM),
    Sim = SIM;
end

for Inr=1:1:InPar.Nnoisereal,
    
    % add noise
    % RN and background
    OutImageN = zeros(size(OutImage));
    if (InPar.AddNoiseRN),
        OutImageN = OutImageN + randn(InPar.ImSize).*InPar.RN;
    end
        
    if (InPar.AddNoiseSource),
        OutImageN = OutImageN + poissrnd(OutImage,InPar.ImSize);
    else
        OutImageN = OutImageN + OutImage;
    end

    if (InPar.AddNoiseBack),
        OutImageN = OutImageN + poissrnd(InPar.Back,InPar.ImSize);
    else
        OutImageN = OutImageN + InPar.Back;
    end
         
    % divide by gain
    OutImageN = OutImageN./InPar.Gain;


    if (~isempty(InPar.OutFileName)),
        ImageName = sprintf('%s%03d',InPar.OutFileName,Inr);
        FITS.write_old(OutImageN,ImageName,InPar.Header,InPar.DataType);
    else
        ImageName = [];
    end
    if (InPar.OutSIM),
        Sim(Inr).(ImageField)    = OutImageN;
        Sim(Inr).(FileNameField) = ImageName;
    else
        Sim = OutImageN;
    end
end

if (InPar.Nnoisereal==0),
    if (~isempty(InPar.OutFileName)),
        ImageName = sprintf('%s',InPar.OutFileName);
        FITS.write_old(OutImage,ImageName,InPar.Header,InPar.DataType);
    else
        ImageName = [];
    end
    if (InPar.OutSIM),
        Sim.(ImageField)    = OutImage;
        Sim.(FileNameField) = ImageName;
    else
        Sim = OutImage;
    end
end


if (nargout>1),
    %VecX = (-InPar.PSFhalfSize:1:InPar.PSFhalfSize)';
    %VecY = VecX;
    if (is_evenint(InPar.ImSize(2))),
        VecX = (0.5-InPar.ImSize(2).*0.5:1:InPar.ImSize(2).*0.5-0.5).';
    else
        VecX = (-floor(InPar.ImSize(2).*0.5):1:floor(InPar.ImSize(2).*0.5))';
    end
    if (is_evenint(InPar.ImSize(1))),
        VecY = (0.5-InPar.ImSize(1).*0.5:1:InPar.ImSize(1).*0.5-0.5).';
    else
        VecY = (-floor(InPar.ImSize(1).*0.5):1:floor(InPar.ImSize(1).*0.5))';
    end
    [MatX,MatY] = meshgrid(VecX,VecY);
    %   [X0, Y0, SigmaX, SigmaY, Rho, Norm, Const],
    PSF.Im = bivar_gauss(MatX,MatY,...
                         [0, 0, InPar.ParPSF(1:3), 1, 0]);
end
