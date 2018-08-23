function [Image_NS,Image,Y]=zerwavefront2image_indiv(n,m,Amp,VecX,VecY,Support,OverSamp,Norm)
% Wavefront from individual Zernike polynomials and calculate the image.
% Package: telescope.Optics
% Description: Construct a wavefront using Zernike polynomials, and
%              calculate the image plane of a point source with a
%              given support (obstraction).
%              However, unlike wavefront.m this function works on
%              individual Zernike polynomials and return a cube in
%              which each the the third dimension corresponds to
%              the image plane of an individual Zernike function.
% Input  : - Row vector of Zernike n indices. Alternatively, if the
%            second input argument is empty, then this is the Noll index
%            (see noll_index.m).
%            Default is (1:1:10).
%          - Row vector of Zernike m indices. Default is empty.
%          - Amplitude corresponding to each [n,m] value.
%            Default is ones(size(n)).
%          - Vector of X in the X-Y grid. Default is (-1:1./63.5:1).'.
%            If scalar then this is the number of elements in the vector.
%          - Vector of Y in the X-Y grid. Default is (-1:1./63.5:1).'.
%            If scalar then this is the number of elements in the vector.
%          - Telescope support. This is a matrix of the telescope
%            support (obscurations). The size of this matrix should be
%            identical to the size of the X-Y grid. Alternatively,
%            this can be one of the following strings:
%            'circ' - generate a circular support. Default.
%          - Over sampling. Default is 2 (i.e., Nyquist sampling).
%            Always provide sampling larger than 2.
%            The final image pixel scale is lambda/d/OverSamp.
%          - Normalize the final image plans (Image and Image_NS).
%            Options are:
%            'no'   - do not normalize.
%            'sum'  - Normalize by the final image cube sum. 
%            'sum1' - Normalize by the individual images sum. Default. 
% Output : - Image plane before the fft shift.
%          - Image plane after the fft shift.
%          - The sum of the Zernike polynomials weighted by their
%            amplitudes.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: zerwavefront2image.m, wavefront2image.m
% Example: [Image_NS,Image,SumY]=zerwavefront2image_indiv;
%          % To generate a realistic speckle image generated
%          % by the atmosphere:
%          % Note that the size of the seeing disk in this example is
%          (D/r0)*OverSamp:
%          [AmpC,J,C]=zer_cj_variance(100,'Nrand',1);
%          [Image_NS,Image,SumY]=zerwavefront2image_indiv(J,[],ones(1,100));
%          pcolor(log10(Image)), shading interp
%          axis square, colorbar
% Reliable: 2
%--------------------------------------------------------------------------

import telescope.Optics.*

Parallel=false;   % indicating if to calculate fft2 with loop or no-loop

Def.n        = (1:1:10);
Def.m        = [];
Def.Amp      = [];
Def.VecX     = (-1:1./63.5:1).';
Def.VecY     = (-1:1./63.5:1).';
Def.Support  = 'circ';
Def.OverSamp = 2;  % in lambda/d units...
Def.Norm     = 'sum1';

if (nargin==0),
    n       = Def.n;
    m       = Def.m;
    Amp     = Def.Amp;
    VecX    = Def.VecX;
    VecY    = Def.VecY;
    Support = Def.Support;
    OverSamp= Def.OverSamp;
    Norm    = Def.Norm;
elseif (nargin==1),
    m       = Def.m;
    Amp     = Def.Amp;
    VecX    = Def.VecX;
    VecY    = Def.VecY;
    Support = Def.Support;
    OverSamp= Def.OverSamp;   
    Norm    = Def.Norm;    
elseif (nargin==2),
    Amp     = Def.Amp;
    VecX    = Def.VecX;
    VecY    = Def.VecY;
    Support = Def.Support;
    OverSamp= Def.OverSamp;
    Norm    = Def.Norm;
elseif (nargin==3),
    VecX    = Def.VecX;
    VecY    = Def.VecY;
    Support = Def.Support;
    OverSamp= Def.OverSamp;
    Norm    = Def.Norm;
elseif (nargin==4),
    VecY    = Def.VecY;
    Support = Def.Support;
    OverSamp= Def.OverSamp;
    Norm    = Def.Norm;
elseif (nargin==5),
    Support = Def.Support;
    OverSamp= Def.OverSamp;
    Norm    = Def.Norm;
elseif (nargin==6),
    OverSamp= Def.OverSamp;
    Norm    = Def.Norm;
elseif (nargin==7),
    Norm    = Def.Norm;
elseif (nargin==8),
    % do nothing
else
    error('Illegal number of input arguments');
end


if (isempty(Amp)),
    Amp = ones(size(n));
end


if (isempty(VecX)),
    VecX = Def.VecX;
end
if (isempty(VecY)),
    VecY = Def.VecY;
end
if (isempty(Support)),
    Support = Def.Support;
end


if (numel(VecX)==1),
    % assume that VecX contains the number of steps
    VecX = (-1:2./(VecX-1):1).';
end
if (numel(VecY)==1),
    % assume that VecY contains the number of steps
    VecY = (-1:2./(VecY-1):1).';
end


if (ischar(Amp)),
    switch lower(Amp)
        case 'rand'
        otherwise
            
    end
end

[MatX,MatY] = meshgrid(VecX,VecY);

if (ischar(Support)),
    switch lower(Support)
        case 'circ'
            MatR        = sqrt(MatX.^2 + MatY.^2);
            Support     = ones(size(MatR));
            Support(MatR>1) = 0;
        otherwise
            error('Unknown Support option');
    end
end


% Calculate a cube of the wavefronts per Zernike functions:
NormN = sum(Support(:));   % true
Y = bsxfun(@times, zernike_xy(n,m,MatX,MatY,NormN) ,reshape(Amp,1,1,length(Amp)));

Size = size(Y);
Y(isnan(Y)) = 0;

% no shift...
if (Parallel)
    % This method is somewhat faster, but for some reason it have some
    % rounidng errors
    FFT_NS   = fft2(bsxfun(@times,Support, exp(-i.*Y)),OverSamp.*Size(1),OverSamp.*Size(2) );
else
    % slower but no rounding errors
    Tmp = bsxfun(@times,Support, exp(-1i.*Y));
    FFT_NS = zeros(Size(1).*OverSamp,Size(2).*OverSamp,Size(3));
    for I=1:1:Size(3),
        FFT_NS(:,:,I) = fft2(Tmp(:,:,I),OverSamp.*Size(1),OverSamp.*Size(2) );
    end
end
    
Image_NS = abs(FFT_NS).^2;


%Image = abs(fftshift(fftshift(fft2(Support.*exp(-i.*SumY),OverSamp.*Size(1),OverSamp.*Size(2) ), 1),2)).^2;
if (nargout>1),
   Image    = abs(fftshift(fftshift(FFT_NS, 1),2)).^2;
end

switch lower(Norm)
    case 'no'
        % do nothing
    case 'sum1'
        SizeIm    = size(Image_NS);
        Im2       = reshape(Image_NS,SizeIm(1).*SizeIm(2),SizeIm(3));
        Norm      = sum(Im2,1);
        Image_NS  = bsxfun(@times,Image_NS,reshape(1./Norm,1,1,SizeIm(3)));
        if (nargout>1),
            Image     = bsxfun(@times,Image,reshape(1./Norm,1,1,SizeIm(3)));
        end
    case 'sum'
        Norm     = sum(Image_NS(:));
        Image_NS = Image_NS./Norm;
        if (nargout>1),
            Image    = Image./Norm;
        end
    otherwise
        error('Unknown Norm option');
end
