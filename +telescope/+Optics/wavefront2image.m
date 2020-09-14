function [Image_NS,Image]=wavefront2image(SumY,Support,OverSamp,Norm,Pars)
% Construct an image from the wavefront
% Package: telescope.Optics
% Description: Construct a wavefront from a matrix representation of the
%              wavefront, and calculate the image plane of a point source
%              with a given support (obstraction). 
% Input  : - Matrix which represent the amplitude of the wavefront in each
%            point in the pupil plane. The pupil plane range is always
%            between -1 and 1 in both axes.
%            For example, this is the SumY matrix returned by
%            zerwavefront2image.m.
%            Default is ones(128,128).
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
%            'no'  - do not normalize.
%            'sum' - Normalize by the final image sum. Default. 
%          - Optional supports parameters. See telescope_support.m
%            Default is 0.15.
% Output : - Image plane before the fft shift.
%          - Image plane after the fft shift (source is in the center).
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: zerwavefront2image_indiv.m, zerwavefront2image.m
% Example: % generate Airy function
%          [Image_NS,Image]=wavefront2image([],'circ',4);
%          pcolor(log10(Image)), shading interp; axis square, colorbar
% Reliable: 2
%--------------------------------------------------------------------------

open telescope.Optics.*

Def.SumY     = ones(128,128);
Def.Support  = 'circ';
Def.OverSamp = 2;  % in lambda/d units...
Def.Norm     = 'sum';
Def.Pars     = 0.15;

if (nargin==0)
    SumY    = Def.SumY;
    Support = Def.Support;
    OverSamp= Def.OverSamp;
    Norm    = Def.Norm;
    Pars    = Def.Pars;
elseif (nargin==1)
    Support = Def.Support;
    OverSamp= Def.OverSamp;   
    Norm    = Def.Norm;
    Pars    = Def.Pars;
elseif (nargin==2)
    OverSamp= Def.OverSamp;
    Norm    = Def.Norm;
    Pars    = Def.Pars;
elseif (nargin==3)
    Norm    = Def.Norm;
    Pars    = Def.Pars;
elseif (nargin==4)
    Pars    = Def.Pars;
elseif (nargin==5)    
    % do nothing
else
    error('Illegal number of input arguments');
end

if (isempty(SumY))
    SumY = Def.SumY;
end
  

Size = size(SumY);

if (ischar(Support))
    [Support,VecX,VecY]=telescope_support(Size,Support,Pars);
end


% Calculate the wavefront summed over all Zernike functions:
%SumY = sum(bsxfun(@times, zernike_xy(n,m,VecX,VecY,true) ,reshape(Amp,1,1,length(Amp))), 3);
%SumY = sum(bsxfun(@times, zernike_xy(n,m,MatX,MatY,true) ,reshape(Amp,1,1,length(Amp))), 3);

SumY(isnan(SumY)) = 0;
% no shift...
FFT_NS   = fft2(Support.*exp(-1i.*SumY),OverSamp.*Size(1),OverSamp.*Size(2) );
Image_NS = abs(FFT_NS).^2;
%Image = abs(fftshift(fftshift(fft2(Support.*exp(-i.*SumY),OverSamp.*Size(1),OverSamp.*Size(2) ), 1),2)).^2;
if (nargout>1),
    Image    = abs(fftshift(fftshift(FFT_NS, 1),2)).^2;
end

switch lower(Norm)
    case 'no'
        % do nothing
    case 'sum'
        Norm     = sum(Image_NS(:));
        Image_NS = Image_NS./Norm;
        if (nargout>1),
            Image    = Image./Norm;
        end
    otherwise
        error('Unknown Norm option');
end
