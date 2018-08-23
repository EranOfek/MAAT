function InterpImage=iminterp(Image,Mask,varargin)
% Interpolate a 2D image over NaNs or mask image
% Package: ImUtil.Im
% Description: In a 2-D image interpolate over NaNs or over pixels in
%              which a bit mask is set.
% Input  : - Image.
%          - Mask Image. This is either a bit mask, or an array of
%            logicals (see 'Bit' for details). Default is empty.
%            Will interpolate in any pixels >0.
%            If empty, then will interpolate over NaN pixels.
%          * Arbitrary number of pairs of input arguments ...,key,val,...
%            The following keywords are available:
%            'Bit' - The index of a bit in the bit mask (second input
%                    argument) which to select pixels over which
%                    interpolation is required.
%                    Default is empty. If empty then assume that the Mask
%                    is a 2-D array of logicals set to true in pixels
%                    over which interpolation is required.
%            'IntMethod'- inpaint_nans.m interpolation method. Default is 0.
% Output : - Interpolated image
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [MatX,MatY]=meshgrid((1:1:5),(1:1:5)); A=MatX.^2+sqrt(MatY);
%          A(3,3)=NaN; A(2,3)=NaN;
%          InterpImage=ImUtil.Im.iminterp(A);
%          [MatX,MatY]=meshgrid((1:1:5),(1:1:5)); A=MatX.^2+sqrt(MatY);
%          Flag=false(5,5); A(3,3)=true;
%          InterpImage=ImUtil.Im.iminterp(A,Flag);
%          BitMask=zeros(5,5,'uint16'); BitMask(3,4)=1;
%          InterpImage=ImUtil.Im.iminterp(A,BitMask,'Bit',1);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1)
    Mask = [];
end

DefV.Bit       = [];
DefV.IntMethod = 0;
%InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (isempty(InPar.Bit))
    % assume Mask is a flag image
else
    Mask = logical(bitget(Mask,InPar.Bit));
end

if (isempty(Mask))
    % interpolate over NaN
    ImageNaN = Image;
else
    ImageNaN       = Image;
    ImageNaN(Mask>0) = NaN;
end

InterpImage = Util.external.inpaint_nans(ImageNaN,InPar.IntMethod);
