function M=mask_init(M,Type,Size)
% Initizalize a MASK object.
% Package: @MASK
% Description: Initizalize a MASK object.
% Input  : - A MASK object.
%          - Data type for mask object. Default is 'uint32'.
%          - Size of Mask image. 
%            If the input object is a SIM object and requested mask size is
%            [0 0], then will set the size to the SIM image size.
%            Default is [0 0].
% Output : - An initilazied Mask object.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: M=mask_init(M);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Type = 'uint32';
Def.Size = [0 0];
if (nargin==1)
    Type   = Def.Type;
    Size   = Def.Size;
elseif (nargin==2)
    Size   = Def.Size;
elseif (nargin==3)
    % do nothing
else
    error('Illegal number of input arguments');
end

MaskField    = 'Mask';
ImageField   = 'Im';

Nm = numel(M);
for Im=1:1:Nm
    if (SIM.issim(M) && sum(Size)==0)
        SizeIm = size(M(Im).(ImageField));
    else
        SizeIm = Size;
    end
    M(Im).(MaskField) = zeros(SizeIm,Type);
end
