function ZP=zp(Band)
% Get the GALEX photometric zeropoint.
% Package: VO.GALEX
% Description: Return GALEX zeropoint for conversion of counts/second
%              to AB mag.
% Input  : - Band: 'NUV'|'n' | 'FUV'|'f'. Default is 'NUV'.
% Output : - Photometric zero point: AB mag which produce 1
%            count per second.
% Example: VO.GALEX.zp
% Reliable: 1

if (nargin==0)
    Band = 'NUV';
end
switch lower(Band)
    case {'nuv','n'}
        ZP = 20.08;
    case {'fuv','f'}
        ZP = 18.82;
    otherwise
        error('Unknown band option');
end

