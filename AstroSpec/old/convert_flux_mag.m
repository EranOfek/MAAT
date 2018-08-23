function Output=convert_flux_mag(Input,Dir,ZeroPoint);
%--------------------------------------------------------------------------
% convert_flux_mag function                                      AstroSpec
% Description: Convert flux (and errors) units to magnitudes and
%              vise-versa. Also handles nJy, AB mag, ergs/cm^2/s/Hz.
% Input  : - magnitude or flux and optional error in second column.
%          - Direction:
%            'f'   - magnitude to flux.
%            'm'   - flux to magnitude.
%            'J'   - AB magnitude to nJy (zero-point=31.40).
%            'A'   - nJy to AB magnitude (zero-point=31.40).
%            'E'   - AB magnitude to ergs/cm^2/s/Hz (zero-point=-48.60).
%            'e'   - ergs/cm^2/s/Hz to AB magnitude (zero-point=-48.60).
%          - zero-point, (if 'f' | 'm' default is 22).
% Output : - flux or magnitude (and their errors, in the second column).
% Tested : Matlab 5.3
%     By : Eran O. Ofek                     April 2002
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------------
ColM             = 1;
ColE             = 2;
ZeroPointDef     = 22;
if (nargin==2),
   ZeroPoint = ZeroPointDef;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end

[Nr,Nc] = size(Input);
if (Nc==1),
   Input = [Input, zeros(Nr,1)];
end

switch Dir
 case {'J','A'}
    ZeroPoint = 31.40;
 case {'E','e'}
    ZeroPoint = -48.60;
 otherwise
    % do nothing
end

switch Dir
case {'f','J','E'}
    %--- Convert magnitude to flux ---
    Mag      = Input(:,ColM);
    MagErr   = Input(:,ColE);
    Flux     = 10.^(0.4.*(ZeroPoint - Mag));
    FluxErr  = abs(-log(10).*0.4.*Flux.*MagErr);

    Output   = [Flux, FluxErr];

 case {'m','A','e'}
    %--- Convert flux to magnitude ---
    Flux     = Input(:,ColM);
    FluxErr  = Input(:,ColE);
    Mag      = ZeroPoint - 2.5.*log10(Flux);
    MagErr   =  (2.5./(log(10).*Flux)).*FluxErr;

    Output   = [Mag, MagErr];


 otherwise
    error('Unknown Dir parameter');
end
