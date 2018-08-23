function [Support,VecX,VecY]=telescope_support(varargin)
% Construct an image of a telescope entrance pupil (support).
% Package: telescope.Optics
% Description: Construct an image of a telescope entrance pupil (support).
%              For example, a clear circular pupil.
% Input  : - Size of matrix to construct (e.g., [128 128]).
%            Default is [128 128].
%          - One of the following support types:
%            'circ'    - Circular support (default).
%            'circobs' - Circular support with central obscuration.
%            'rect'    - Rectangular support.
%          - Additional parameters describing the support.
%            For 'circ' - No parameters.
%            For 'circobs' - A parameter describing the radius of the
%                            central obscuration in units of the
%                            aperture. 
%            For 'rect' - The ratio between the axis.   
%            Default is 0.15.
% Output : - A matrix of the support.
%          - Vector of X-axis grid points.
%          - Vector of Y-axis grid points.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Support,VecX,VecY]=telescope_support([128 128]);
% Reliable: 2
%--------------------------------------------------------------------------

NumVarargs = length(varargin);
if NumVarargs > 3
     errId = 'telescope_support:TooManyInputArguments';
     errMsg = 'Size, [Support, Pars]';
     error(errId, errMsg);
end
Gaps = cellfun(@isempty, varargin);
DefArgs = {[128 128] 'circ' 0.15};    % default input arguments
Suboptargs = DefArgs(1 : NumVarargs);
varargin(Gaps) = Suboptargs(Gaps);
DefArgs(1 : NumVarargs) = varargin;
[Size, Support, Pars] = DefArgs{:};

% 
% Def.Support = 'circ'; %{'circ'|'circobs'}
% Def.Pars    = 0.15;
% if (nargin==1),
%     Support = Def.Support;
%     Pars    = Def.Pars;
% elseif (nargin==2),
%     Pars    = Def.Pars;
% elseif (nargin==3),
%     % do nothing
% else
%     error('Illegal number of input arguments');
% end



if (ischar(Support)),
    switch lower(Support)
        case 'circ'
            VecX = (-1:2./(Size(2)-1):1).';
            VecY = (-1:2./(Size(1)-1):1).';

            [MatX,MatY] = meshgrid(VecX,VecY);

            MatR        = sqrt(MatX.^2 + MatY.^2);
            Support     = ones(size(MatR));
            Support(MatR>1) = 0;
        case 'circobs'
            VecX = (-1:2./(Size(2)-1):1).';
            VecY = (-1:2./(Size(1)-1):1).';

            [MatX,MatY] = meshgrid(VecX,VecY);

            MatR        = sqrt(MatX.^2 + MatY.^2);
            Support     = ones(size(MatR));
            Support(MatR>1 | MatR<Pars) = 0;
        case 'rect'
            VecX = (-1:2./(Size(2)-1):1).';
            VecY = (-1:2./(Size(1)-1):1).';

            [MatX,MatY] = meshgrid(VecX,VecY);
            Support     = ones(size(MatX));
            Support(abs(MatY)>Pars) = 0;
        otherwise
            error('Unknown Support option');
    end
else
    error('Support msut be a char array');
end
