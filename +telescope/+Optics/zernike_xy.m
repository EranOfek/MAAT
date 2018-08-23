function Zxy=zernike_xy(n,m,VecX,VecY,Norm,PutNaN)
% A driver function for the external zerfun.m function.
% Package: telescope.Optics
% Description: A driver function for the external zerfun.m function.
%              Given a list of zernike indices and a grid in the X-Y plan
%              calculate a cube of the zernike functions. The first
%              and second dimensions of the cube corresponds to the
%              X and Y axes, while the third dimension to the zernike
%              index.
% Input  : - Row vector of Zernike n indices. Alternatively, if the
%            second input argument is empty, then this is the Noll index
%            (see noll_index.m).
%            Default is (1:1:10).
%          - Row vector of Zernike m indices. Default is empty.
%          - Vector of X in the X-Y grid. Default is (-1:0.01:1).'.
%            Alternatively, this can be a matrix of X values returned
%            by meshgrid.
%          - Vector of Y in the X-Y grid. Default is (-1:0.01:1).'.
%            Alternatively, this can be a matrix of Y values returned
%            by meshgrid.
%          - Normalization method:
%            'no'   - no normalization.
%            'norm' - Normalized Zernike functions.
%                     See zernfun.m for details.
%            Alternatively, if this is a number (of active pixels) than
%            each Zernike function is normalized such that
%            \sum_{x,y}{z^2(x,y)} = N, where N is the number
%            of active pixels. Default.
%          - A flag {true|false} indicating if to use NaN's outside
%            the circle or zero. Default is true (i.e., put NaN's).
% Output : - A cube of the Zernike functions. The first
%            and second dimensions of the cube corresponds to the
%            X and Y axes, while the third dimension to the zernike
%            index.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Zxy=zernike_xy; surface(sum(Zxy,3)); shading interp
% Reliable: 2
%--------------------------------------------------------------------------

import telescope.Optics.*

Def.n    = (1:1:10);
Def.m    = [];
Def.VecX = (-1:0.01:1).';
Def.VecY = (-1:0.01:1).';
Def.Norm = 'z2';
Def.PutNaN  = true;

if (nargin==0)
    n     = Def.n;
    m     = Def.m;
    VecX  = Def.VecX;
    VecY  = Def.VecY;
    Norm  = Def.Norm;
    PutNaN= Def.PutNaN;
elseif (nargin==1)
    m     = Def.m;
    VecX  = Def.VecX;
    VecY  = Def.VecY;
    Norm  = Def.Norm;
    PutNaN= Def.PutNaN;
elseif (nargin==2)
    VecX  = Def.VecX;
    VecY  = Def.VecY;
    Norm  = Def.Norm;  
    PutNaN= Def.PutNaN;
elseif (nargin==3)
    VecY  = Def.VecY;
    Norm  = Def.Norm;  
    PutNaN= Def.PutNaN;
elseif (nargin==4)
    Norm  = Def.Norm;  
    PutNaN= Def.PutNaN;
elseif (nargin==5)
    PutNaN= Def.PutNaN;
elseif (nargin==6)
    % do nothing
else
    error('Illegal number of input arguments');
end

NormPar    = cell(0,1);
if (ischar(Norm))
    switch lower(Norm)
        case 'norm'
            NormPar{1} = 'norm';
        otherwise
            % do nothing
    end
end

if (isempty(m))
    % convert j to [n,m]
    [n,m] = noll_index(n.');
    n = n.';
    m = m.';
end

Nj   = length(n);

if (numel(VecX)>max(size(VecX)))
    % VecX and VecY are in matrix form already
    MatX    = VecX;
    MatY    = VecY;
    [Nx,Ny] = size(MatX);
    
else
    
    Nx   = length(VecX);
    Ny   = length(VecY);

    [MatX,MatY] = meshgrid(VecX,VecY);
end
[Theta,R]   = cart2pol(MatX,MatY);
%Idx = R<=1;
Z          = nan(size(MatX));
if (PutNaN)
    R(R>1)     = NaN;
else
    R(R>1)     = 0;
end
ZZ         = Util.external.zernfun(n,m,R(:),Theta(:)); %,NormPar);
if (isnumeric(Norm))
    % noramlize such that sum(z^2)=N
    ZZ = bsxfun(@times,ZZ,sqrt(Norm./nansum(ZZ.^2,1)));
end
%nansum(ZZ.^2,1) should equal Norm
Zxy        = reshape(ZZ,Ny,Nx,Nj);


% (Normalization) verify that:
%std(SumY(find(SumY~=0))) = sqrt(sum(A(2:end).^2))


