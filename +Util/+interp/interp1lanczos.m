function YI=interp1lanczos(V,XI,A)
% 1-D Lanczos interpolation
% Package: Util.interp
% Description: 1D Lanczos interpolation.
% Input  : - Equally spaced vector to interpolate.
%          - Positions to which to interpolate.
%          - Lanczos order (2|3). Default is 2.
% Output : - Interpolated vector.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: YI=Util.interp.interp1lanczos([1;1;1.05;1.1;1.05;1;1],1.5)
% Reliable: 
%--------------------------------------------------------------------------

if (nargin==2)
    A = 2;
end


% Make V a row vector
V  = V(:).';
NV = numel(V);

X  = (1:1:NV);

if (A==2)
    MatV = zeros(5,NV+4);
    MatV(1,:) = [V, 0, 0, 0, 0];
    MatV(2,:) = [0, V, 0, 0, 0];
    MatV(3,:) = [0, 0, V, 0, 0];
    MatV(4,:) = [0, 0, 0, V, 0];
    MatV(5,:) = [0, 0, 0, 0, V];

    MatX = zeros(5,NV+4);
    MatX(1,:) = [X, 0, 0, 0, 0];
    MatX(2,:) = [0, X, 0, 0, 0];
    MatX(3,:) = [0, 0, X, 0, 0];
    MatX(4,:) = [0, 0, 0, X, 0];
    MatX(5,:) = [0, 0, 0, 0, X];
elseif (A==3)
    MatV = zeros(7,NV+6);
    MatV(1,:) = [V, 0, 0, 0, 0, 0, 0];
    MatV(2,:) = [0, V, 0, 0, 0, 0, 0];
    MatV(3,:) = [0, 0, V, 0, 0, 0, 0];
    MatV(4,:) = [0, 0, 0, V, 0, 0, 0];
    MatV(5,:) = [0, 0, 0, 0, V, 0, 0];
    MatV(6,:) = [0, 0, 0, 0, 0, V, 0];
    MatV(7,:) = [0, 0, 0, 0, 0, 0, V];

    MatX = zeros(7,NV+6);
    MatX(1,:) = [X, 0, 0, 0, 0, 0, 0];
    MatX(2,:) = [0, X, 0, 0, 0, 0, 0];
    MatX(3,:) = [0, 0, X, 0, 0, 0, 0];
    MatX(4,:) = [0, 0, 0, X, 0, 0, 0];
    MatX(5,:) = [0, 0, 0, 0, X, 0, 0];
    MatX(6,:) = [0, 0, 0, 0, 0, X, 0];
    MatX(7,:) = [0, 0, 0, 0, 0, 0, X];

    
else
    error('Lanczos order must be 2 | 3');
end

I   = floor(XI);
Eta = bsxfun(@minus,MatX(:,I+A),XI);
L   = lanczos(Eta,A);
YI = sum(MatV(:,I+A).*L,1);
