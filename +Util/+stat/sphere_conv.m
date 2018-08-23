function ConvMat=sphere_conv(Long,Lat,Mat,Weight,Kernel)
%--------------------------------------------------------------------------
% sphere_conv function                                           AstroStat
% Description: Convolve a 2D matrix which point coordinates are on the
%              celestial sphere, with a kernel.
% Input  : - Vector of longitudes [rad].
%          - Vector of latitudes [rad].
%          - Matrix [Long,Lat].
%          - Matrix of weights (or area per cell).
%          - Kernel [radius(rad), amplitude].
% Output : - Convolved matrix.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ConvMat=sphere_conv(Long,Lat,Mat,Weight,Kernel)
% Reliable: bug (not normalized by area)
%--------------------------------------------------------------------------
InterpMethod = 'linear';

[MatLong, MatLat] = meshgrid(Long,Lat);
N = numel(MatLong);
ConvMat = zeros(size(Mat));
for I=1:1:N,
    D  = sphere_dist(MatLong(I),MatLat(I),MatLong(:),MatLat(:));
    KV = Weight(:).*interp1(Kernel(:,1),Kernel(:,2),real(D),InterpMethod);
    ConvMat(I) = meannd(Mat(:).*KV);
end
    