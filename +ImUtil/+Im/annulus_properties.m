function [BackVal,ErrVal]=annulus_properties(Image,X,Y,AnnRad,BackFun,ErrFun,MatX,MatY)
% Annulus background and noise around coordinates in an image.
% Package: ImUtil.Im
% Description: Annulus background and noise around coordinates in an image.
% Input  : - Image matrix.
%          - X coordinate of annulus center.
%          - Y coordinate of annulus center.
%          - A vector of inner and outer annulus radii.
%          - Annulus background estimator function handle.
%            Default is @median.
%          - Annulus std estimator function handle.
%            Default is @std.
%          - Matrix of image X coordinates.
%            Default is to calculate it using meshgrid.
%          - Matrix of image Y coordinates.
% Output : - Annulus background estimator.
%          - Annulus std/sqrt(N) estimator.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [BackVal,ErrVal]=ImUtil.Im.annulus_properties(S.Im,X,Y,[4 6])
% Reliable: 2
%--------------------------------------------------------------------------

Def.BackFun     = @median;
Def.ErrFun      = @std;
Def.MatX        = [];
Def.MatY        = [];
if (nargin==4)
    BackFun     = Def.BackFun;
    ErrFun      = Def.ErrFun;
    MatX        = Def.MatX;
    MatY        = Def.MatY;
elseif (nargin==5)
    ErrFun      = Def.ErrFun;
    MatX        = Def.MatX;
    MatY        = Def.MatY;
elseif (nargin==6)
    MatX        = Def.MatX;
    MatY        = Def.MatY;
else
    % do nothing
end

[NY,NX] = size(Image);

if (isempty(MatX))
    [MatX,MatY] = meshgrid((1:1:NX),(1:1:NY));
end

R1 = AnnRad(1);
R2 = AnnRad(2);

X1 = max(1, floor(X - R2));
X2 = min(NX, ceil(X + R2));
Y1 = max(1, floor(Y - R2));
Y2 = min(NY, ceil(Y + R2));

% trim image for faster execuation
SubImage = Image(Y1:Y2,X1:X2);
SubMatX  = MatX(Y1:Y2,X1:X2);
SubMatY  = MatY(Y1:Y2,X1:X2);
SubMatR2 = (SubMatX - X).^2 + (SubMatY - Y).^2;
FlagAnn  = SubMatR2>R1.^2 & SubMatR2<=R2.^2;

% calculate annulus background level
%BackVal  = BackFun(SubImage(FlagAnn));
BackVal  = double(BackFun(SubImage(FlagAnn))); % Na'ama, 20180828

if (nargout>1)
    % calculate annulus noise level
    Npix     = sum(FlagAnn(:));
    %ErrVal   = ErrFun(SubImage(FlagAnn))./sqrt(Npix);
    ErrVal   = ErrFun(double(SubImage(FlagAnn)))./sqrt(Npix); % Na'ama, 20180828
end





