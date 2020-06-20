function [Mom,Mom2,Aper]=im_moments(Image,X,Y,Radius,Sigma,MaxIter,RadAnn)
% Calculate 1st and 2nd moments in a specific locations in an image
% Package: ImUtil.Im
% Description: Calculate 1st and 2nd moments in a specific locations in
%              an image. The moments are calculated iteratively,
%              with Gaussian weights.
% Input  : - Image (matrix).
%          - Vector of intial X coordinates in which to calculate the
%            moments.
%          - Vector of intial Y coordinates in which to calculate the
%            moments.
%          - Radius in which to calculate the moments and properties.
%          - Sigma of Gaussian by which to weight the pixels.
%            Default is 1.5.
%          - Maximum number of 1st moment (centeroiding) iterations.
%            If 1 then use initial coordinates.
%            Default is 3.
%          - Annulus radius. Default is [7 9].
% Output : - Structure array of 1st moments for each source.
%          - Structure array of 2nd moments for each source.
%          - Structure array of aperure photometry for each source.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Image = rand(1024,1024); [M,M2,Aper]=ImUtil.Im.im_moments(Image, rand(100,1), rand(100,1),6)
% Reliable: 2
%--------------------------------------------------------------------------
PosConv = 1e-4;

Def.Sigma   = 1.5;
Def.MaxIter = 3;
Def.RadAnn  = [7 9];

if (nargin==4)
    Sigma   = Def.Sigma;
    MaxIter = Def.MaxIter;
    RadAnn  = Def.RadAnn;
elseif (nargin==5)
    MaxIter = Def.MaxIter;
    RadAnn  = Def.RadAnn;
elseif (nargin==6)
    RadAnn  = Def.RadAnn;
elseif (nargin==7)
    % do nothing
else
    error('Illegal number of input arguments: im_moments(Image,X,Y,Radius,Sigma,MaxIter)');
end


Radius2 = Radius.^2;
%Area    = ceil(pi.*(Radius+0.5).^2);

TwoSig2 = 2.*Sigma.^2;

Size = size(Image);
%[MatX,MatY] = meshgrid((1:1:Size(2)),(1:1:Size(1)));

Nsrc = numel(X);

Mom.X     = X;
Mom.Y     = Y;
RoundX    = round(X);
RoundY    = round(Y);
Mom2.X2   = zeros(size(X));
Mom2.Y2   = zeros(size(X));
Mom2.XY   = zeros(size(X));
Aper.Phot = zeros(size(X));
Aper.Npix = zeros(size(X));
Aper.NpixBack = zeros(size(X));
Aper.Back = zeros(size(X));



%XV = zeros(Area,Nsrc).*NaN;
%YV = zeros(Area,Nsrc).*NaN;
%IV = zeros(Area,Nsrc).*NaN;
%W  = zeros(Area,Nsrc).*NaN;

StampSize = Radius+5;
VecX      = (-StampSize:1:StampSize);
VecY      = VecX;
[MatX, MatY] = meshgrid(VecX,VecY);
% for each source
for Isrc=1:1:Nsrc
    StX = max(RoundX(Isrc)+VecX,1);
    StX = min(StX,Size(2));
    StY = max(RoundY(Isrc)+VecY,1);
    StY = min(StY,Size(1));
    
    Stamp = Image(StY,StX);
    
    %for Iiter=1:1:Niter,
    StepX = Inf;
    StepY = Inf;
    Iiter = 0;
    
    Dist2 = (RoundX(Isrc)+MatX(:) - Mom.X(Isrc)).^2 + (RoundY(Isrc)+MatY(:) - Mom.Y(Isrc)).^2;

    Flag = Dist2<Radius2;

    StampFlag = Stamp(Flag);
    W         = exp(-Dist2(Flag)./TwoSig2);

    StampFlagW= StampFlag.*W;
    InvSumStampFlagW = 1./sum(StampFlagW);
    
    % converge on first moment
    while ((abs(StepX)>PosConv || abs(StepY)>PosConv) && Iiter<MaxIter)
        Iiter = Iiter + 1;
        Dist2 = (RoundX(Isrc)+MatX(:) - Mom.X(Isrc)).^2 + (RoundY(Isrc)+MatY(:) - Mom.Y(Isrc)).^2;

        Flag = Dist2<Radius2;

        StampFlag = Stamp(Flag);
        W         = exp(-Dist2(Flag)./TwoSig2);

        StampFlagW= StampFlag.*W;
        InvSumStampFlagW = 1./sum(StampFlagW);
        StepX = sum((RoundX(Isrc)+MatX(Flag)-Mom.X(Isrc)).*StampFlagW).*InvSumStampFlagW;
        StepY = sum((RoundY(Isrc)+MatY(Flag)-Mom.Y(Isrc)).*StampFlagW).*InvSumStampFlagW;
        Mom.X(Isrc) = Mom.X(Isrc) + 2.*StepX;
        Mom.Y(Isrc) = Mom.Y(Isrc) + 2.*StepY;
        % Iiter
        % Mom
    end
    Mom.Iiter(Isrc) = Iiter;
    
    % calculate additional properties
    if (nargout>1)
        % second moments
        Mom2.X2(Isrc) = sum(StampFlagW.*(RoundX(Isrc)+MatX(Flag) - Mom.X(Isrc)).^2).*InvSumStampFlagW;
        Mom2.Y2(Isrc) = sum(StampFlagW.*(RoundY(Isrc)+MatY(Flag) - Mom.Y(Isrc)).^2).*InvSumStampFlagW;
        Mom2.XY(Isrc) = sum(StampFlagW.*(RoundX(Isrc)+MatX(Flag) - Mom.X(Isrc)).*(RoundY(Isrc)+MatY(Flag) - Mom.Y(Isrc))  ).*InvSumStampFlagW;
        
        
        if (nargout>2)
            % aperture photometry
            Aper.Npix(Isrc) = sum(Flag);
            %Aper.Phot(Isrc) = sum(StampFlag);
            
            FlagAnn = Dist2>RadAnn(1).^2 & Dist2<RadAnn(2).^2;          
            Aper.NpixBack(Isrc) = sum(FlagAnn);
            Aper.Back(Isrc)     = sum(Stamp(FlagAnn));
            Aper.Phot(Isrc)     = sum(StampFlag) - Aper.Back(Isrc).*Aper.Npix(Isrc)./Aper.NpixBack(Isrc);
            
        end
    end
end

