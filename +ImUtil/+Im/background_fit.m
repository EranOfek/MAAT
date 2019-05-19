function [BackIm,NoiseIm,Back]=background_fit(Image,varargin)
% SHORT DESCRIPTION HERE
% Package: ImUtil
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [BackIm,NoiseIm,Back]=ImUtil.Im.background_fit(Image);
% Reliable: 
%--------------------------------------------------------------------------


DefV.BinSize              = [32 32]; %[64 64];   % X Y
DefV.StepFraction         = 0.5;
DefV.InterpMethod         = 'spline';
DefV.Sampling             = 2;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (numel(InPar.BinSize)==1)
    InPar.BinSize = [InPar.BinSize, InPar.BinSize];
end

[SizeY, SizeX] = size(Image);

%Median = median(Image(:));
%RStd   = Util.stat.rstd(Image(:));

HalfBinSize = InPar.BinSize.*0.5;

Step     = InPar.BinSize.*InPar.StepFraction;
BinCenX  = unique([(1:Step(1):SizeX).'; SizeX]);
BinCenY  = unique([(1:Step(2):SizeY).'; SizeY]);
Nx       = numel(BinCenX);
Ny       = numel(BinCenY);

% remove obvious stars from the image
%PSF = Kernel2.gauss(2,2);
%[FiltIm,~,CorrF] = 


if (InPar.Sampling == 1)
    [Mode,Noise] = Util.stat.mode_fit(Image);
else
    [Mode,Noise] = Util.stat.mode_fit(Image(1:InPar.Sampling:end,1:InPar.Sampling:end));
end
Variance = Noise.^2;

OneSigLow = normcdf(-1,0,1); 

Back.LowQuantile = zeros(Ny,Nx);
Back.MidQuantile = zeros(Ny,Nx);
Back.Mode        = zeros(Ny,Nx);
Back.Noise       = zeros(Ny,Nx);
Back.Gain        = zeros(Ny,Nx);

for Ix=1:1:Nx
    for Iy=1:1:Ny
        J1 = BinCenX(Ix)-HalfBinSize(1);
        J2 = BinCenX(Ix)+HalfBinSize(1);
        I1 = BinCenY(Iy)-HalfBinSize(2);
        I2 = BinCenY(Iy)+HalfBinSize(2);
        
        I1 = max(I1,1);
        I2 = min(I2,SizeY);
        J1 = max(J1,1);
        J2 = min(J2,SizeX);
        
        SubImage = Image(I1:I2,J1:J2);

        % faster than prctile
        Ns = numel(SubImage);
        SSIm = sort(SubImage(:));
        Back.LowQuantile(Iy,Ix) = SSIm(floor(Ns.*OneSigLow));
        %Back.LowQuantile(Iy,Ix) = prctile(SubImage(:),OneSigLow.*100);
        % faster than median
        Back.MidQuantile(Iy,Ix) = SSIm(floor(Ns.*0.5));
        %Back.MidQuantile(Iy,Ix) = median(SubImage(:));
        Back.Noise(Iy,Ix)       = Back.MidQuantile(Iy,Ix) - Back.LowQuantile(Iy,Ix);
        Back.Gain(Iy,Ix)        = Back.LowQuantile(Iy,Ix)./Variance;
        
    end
end
MedianGain = median(Back.Gain(:));
Back.Mode  = Back.LowQuantile + Noise.*Back.Gain./MedianGain;

% extend to full pixels
[MatX,MatY] = meshgrid((1:1:SizeX),(1:1:SizeY));
BackIm  = interp2(BinCenX,BinCenY,Back.Mode,MatX,MatY,InPar.InterpMethod);
NoiseIm = interp2(BinCenX,BinCenY,Back.Noise,MatX,MatY,InPar.InterpMethod);


