function Res=moment2d(MatX,MatY,MatVal,X,Y,varargin)
% Calculate the moments and flux around given coordinates in a 2D image.
% Package: ImUtil.Im
% Description: Given an image, calculate the moments and flux around
%              given set of coordinates.
%              OBSOLETE: Use ImUtil.Im.im_moments.m instead.
% Input  : - Either a matrix of X coordinates corresponding to the
%            value matrix (third input argument), or a vector of X
%            coordinates. If empty then assume 1:size grid.
% 	       - Either a matrix of Y coordinates corresponding to the
%            value matrix (third input argument), or a vector of Y
%            coordinates. If empty then assume 1:size grid.
%          - Matrix of values (image).
%          - Vector of X coordinates around to calculate the moments,
%            and flux.
%          - Vector of Y coordinates around to calculate the moments,
%            and flux.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            Keywords and values can be one of the followings:
%            'Radius'  - Radius in which to calculate the first and
%                        second moments. Default is 10.
%            'AperR'   - Vector of aperture photometry radii.
%                        Default is [4,6,8,10].
%            'Niter'   - Number of centering iterations. Default is 3.
%            'Annulus' - Inner and outer background annulus radii.
%                        Default is [20 30].
%            'BackFun' - Handle to a function for background calculation.
%                        Default is @median.
%            'Gain'    - CCD gain. Default is 1.
%            'RN'      - CCD readout noise [e-]. Default is 6 e-.
%            'ZP'      - Zero point to calculate magnitudes. Default is 22.
% Output : Structure array, with element per X,Y input coordinates,
%          with the following fields:
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=ImUtil.Im.moment2d(MatX,MatY,MatVal,[],[]);
% Reliable: 2
%--------------------------------------------------------------------------

Def.X = [];
Def.Y = [];
Def.Radius = 10;
Def.Niter  = 3;
if (nargin==3)
    X = Def.X;
    Y = Def.Y;
elseif (nargin==4)
    Y = Def.Y;
else
    % do nothing
end

DefV.Radius  = 10;
DefV.AperR   = [4,6,8,10];
DefV.Niter   = 3;
DefV.Annulus = [20 30];
DefV.BackFun = @median;
DefV.Gain    = 1;
DefV.RN      = 6;
DefV.ZP      = 22;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Naper = numel(InPar.AperR);

[SizeY, SizeX] = size(MatVal);
if (isempty(MatX) || isempty(MatY))
    VecX = (1:1:SizeX).';
    VecY = (1:1:SizeY).';
    [MatX,MatY] = meshgrid(VecX,VecY);
elseif (min(size(MatX))==1 || min(size(MatY))==1)
    [MatX,MatY] = meshgrid(MatX,MatY);
else
    % do nothing
end

if (isempty(X))
    X = mean(MatX(1,:));
end
if (isempty(Y))
    Y = mean(MatY(:,1));
end
Res.InitX = X;
Res.InitY = Y;

Ncoo = numel(X);
MeanBack = InPar.BackFun(MatVal(:));
MatValBS = MatVal - MeanBack;

for Icoo=1:1:Ncoo
    
    for Iiter=1:1:InPar.Niter

        MatR    = sqrt( (MatX-X(Icoo)).^2 + (MatY-Y(Icoo)).^2);
        FlagR   = MatR<InPar.Radius;
        SumR    = sum(MatValBS(FlagR));
        X(Icoo) = sum(MatX(FlagR).*MatValBS(FlagR))./SumR;
        Y(Icoo) = sum(MatY(FlagR).*MatValBS(FlagR))./SumR;
    end

    X2 =    sum(MatX(FlagR).^2.*MatValBS(FlagR))./SumR - X(Icoo).^2;
    Y2 =    sum(MatY(FlagR).^2.*MatValBS(FlagR))./SumR - Y(Icoo).^2;
    XY =    sum(MatX(FlagR).*MatY(FlagR).*MatValBS(FlagR))./SumR - X(Icoo).*Y(Icoo);

    Res(Icoo).X  = X(Icoo);
    Res(Icoo).Y  = Y(Icoo);
    Res(Icoo).X2 = X2;
    Res(Icoo).Y2 = Y2;
    Res(Icoo).XY = XY;
    [Max,MaxI] = max(MatVal(:));
    Res(Icoo).PeakX  = MatX(MaxI);
    Res(Icoo).PeakY  = MatY(MaxI);
    Res(Icoo).PeakVal= Max;
    
    FlagB   = MatR>InPar.Annulus(1) & MatR<InPar.Annulus(2);
    Res(Icoo).AreaBack = numel(find(FlagB));
    Res(Icoo).Back     = InPar.BackFun(MatVal(FlagB));

    for Iaper=1:1:Naper
        FlagR   = MatR<InPar.AperR(Iaper);
        Res(Icoo).Sum(Iaper)    = sum(MatVal(FlagR)); 
        Res(Icoo).Area(Iaper)   = numel(find(FlagR));
        Res(Icoo).Flux(Iaper)   = Res(Icoo).Sum(Iaper) - Res(Icoo).Back.*Res(Icoo).Area(Iaper);
        Res(Icoo).SN(Iaper)     = Res(Icoo).Flux(Iaper).*InPar.Gain./sqrt(Res(Icoo).Sum(Iaper).*InPar.Gain + InPar.RN.^2);
        Res(Icoo).Mag(Iaper)    = InPar.ZP - 2.5.*log10(Res(Icoo).Flux(Iaper));
        Res(Icoo).MagErr(Iaper) = 1.086./Res(Icoo).SN(Iaper);
    end
end
Res.AperR   = InPar.AperR;
Res.Annulus = InPar.Annulus; 
Res.ZP      = InPar.ZP;



