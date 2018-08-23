function predict_images(ModelPars,ModelType,StartMethod,ImPos,ThetaX,ThetaY,Dls_Ds,Threshold,Color);
%----------------------------------------------------------------------
% predict_images function                                        glens
% Description:
% Input  : - Model parameters or AlphaX
%          - Model Type       or AlphaY
%          - Start position:
%            'I' - InPos is image position.
%            'S' - ImPos is source position.
%          - Image position [X, Y]
%          - ThetaX matrix to search in
%          - ThetaY matrix to search in
%          - Dls/Ds
%          - Image search threshold (distance between sources in source plane).
%          - Color
% Output : -
% Tested : Matlab 6.5
%     By : Eran O. Ofek       March 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------------
PlotSource = 'yes';
if (nargin==7),
   Color = [];
end
InterpMethod = 'linear';

%--- ModelType vector ---
ColModelType   = 1;   % Model type
%--- ModelPars matrix ---
ColX0          = 1;   % Potential center X [?]
ColY0          = 2;   % Potential center Y [?]
ColE           = 3;   % eccentricity
ColPA          = 4;   % PA [rad]
ColCore        = 5;   % Core radius [?]
ColRs          = 6;   % ScDls_Ds,ale radius [?]
ColGamma       = 7;   % power-law (1-for NFW)
ColNorm        = 8;   % density normalization [?]
ColPert        = 9;   % perturbation normalization [?]
% chebyshec 2
ColT10         = 10;
ColT20         = 11;
ColT01         = 12;
ColT02         = 13;
ColT11         = 14;
% chebyshec 3
ColT30         = 15;
ColT03         = 16;
ColT21         = 17;
ColT12         = 18;

if (min(size(ModelType))>1),
   % assume AlphaX and AlphaY are given instead of ModelPars and ModelType
   AlphaX = ModelPars;
   AlphaY = ModelType;
else

   %--- Number of components in the model ---
   Ncomp = size(ModelPars,1);

   %SizeT  = size(ThetaX);
   %ThetaX = Util.array.mat2vec(ThetaX);
   %ThetaY = Util.array.mat2vec(ThetaY);


   [AlphaX,AlphaY,A_11,A_22,A_12] = calc_alpha(ThetaX,ThetaY,ModelPars,ModelType);


   save AlphaX.mat AlphaX
   save AlphaY.mat AlphaY

end

Nd = length(Dls_Ds);
ColorMap = colormap('jet');
colormap('gray');
SizeCM   = size(ColorMap,1);
ColorInd = floor([1:SizeCM./Nd:SizeCM].');
for Id=1:1:Nd,

   switch StartMethod
    case 'I'
       %--- start with image position ---
       %--- deflection at THE image position ---
       AlphaX1 = interp2(ThetaX,ThetaY,AlphaX,ImPos(1),ImPos(2),InterpMethod);
       AlphaY1 = interp2(ThetaX,ThetaY,AlphaY,ImPos(1),ImPos(2),InterpMethod);

       %--- source position for THE image position ---
       Beta1X = ImPos(1) - Dls_Ds(Id).*AlphaX1;
       Beta1Y = ImPos(2) - Dls_Ds(Id).*AlphaY1;
    case 'S'
       %--- start with source position ---
       Beta1X = ImPos(1);
       Beta1Y = ImPos(2);
         
    otherwise
       error('Unknown StartMethod Option');
   end

   %--- source position for the rest of image plane ---
   BetaX  = ThetaX   - Dls_Ds(Id).*AlphaX;
   BetaY  = ThetaY   - Dls_Ds(Id).*AlphaY;

   DiffMat = sqrt( (Beta1X-BetaX).^2 + (Beta1Y-BetaY).^2 );

   %--- search for images using threshold ---
   It = find(DiffMat<Threshold);

   if (isempty(Color)==1),
      plot(ThetaX(It),ThetaY(It),'.','Color',ColorMap(ColorInd(Id),:) );
      hold on;
   else
      plot(ThetaX(It),ThetaY(It),'.','Color',Color );
      hold on;
   end

   switch PlotSource
    case 'yes'
       %--- Plot source position ---
       if (isempty(Color)==1),
          plot(ThetaX(It),ThetaY(It),'.','Color',ColorMap(ColorInd(Id),:) );
       else
          plot(Beta1X,Beta1Y,'x','Color',Color );
       end
    otherwise
      %--- do nothing ---
   end
end
