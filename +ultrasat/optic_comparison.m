function optic_comparison


%%
RAD = 180./pi;

DetSize = 90;
PixSize = 9e-3;
AssemblyErr = 1.3;
Jitter      = 10.8./3 .*2.35;  % [arcsec] for FWHM
PixDigi     = 0.0;

InnerArea   = 180;  % [deg^2] statistics for inner area
InnerRad    = sqrt(InnerArea./pi);

I = 0;
I = I + 1;
Tel(I).FWHM = Util.IO.load2('home/yossi/matlab/MAAT/+ultrasat/opt8_f327_1asphere_EE50.txt');
Tel(I).FL   = 327;
Tel(I).ScaleAS = 1;


I = I + 1;
Tel(I).FWHM = Util.IO.load2('home/yossi/matlab/MAAT/+ultrasat/opt8_f357_1asphere_EE50.txt');
Tel(I).FL   = 357;
Tel(I).ScaleAS = 1;

I = I + 1;
Tel(I).FWHM = Util.IO.load2('home/yossi/matlab/MAAT/+ultrasat/opt11_f360_2asphere_EE50.txt');
Tel(I).FL   = 360;
Tel(I).ScaleAS = 1.138;

N = numel(Tel);
for I=1:1:N
    Tel(I).FWHM  = Tel(I).FWHM.*Tel(I).ScaleAS;
    
    Tel(I).FOV = DetSize./Tel(I).FL.*RAD;   % [deg] 
    Tel(I).PixScale = PixSize./Tel(I).FL.*RAD.*3600;  % [arcsec]

    Tel(I).EffFWHM = sqrt( (Tel(I).FWHM.*AssemblyErr).^2 + (PixDigi.*Tel(I).PixScale).^2 + Jitter.^2);
   
    Tel(I).WMeanEffFWHM = mean(Tel(I).EffFWHM(:));
    Tel(I).WMeanFWHM = mean(Tel(I).FWHM(:));
    
    Tel(I).MinEffFWHM = min(Tel(I).EffFWHM(:));
    Tel(I).MaxEffFWHM = max(Tel(I).EffFWHM(:));
    
    Tel(I).MatDim = size(Tel(I).FWHM,1);
    
    Step = (2.*Tel(I).FOV)./(Tel(I).MatDim-1);  % map step size [deg]
    
    [MatX,MatY] = meshgrid((-Tel(I).FOV:Step:Tel(I).FOV) , (-Tel(I).FOV:Step:Tel(I).FOV) );
    MatR = sqrt(MatX.^2 + MatY.^2);
    Flag = MatR<InnerRad;
    Tel(I).WMeanFWHMin = mean(Tel(I).FWHM(Flag));
    Tel(I).WMeanEffFWHMin = mean(Tel(I).EffFWHM(Flag));
    
    
    

    
    Tel(I).Grasp = sum(Tel(I).EffFWHM(:).^(-3./2) .* Step.^2   );
end

