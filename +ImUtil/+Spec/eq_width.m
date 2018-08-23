function [LinesProp,ConPos]=eq_width(Spec,ConPos)
%--------------------------------------------------------------------------
% eq_width function                                                 ImSpec
% Description: Given a spectrum calculate the equivalent width, flux, and 
%              center of selected spectral lines.
%              The script has interactive and noninteractive modes.
%              * Interactive mode: The user select (left click) pairs of
%              wavelength positions marking the continuum around each line
%              (right click to abort).
%              * Noninteractive mode: If a second input argument is given 
%              then the lines properties are calculated noninteractively.
%              This function is obsolete, use: fit_specline.m instead.
% Input  : - Spectrum [wavelength, Intensity].
%          - Continuum position [X1 X2] around each line.
% Output : - Structure array of properties for each line.
%            The following fields are avialable:
%            .WaveMaxInt  - Wavelength at extramum intensity
%            .LineCenter  - Line center (first moment),
%            .LineStd     - Line second moment,
%            .FWHM        - Line FWHM
%            .Peak        - Extramum value,
%            .Flux        - Flux of line,
%            .EW          - Line eqwivalent width].
% See also: fit_specline.m
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [LinesProp,ConPos]=eq_width(get_spectra('QSO_LBQS'));
% Reliable: 2
%--------------------------------------------------------------------------
ColW       = 1;
ColI       = 2;
SpecMarker = 'b-';
ContMarker = 'r--';
Ncol       = 7;

I = find(isnan(Spec(:,2))==0);
Spec = Spec(I,:);

if (nargin<2),
   %--- interactive ---
   LinesProp = zeros(0,Ncol);
   ConPos    = zeros(0,2);

   if (nargin==1),
      Hfig = figure;
      stairs(Spec(:,ColW),Spec(:,ColI)); %,SpecMarker);
   else
      % get data from figure
      H    = get(gca,'Children');
      Spec = [get(H,'XData').', get(H,'YData').'];
   end
   hold on;

   Stop = 0;
   while (Stop==0), 
      zoom on;
      R=input('Zoom in/out, press any key to continue');
      zoom off;
      disp(sprintf('Select line by right click followed by left click'))
      disp(sprintf('Use double left click to abort'))
      figure(Hfig);
      [X,Y] = getline;
      %[X,Y] = ginput(2);
      %Stop=1; 
      if (length(X)==1),
         Stop = 1;
      else
         figure(Hfig);
         [ContX,ContY,Prop] = measure_line_prop(Spec,X(1),X(2));
         plot(ContX,ContY,ContMarker);
         LinesProp = [LinesProp; Prop];
         ConPos    = [ConPos; [X(1), X(2)]];
      end
   end
else
   %--- noninteractive ---
   Nl        = size(ConPos,1);
   LinesProp = zeros(Nl,Ncol);
   for Il=1:1:Nl,
      [ContX,ContY,Prop] = measure_line_prop(Spec,ConPos(Il,ColW),ConPos(Il,ColI));

      LinesProp(Il,:) = Prop;
   end
end


%--------------------------
%--- Measure Properties ---
%--------------------------
function [ContX,ContY,Prop]=measure_line_prop(Spec,X1,X2);
%--------------------------
ColW = 1;
ColI = 2;
InterpMethod = 'linear';

Xmin = min([X1;X2]);
Xmax = max([X1;X2]);

%--- spectrum intensity at Xmin and Xmax wavelength ---
SpecImin   = interp1(Spec(:,ColW),Spec(:,ColI),Xmin,InterpMethod);
SpecImax   = interp1(Spec(:,ColW),Spec(:,ColI),Xmax,InterpMethod);

Cont1      = [Xmin, SpecImin];  
Cont2      = [Xmax, SpecImax];  
ContX      = [Xmin; Xmax];
ContY      = [SpecImin; SpecImax];

%--- indices of line region ---
[Min,Ind1] = min(abs(Spec(:,ColW) - Xmin)); 
[Min,Ind2] = min(abs(Spec(:,ColW) - Xmax)); 
LineRegion = [Ind1:1:Ind2].';

[Max,MaxI] = max(abs(Spec(LineRegion,ColI)));

WaveExt    = Spec(LineRegion(MaxI),ColW);
Moment1    = sum(Spec(LineRegion,ColI).*Spec(LineRegion,ColW))./sum(Spec(LineRegion,ColI));
Moment2    = sum(Spec(LineRegion,ColI).*(Spec(LineRegion,ColW)-Moment1).^2)./sum(Spec(LineRegion,ColI));
ExtramumVal= Spec(LineRegion(MaxI),ColI);
TotFlux    = trapz(Spec(LineRegion,ColW),Spec(LineRegion,ColI));
ContFlux   = trapz(ContX,ContY);
LineFlux   = TotFlux - ContFlux;

% linearly interpolate contimuum within Region
ContInterp = ContY(1) + (Spec(LineRegion,ColW) - ContX(1)).*diff(ContY)./diff(ContX); 

EqWidth    = trapz(Spec(LineRegion,ColW), 1 - Spec(LineRegion,ColI)./ContInterp);


HalfMaxLevel = mean(ContY) + 0.5.*(ExtramumVal - mean(ContY));
AllIndHalfMaxCross = find(abs(diff(sign(Spec(LineRegion,ColI) - HalfMaxLevel)))==2)+1;
FWHM  = Spec(LineRegion(max(AllIndHalfMaxCross)),ColW) - Spec(LineRegion(min(AllIndHalfMaxCross)),ColW);




Prop       = [WaveExt, Moment1, Moment2, FWHM, ExtramumVal, LineFlux, EqWidth];

