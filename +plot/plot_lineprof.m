function [Prof,VecR,VecX,VecY]=plot_lineprof(Width,Step,Stat,Interp)
%-----------------------------------------------------------------------------
% plot_lineprof function                                             plotting
% Description: Clicking on two points in current image the script return
%              the intensity as function of position along the line
%              between the two points.
% Input  : - Width of the line (odd integer), default is 1.
%          - Approximate step size along the line [pixels], default is 1.
%            The actual step size is adjusted so nultiply it by an integer
%            will give the length of the line.
%          - Value to calculate:
%            'sum' | 'mean' | 'median' | 'std' | 'min' | 'max',
%            default is 'mean'.
%          - Interpolation method (see interp1), default is 'linear'
% Output : - Vector of profile.
%          - Vector of position along the line [pixels]
%          - Corresponding X position [pixels]
%          - Corresponding Y position [pixels]
% Tested : Matlab 5.3
%     By : Eran O. Ofek                       Feb 2004
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%-----------------------------------------------------------------------------
DefWidth    = 1;
DefStep     = 1;
DefStat     = 'mean';
DefInterp   = 'linear';

if (nargin==0),
   Width    = DefWidth;
   Step     = DefStep;
   Stat     = DefStat;
   Interp = DefInterp;
elseif (nargin==1),
   Step     = DefStep;
   Stat     = DefStat;
   Interp = DefInterp;
elseif (nargin==2),
   Stat     = DefStat;
   Interp = DefInterp;
elseif (nargin==3),
   Interp = DefInterp;
elseif (nargin==4),
   % no default
else
   error('Illegal number of input arguments');
end

if (Width~=floor(Width) | floor(0.5.*(Width+1))~=0.5.*(Width+1)),
   error('Width should be an odd integer');
end

% Read Image to memory
Image  = get(get(gca,'Children'),'CData');
ImX    = get(get(gca,'Children'),'XData');
ImY    = get(get(gca,'Children'),'YData');

[X,Y]  = ginput(2);
DiffX  = diff(X);
DiffY  = diff(Y); 
SlopeL = atan2(DiffY,DiffX);   % slope of the line
SlopeV = atan2(DiffX,DiffY);   % slope of vertical to line

Length = sqrt(DiffX.^2+DiffY.^2);
Nstep  = ceil(Length./Step);
Step   = Length./Nstep;
VecR   = [0:Step:Length].';
VecX   = [X(1):DiffX./Nstep:X(2)].';
VecY   = [Y(1):DiffY./Nstep:Y(2)].';
HalfWidth = 0.5.*(Width - 1);

Data = zeros(Width,Nstep+1);
K    = 0;
for PosWidth=-HalfWidth:1:HalfWidth,
   K       = K + 1;
   CurX    = X + PosWidth.*cos(SlopeV);
   CurY    = Y - PosWidth.*sin(SlopeV);
   
   CurVecX = [CurX(1):DiffX./Nstep:CurX(2)].';
   CurVecY = [CurY(1):DiffY./Nstep:CurY(2)].';
   Data(K,:) = interp2(ImX,ImY,Image,CurVecX,CurVecY,Interp).';

end


switch Stat
 case 'mean'
    Prof = mean(Data,1);
 case 'std'
    Prof = mean(Data,0,1);
 case 'median'
    Prof = median(Data,1);
 case 'sum'
    Prof = sum(Data,1);
 case 'min'
    Prof = min(Data,[],1);
 case 'max'
    Prof = max(Data,[],1);
 otherwise
    error('Unknown statistics option');
end
