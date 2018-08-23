function [MaxX,MaxY,ErrX,Flag]=find_peak(X,Y,X0,Method,SemiWidth,StirSel,StirWin)
%------------------------------------------------------------------------------
% find_peak function                                                   General
% Description: Given a tabulated function [X,Y], find the maximum
%              near a given X0.
% Input  : - X of tabulated function.
%          - Y of tabulated function.
%          - X0 around at which to search for the maximum.
%          - Method for maximum position estimation
%            {'max'|'wmean'|'stirling4'},
%            default is 'stirling4'.
%            The 'stirling4' calls find_local_extramum.m function.
%          - Semi width around X0 in which to search for the maximum,
%            default is 3.
%          - Optional selection method for 'stirling4' - if more than
%            one maximum is found in region, then select the
%            {'nearest'|'highest'}, default is 'heighest'.
%          - Optional half window size for find extramum using the
%            'stirling4' option, default is 10.
% Output : - X position of of the maximum.
%          - Y value at maximum.
%          - Error in X position (only for 'wmean' option).
%          - Maximum type flag:
%             0  : Maximum is found inside search region
%             1  : Maximum is found on boundry of left side region,
%                  however it is higher than the value Y to its left.
%             2  : Maximum is found on boundry of left side region,
%                  however it is lower than the value Y to its left.
%             3: : Search region is out of vector (left) range.
%            -1  : Maximum is found on boundry of right side region,
%                  however it is higher than the value Y to its right.
%            -2  : Maximum is found on boundry of right side region,
%                  however it is lower than the value Y to its right.
%            -3: : Search region is out of vector (right) range.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                     Mar 2007
%    URL : http://wweizmann.ac.il/home/eofek/matlab/
% Example: X=(1:1:10).'; Y=ones(size(X)); Y(4:6)=[2 3 2.1];
%          [MaxX,MaxY,ErrX,Flag]=find_peak(X,Y,4);
% Reliable: 2
%------------------------------------------------------------------------------
DefMethod      = 'stirling4';
DefSemiWidth   = 3;
DefStirSel     = 'highest';
DefStirWin     = 10;

if (nargin==3),
   Method    = DefMethod;
   SemiWidth = DefSemiWidth;
   StirSel   = DefStirSel;
   StirWin   = DefStirWin;
elseif (nargin==4),
   SemiWidth = DefSemiWidth;
   StirSel   = DefStirSel;
   StirWin   = DefStirWin;
elseif (nargin==5),
   StirSel   = DefStirSel;
   StirWin   = DefStirWin;
elseif (nargin==6),
   StirWin   = DefStirWin;
elseif (nargin==7),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (min(X)>(X0-SemiWidth)),
   % Out of range (left)
   MaxX = NaN;
   MaxY = NaN;
   ErrX = NaN;
   Flag = 3;
elseif (max(X)<(X0+SemiWidth)),
   % Out of range (right)
   MaxX = NaN;
   MaxY = NaN;
   ErrX = NaN;
   Flag = -3;
else
   Ix = find(abs(X-X0)<=SemiWidth);
   switch lower(Method)
    case 'max'
       [MaxY,Ind] = max(Y(Ix));
       MaxX       = X(Ix(Ind));
       ErrX       = NaN;
    case 'wmean'
       [MaxX,ErrX] = wmean([X(Ix),1./Y(Ix)]);
       [MaxY,Ind] = max(Y(Ix));
%    case 'spline'
%       PP = spline(X(Ix),Y(Ix));
%       ppval(PP,

    case 'stirling4'

       Ixs = find(abs(X-X0)<=StirWin);

       Extram = find_local_extramum(X(Ixs),Y(Ixs),4);

       % select maxima in range
       I=find( (Extram(:,1)-X0)<=SemiWidth & Extram(:,3)<0 );

       switch StirSel
        case 'nearest'
           [M,Ind] = min(abs(X0-Extram(I,1)));
        case 'highest'
           [M,Ind] = max(Extram(I,2));

        otherwise
	   error('Unknown StirSel option');
       end
       MaxX = Extram(I(Ind),1);
       MaxY = Extram(I(Ind),2);
       ErrX = NaN;
       Flag = 0;

       
    otherwise
       error('Unknown Method option');
   end
   
   if (MaxX==min(X(Ix))),
      % MaxX is on (left) edge of serached region
   
      if (MaxY>X(Ix-1)),
         % however a real maximum
         Flag = 1; 
      else
         % unreal maximum
         Flag = 2; 
      end
   elseif (MaxX==max(X(Ix))),
      % MaxX is on (right) edge of serached region
      if (MaxY>X(Ix+1)),
         % however a real maximum
         Flag = -1; 
      else
         % unreal maximum
         Flag = -2; 
      end
   
   else
      % MaxX is not on edge of serached region
      Flag = 0;
   end
end
