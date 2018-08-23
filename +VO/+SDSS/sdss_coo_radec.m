function [RA,Dec]=sdss_coo_radec(FrameRow,FrameCol,Matrix,Filter,Color,Margin,Col)
% Convert the SDSS great circles coordinate system to J2000 RA and Dec.
% Package: VO.SDSS
% Description: Convert the SDSS great circles coordinate system to J2000
%              right ascension and declination.
% Input  : - Matrix of X Coordinates (frame row).
%            Each line for each field, where lines can contain
%            a single coordinate or multiple coordinates.
%            If empty, then default is X position of the four corners
%            of the field [65:1425,1:2048].
%          - Matrix of Y Coordinates (frame col).
%            Each line for each field, where lines can contain
%            a single coordinate or multiple coordinates.
%            If empty, then default is Y position of the four corners
%          - Matrix of SDSS coordinates and coef.
%            (e.g., SDSS_DR5_Fields_* files)
%          - Filter of image in which the row and col are specified
%            {'g'|'r'|'i'}, default is 'g'.
%          - Color of object for which to calculate poition,
%            default is 0 (i.e., do not apply color correction).
%            color should be g-r for g-band, and r-i for r and i-bands.
%            UNSPORTED!
%          - If frame row and frame col are empty, then subtract/add
%            margin in pixels to the position of corners, default is 0.
%          - Columns info structure:
%            .A - columns in matrix containing alpha (Node),
%                 default is 69.
%            .I - column in matrix containing inclination.
%                 default is 68.
%            .GH- columns containing
%                 [DROW0_G, DROW1_G, DROW2_G, DROW3_G, DCOL0_G, DCOL1_G, DCOL2_G, DCOL3_G],
%                 default is [14:1:21].
%            .AF- columns containing [A_G, B_G, C_G, D_G, E_G, F_G],
%                 default is [24:1:29].
% Output : - RA of corners [rad].
%            Dec of corners [rad].
% Tested : Matlab 7.0                     
%     By : Eran O. Ofek / Dovi Poznanski   Jan 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: TargetFields_SDSS_DR4.head
% Reference: http://www.sdss.org/dr5/products/general/astrometry.html
% Example: % Calculating the RA/Dec of corners of frames:
%            [RA,Dec]=VO.SDSS.sdss_coo_radec([],[],SDSS_DR5_Fields_Best(1:10,:));
% Reliable: 2
%------------------------------------------------------------------------------


RAD   = 180./pi;

%Yul   = 2048;
%Xul   = 1489;   % entire field including overlap (with next field) region.
%Xul   = 1361;   % unique field

% unique field definition:
Yll   = 1;
Yul   = 2048;
Xll   = 65;
Xul   = 1425;


DefFilter = 'g';
DefColor  = 0;
DefMargin = 0;

if (nargin==3)
   Filter     = DefFilter;
   Color      = DefColor;
   Margin     = DefMargin;
elseif (nargin==4)
   Color      = DefColor;
   Margin     = DefMargin;
elseif (nargin==5)
   Margin     = DefMargin;
elseif (nargin==6)
   % do nothing
elseif (nargin==7)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (nargin~=7)
   switch Filter
    case 'g'
       Col.GH = [14:1:21];
       Col.AF = [24:1:29];
       Col.CCrow = 22;
       Col.CCcol = 23;
       Col.CSrow = NaN;    
       Col.CScol = NaN;
       Col.I  = 68;
       Col.A  = 69;
       Color0 = 1.5;   % g-r
    case 'r'
       Col.GH = [32:1:35];
       Col.AF = [36:1:39];
       Col.CCrow = 40;
       Col.CCcol = 41;
       Col.CSrow = NaN;    
       Col.CScol = NaN;
       Col.I  = 68;
       Col.A  = 69;
       Color0 = Inf;   % r-i
    case 'i'
       Col.GH = [50:1:53];
       Col.AF = [54:1:57];
       Col.CCrow = 58;
       Col.CCcol = 59;
       Col.CSrow = NaN;    
       Col.CScol = NaN;
       Col.I  = 68;
       Col.A  = 69;
       Color0 = Inf;   % r-i
    otherwise
       error('Unknown Filter option');
   end
end

if (isempty(FrameRow)==1)
   % set X to corners
   X = [Xll+Margin, Xul-Margin, Xul-Margin, Xll+Margin];
end
if (isempty(FrameCol)==1)
   % set Y to corners
   Y = [Yll+Margin, Yll+Margin, Yul-Margin, Yul-Margin];
end

N   = size(Matrix,1);
Nx  = size(X,2);
Ny  = size(Y,2);

if (size(X,1)==1)
   X = ones(N,1)*X;
end
if (size(Y,1)==1)
   Y = ones(N,1)*Y;
end
if (size(Color,1)==1)
   Color = ones(N,1)*Color;
end

RA  = zeros(N,Nx);
Dec = zeros(N,Ny);
for I=1:1:N
   In = Matrix(I,:);
   % x' = x + g0 + g1 y + g2 y^2 + g3 y^3
   % y' = y + h0 + h1 y + h2 y^2 + h3 y^3
   X_tag = zeros(1,Nx);
   Y_tag = zeros(1,Ny);
   X_tag = X(I,:) + In(Col.GH(1)) + In(Col.GH(2)).*Y(I,:) + ...
                                    In(Col.GH(3)).*Y(I,:).^2 + ...
                                    In(Col.GH(4)).*Y(I,:).^3;
   Y_tag = Y(I,:) + In(Col.GH(5)) + In(Col.GH(6)).*Y(I,:) + ...
                                    In(Col.GH(7)).*Y(I,:).^2 + ...
                                    In(Col.GH(8)).*Y(I,:).^3;

   if (Color(I)==0)
      % donot apply color correction
   else
      if (Color(I)<Color0)
         X_tag = X_tag + In(Col.CSrow).*Color;
      else
         X_tag = X_tag + In(Col.CCrow);
      end
      if (Color(I)>Color0)
         Y_tag = Y_tag + In(Col.CScol).*Color;
      else
         Y_tag = Y_tag + In(Col.CCcol);
      end
   end

   % miu = a + b x' + c y'
   % niu = d + e x' + f y'
   Miu = In(Col.AF(2)).* X_tag + In(Col.AF(3)).* Y_tag + In(Col.AF(1));
   Niu = In(Col.AF(5)).* X_tag + In(Col.AF(6)).* Y_tag + In(Col.AF(4));

   Miu=Miu./RAD;
   Niu=Niu./RAD;
   %                      sin(miu-alpha0) * cos(niu)  * cos(i)   -   sin(niu) * sin(i)
   % tan (RA-alpha0) =     -----------------------------------------------------------
   %                                      cos(miu-alpha0) *cos(niu)
   %
   % sin (DEC) =  sin(miu-alpha0) * cos(niu)*sin(i)  +  sin(niu) * cos(i)

   RA_rhs_nom   = sin(Miu-In(Col.A)./RAD) .* cos(Niu) .* cos(In(Col.I)./RAD)   - sin(Niu) .* sin(In(Col.I)./RAD);
   RA_rhs_denom = cos(Miu-In(Col.A)./RAD) .* cos(Niu);

   RA(I,:)   = [In(Col.A)./RAD + atan2(RA_rhs_nom,RA_rhs_denom)];
   Dec(I,:)  = [asin( sin(Miu-In(Col.A)./RAD) .* cos(Niu) .* sin(In(Col.I)./RAD) + sin(Niu) .* cos(In(Col.I)./RAD) )];
end


