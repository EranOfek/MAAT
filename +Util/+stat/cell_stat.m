function [X,Y,Mn,Mnin,Mmean,Mmedian,Mwmean,Mwerr,Mstd,Ment]=cell_stat(Data,CellSize,SpaceLimit,MinPoints)
%----------------------------------------------------------------------------
% cell_stat function                                               AstroStat
% Description: Given a list of x,y coordinates (with optional property 
%              columns), count the number of points in each cell in the
%              x-y plane and calculate the statistics (e.g., mean, medain
%              [ignore NaNs]) of the optional property in each cell.
% Input  : - matrix of events in 2-d space [X, Y, Value, Err]
%            Value and Err are optional.
%            If Err is not given it is taken to be 1.
%          - Cell size [StepX, StepY].
%          - Space boundry [Xmin, Xmax, Ymin, Ymax].
%          - Minimum number of points for which to calculate statistics,
%            default is 0.
% Output : - Vector of cells X-coordinate position.
%          - Vector of cells Y-coordinate position.
%          - Matrix of number of all points in each cell.
%            If value is not given in input then return NaN.
%          - Matrix of number of non-NaNs points in each cell.
%            If value is not given in input then return NaN.
%          - Matrix of mean value of point in each cell.
%            If value is not given in input then return NaN.
%          - Matrix of median value of point in each cell.
%            If value is not given in input then return NaN.
%          - Matrix of values weighted mean in each cell.
%            If value is not given in input then return NaN.
%          - Matrix of values weighted error in each cell.
%            If value is not given in input then return NaN.
%          - Matrix of Values StD in each cell.
%            If value is not given in input then return NaN.
%          - Matrix of Values information entropy
%            [-sum(V*ln(V))] in each cell.
%            If value is not given in input then return NaN.
% Tested : Matlab 5.3
%     By : Eran O. Ofek         July 2003
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------------------
if (nargin==3),
   MinPoints = 0;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end
ColX     = 1;
ColY     = 2;
ColV     = 3;
ColE     = 4;


SizeData = size(Data);
N        = SizeData(1);
Col      = SizeData(2);

if (Col==2),
   Data = [Data, NaN.*zeros(N,2)];
elseif (Col==3),
   Data = [Data, NaN.*zeros(N,1)];
else
   % do nothing
end


%--- cells coordinate center ---
StepX   = CellSize(1);
StepY   = CellSize(2);
X       = [SpaceLimit(1)+0.5.*StepX:StepX:SpaceLimit(2)].';
Y       = [SpaceLimit(3)+0.5.*StepY:StepY:SpaceLimit(4)].';

Nx      = length(X);
Ny      = length(Y);
Mn      = zeros(Ny,Nx).*NaN;
Mnin    = zeros(Ny,Nx).*NaN;
Mmean   = zeros(Ny,Nx).*NaN;
Mmedian = zeros(Ny,Nx).*NaN;
Mwmean  = zeros(Ny,Nx).*NaN;
Mwerr   = zeros(Ny,Nx).*NaN;
Mstd    = zeros(Ny,Nx).*NaN;
Ment    = zeros(Ny,Nx).*NaN;


for I=1:1:Ny,
   for J=1:1:Nx, 
      K      = find( Data(:,ColX)>=(X(J)-0.5.*StepX) & Data(:,ColX)<(X(J)+0.5.*StepX) & Data(:,ColY)>=(Y(I)-0.5.*StepY) & Data(:,ColY)<(Y(I)+0.5.*StepY));

      Mn(I,J)        = length(K);
      if (Mn(I,J)>MinPoints),
         Inonan         = find(isnan(Data(K,ColV))==0);
         if (length(Inonan)>0),
            Kn             = K(Inonan);
            Mnin(I,J)      = length(Inonan);
            Mmean(I,J)     = mean(Data(Kn,ColV));
            Mmedian(I,J)   = median(Data(Kn,ColV));
            [WM,WE]        = Util.stat.wmean(Data(Kn,[ColV,ColE]));
            if (isempty(WM))
                Mwmean(I,J)    = NaN;
            else
                Mwmean(I,J)    = WM;
            end
            if (isempty(WE))
                Mwerr(I,J)     = NaN;
            else
                Mwerr(I,J)     = WE;
            end
            Mstd(I,J)      = std(Data(Kn,ColV));
            Ment(I,J)      = -sum(Data(Kn,ColV).*log(Data(Kn,ColV)));
         end
      end
   end
end



