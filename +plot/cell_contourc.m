function Levels=cell_contourc(varargin);
%--------------------------------------------------------------------------
% cell_contourc function                                          plotting
% Description: A contourc-like program that return a cell array of the
%              contour levels.
% Input  : * See contourc for available input arguments.
% Output : - Structure containing two elements:
%            .L - containing the contour levels
%            .C - Cell array, for each level, containing [X,Y] for the
%                 contours.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                   January 2008
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%--------------------------------------------------------------------------

CS = contourc(varargin{:});
N        = size(CS,2);

IndStart = 2;
Number   = CS(2,1);
IndEnd   = IndStart + Number - 1; 
Level    = CS(1,1);

Levels.L{1} = Level;
Levels.C{1} = [CS(1,IndStart:IndEnd).', CS(2,IndStart:IndEnd).'];

I = 1;

while (IndEnd<N),
   I = I + 1;
   IndStart = IndEnd + 2;
   Number   = CS(2,IndEnd+1);
   IndEnd   = IndStart + Number - 1; 
   Level    = CS(1,IndStart-1);

   Levels.L{I} = Level;
   Levels.C{I} = [CS(1,IndStart:IndEnd).', CS(2,IndStart:IndEnd).'];   
end
