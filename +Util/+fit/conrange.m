function [MinVal,MinX,LowerX,UpperX,LowerE,UpperE]=conrange(X,Y,Offsets)
%--------------------------------------------------------------------------
% conrange function                                                 FitFun
% Description: Given two vectors of Y and X, calculate the range
%              in X that satisfies the constrain Y<=(min(Y)+Offset).
%              This is useful for calculating likelihood/chi^2 errors for
%              1-D data.
% Input  : - X vector.
%          - Y vector, which is Y(X).
%          - Vector of Offsets for which to calculate the ranges.
% Output : - min(Y).
%          - X at min(Y).
%          - min(X) that satisfy the constrain: Y<=(Min+Offset)
%          - max(X) that satisfy the constrain: Y<=(Min+Offset)
%          - Distance between min(X) that satisfy the constrain
%            Y<=(Min+Offset) and the X at min(Y).
%            This is the lower error.
%          - Distance between min(X) that satisfy the constrain
%            Y<=(Min+Offset) and the X at max(Y).
%            This is the upper error.
% Tested : Matlab 7.0
%     By : Eran O Ofek                     Feb 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X = [-100:1:100].';
%          Y = X.^2;
%          [MinVal,MinX,LowerX,UpperX,LowerE,UpperE]=conrange(X,Y,[1;2;100]);
% Reliable: 2
%--------------------------------------------------------------------------

Data = sortrows([X,Y]);
X    = Data(:,1);
Y    = Data(:,2);

[MinVal,MinInd] = min(Y);
MinX = X(MinInd);



% for each offset calculate bounderies
No = length(Offsets);
LowerX = zeros(No,1);
UpperX = zeros(No,1);
LowerE = zeros(No,1);
UpperE = zeros(No,1);
for Io=1:1:No,
   I = find( Y<=(MinVal+Offsets(Io)) );

   LowerI = min(I);
   UpperI = max(I);
   
   LowerX(Io) = X(LowerI);
   UpperX(Io) = X(UpperI);

   LowerE(Io) = MinX - LowerX(Io);
   UpperE(Io) = UpperX(Io) - MinX;
   
end
