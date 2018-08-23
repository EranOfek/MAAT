function [CosB,Range]=cosbell(PrecentFlat,Range)
% Cosine bell function (Obsolete: use timeseries.taper instead)
% Package: timeseries
% Description: cosine bell function
%                    Generating cosine bell function in the range
%                    [Start:End] with its inner PercentFlat part
%                    as flat function.
% Input  : - The precentage of flat part of the cosine bell.
%            (in the range 0..1).
%          - Column vector of independent variable,
%            default is [0:0.01:1]'.
% Output : - Column vector of cosine bell value.
%          - The Range (independent variable) column vector.
% Tested : Matlab 5.0
%     By : Eran O. Ofek                    Dec 1997
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable : 2
%--------------------------------------------------------------------

if nargin==1,
   Range = [0:0.01:1]';
elseif nargin>2,
   error('Number of input arguments should be 1 or 2');
end

[N,M]   = size(Range);
Start   = min(Range);
End     = max(Range);
Total   = End - Start;
TotalCB = Total.*(1 - PrecentFlat).*0.5;
StartCB = Start + Total.*(1 - PrecentFlat).*0.5;
EndCB   = End   - Total.*(1 - PrecentFlat).*0.5;

CosB    = zeros(N,1);
for Ind=1:1:N,
   if Range(Ind)<StartCB,
      CosB(Ind) = 0.5+0.5.*cos(pi.*(Range(Ind)-Start)./TotalCB + pi);
   elseif Range(Ind)>EndCB,
      CosB(Ind) = 0.5+0.5.*cos(pi.*(Range(Ind)-EndCB)./TotalCB);
   else
      CosB(Ind) = 1;
   end
end
