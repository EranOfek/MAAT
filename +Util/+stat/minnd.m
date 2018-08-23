function [Min,MinIndVec]=minnd(Data);
%------------------------------------------------------------------------------
% minnd function                                                     AstroStat
% Description: Return the global minimum of a N-D matrix and its indices.
%              This is equivalent to min(Data(:)), but it also returns the
%              indices of the global minimum.
% Input  : - N-D matrix
% Output : - Value of global Minimum
%          - Vector that contains the indices of the global minimum
%            in the N-D matrix.
% Tested : Matlab 6.5
%     By : Eran O. Ofek                       Feb 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [M,I]=minnd(rand(4,4,4,4))
% Reliable: 1
%------------------------------------------------------------------------------

SizeData = size(Data);
Vec = reshape(Data, [prod(SizeData), 1]);
[Min, MinInd] = min(Vec);

ND = length(SizeData);

StrCom = '[I1 ';
for Idim=2:1:ND,
   StrTemp = sprintf('I%d ',Idim);
   StrCom  = [StrCom,',',StrTemp];
end
StrInd = [StrCom,']'];
StrCom = [StrInd,'=ind2sub(SizeData,MinInd);'];
eval(StrCom);
MinIndVec=eval(StrInd);

%[I1, I2, I3] = ind2sub(SizeData,MinInd)


   
