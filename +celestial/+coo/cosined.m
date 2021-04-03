function Res=cosined(Vec)
% Convert between coordinates and cosine directions
% Package: celestial.coo
% Description: Cosine direction transformation. Convert longitude
%              and latitude to cosine direction and visa versa.
%              See also: coo2cosined.m, cosined2coo.m
% Input  : - Column matrix, with 2 or 3 colums.
%            If 2 colums are given then [longitude latitude] in
%            radians, and the output is the cosine direction.
%            If three columns are given, then assume each line
%            contains cosine direction [X Y Z], and the output will
%            be the longitude and latitude in radians.
% Output : - Cosine direction or longitude and latitude [radian].
%            Each raw in the output corresponds to each raw in
%            the input.
% Tested : Matlab 7.0
% See also: cosined2coo, coo2cosined
%     By :  Eran O. Ofek                   Jul 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.coo.cosined(rand(5,3));    % Lon/lat from cosine-dir.
%          celestial.coo.cosined([0 0;1 1]);    % cosine-dir from long/lat
% Reliable: 1
%------------------------------------------------------------------------------
if (length(Vec(1,:))==2)
   Alpha = Vec(:,1);
   Delta = Vec(:,2);
   Res          = zeros(length(Vec(:,1)),3);
   Res(:,1)     = cos(Alpha).*cos(Delta);
   Res(:,2)     = sin(Alpha).*cos(Delta);
   Res(:,3)     = sin(Delta);
elseif (length(Vec(1,:))==3)
   L1           = Vec(:,1);
   L2           = Vec(:,2);
   L3           = Vec(:,3);
   Res          = zeros(length(Vec(:,1)),2);
   Res(:,1)     = atan2(L2,L1);                 % Alpha
   SLL          = sqrt(L1.^2+L2.^2);
   I0  = find(SLL==0);
   In0 = find(SLL~=0);
   Res(In0,2)     = atan(L3(In0)./SLL(In0));  % Delta
   Res(I0,2)      = sign(L3(I0)).*pi./2;  % Delta
else
   error('only 2/3 columns are allowed');
end
