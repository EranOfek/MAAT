function Flag=in_box(Lon,Lat,BoxCenter,BoxHalfSize)
% Check if celestial coordinates are in a box.
% Package: celestial
% Description: Check if celestial coordinates are in a box defined by four
%              corners and its sides are great circles.
%              See also: celestial.htm.in_polysphere
% Input  : - Vector of longitude to check [rad].
%          - Vector of latitude to check [rad].
%          - [Longitude,Latitiude] of Box center [rad].
%          - [HalfSizeLon, HalfSizeLat] of box [rad].
% Output : - A vector of logical flags indicating if coordinates is in box.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Flag=celestial.coo.in_box(rand(100,1),rand(100,1),[0 0],[0.5 0.5])
% Reliable: 
%--------------------------------------------------------------------------



BoxCenter   = BoxCenter(:).';
BoxHalfSize = BoxHalfSize(:).';

Corners = [BoxCenter(1) - BoxHalfSize(1), BoxCenter(2) - BoxHalfSize(2); ...
           BoxCenter(1) + BoxHalfSize(1), BoxCenter(2) - BoxHalfSize(2); ...
           BoxCenter(1) + BoxHalfSize(1), BoxCenter(2) + BoxHalfSize(2); ...
           BoxCenter(1) - BoxHalfSize(1), BoxCenter(2) + BoxHalfSize(2)];
           

Flag = celestial.htm.in_polysphere([Lon,Lat],Corners);