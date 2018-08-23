function [Long,Lat]=cosined2coo(CD1,CD2,CD3)
% Cosine direction to coordinates
% Package: celestial.coo
% Description: Convert cosine directions to coordinates in the same
%              reference frame. See also: cosined.m, coo2cosined.m
% Input  : - Matrix of first cosine directions.
%          - Matrix of second cosine directions.
%          - Matrix of third cosine directions.
% Output : - Matrix of longitudes [radians].
%          - Matrix of latitudes [radians].
% Tested : Matlab 7.10
%     By :  Eran O. Ofek                   Oct 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [RA,Dec]=celestial.coo.cosined2coo(0.1,0,1)
% Reliable: 1
%--------------------------------------------------------------------------

Long = atan2(CD2,CD1);
SLL  = sqrt(CD1.^2+CD2.^2);
I0   = SLL==0;
In0  = SLL~=0;
Lat  = zeros(size(Long));
Lat(In0) = atan(CD3(In0)./SLL(In0));
Lat(I0)  = sign(CD3(I0)).*pi./2;
% convert Long to 0..2pi range
InL      = Long<0;
Long(InL)= 2.*pi + Long(InL);
 
