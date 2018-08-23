function [EasterDate,EasterJD]=easter_date(Y)
% Calculate the date of Easter
% Package: celestial.time
% Description: Calculate the date of Easter for any Gregorian year.
% Input  : - Year (integer).
% Output : - Date of Easter, [D M Y].
%          - Date of Easter, [JD].
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: Oudin (1940)
% Example: easter_date(2010);  % will return [4 4 2010].
% Reliable: 2
%--------------------------------------------------------------------------
C  = fix(Y./100);
N  = Y - 19.*fix(Y./19);
K  = fix((C - 17)./25);
I  = C - fix(C./4) - fix((C - K)./3) + 19.*N + 15;
I  = I - 30.*fix(I./30);
I  = I - fix(I./28).*(1 - fix(I./28).*fix(29./(I+1)).*fix((21 - N)./11));
J  = Y + fix(Y./4) + I + 2 - C + fix(C./4);
J  = J - 7.*fix(J./7);
L  = I - J;
M  = 3 + fix((L + 40)./44);
D  = L + 28 - 31.*fix(M./4);

EasterDate = [D M Y];
EasterJD   = convert.date2jd(EasterDate);

 
