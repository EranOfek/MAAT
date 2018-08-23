function JD=read_mpc_packed_epoch(PackFormat)
% Convert the MPC packed date format to JD
% Package: celestial.SolarSys
% Description: Convert the MPC packed date format to JD.
% Input  : - Cell array in which each cell containing string of packed
%            date.
% Output : - Vector of JD.
% Reference: http://www.cfa.harvard.edu/iau/info/PackedDates.html
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jan 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

N = length(PackFormat);
JD = zeros(N,1);
for I=1:1:N,
   if (length(PackFormat{I})<5),
      %
      JD(I) = NaN;
   else
      switch PackFormat{I}(1)
       case 'I'
          Year = 1800;
       case 'J'
          Year = 1900;
       case 'K'
          Year = 2000;
       otherwise
          error('Unknown PackFormat');
      end
   
      Year = Year + str2num(PackFormat{I}(2:3));
      switch PackFormat{I}(4)
       case '1'
          Month = 1;
       case '2'
          Month = 2;
       case '3'
          Month = 3;
       case '4'
          Month = 4;
       case '5'
          Month = 5;
       case '6'
          Month = 6;
       case '7'
          Month = 7;
       case '8'
          Month = 8;
       case '9'
          Month = 9;
       case 'A'
          Month = 10;
       case 'B'
          Month = 11;
       case 'C'
          Month = 12;
       otherwise
          error('Unknown PackFormat');
      end   
   
      switch PackFormat{I}(5)
       case '1'
          Day  = 1;
       case '2'
          Day  = 2;
       case '3'
          Day  = 3;
       case '4'
          Day  = 4;
       case '5'
          Day  = 5;
       case '6'
          Day  = 6;
       case '7'
          Day  = 7;
       case '8'
          Day  = 8;
       case '9'
          Day  = 9;
       case 'A'
          Day  = 10;
       case 'B'
          Day  = 11;
       case 'C'
          Day  = 12;
       case 'D'
          Day  = 13;
       case 'E'
          Day  = 14;
       case 'F'
          Day  = 15;
       case 'G'
          Day  = 16;
       case 'H'
          Day  = 17;
       case 'I'
          Day  = 18;
       case 'J'
          Day  = 19;
       case 'K'
          Day  = 20;
       case 'L'
          Day  = 21;
       case 'M'
          Day  = 22;
       case 'N'
          Day  = 23;
       case 'O'
          Day  = 24;
       case 'P'
          Day  = 25;
       case 'Q'
          Day  = 26;
       case 'R'
          Day  = 27;
       case 'S'
          Day  = 28;
       case 'T'
          Day  = 29;
       case 'U'
          Day  = 30;
       case 'V'
          Day  = 31;
       otherwise
          error('Unknown PackFormat');
      end
      JD(I) = celestial.time.julday([Day Month Year 0]);
   end
end
