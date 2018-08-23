function [FullName,ShortName]=month_name(Month)
% Convert month number to name
% Package: celestial.time
% Description: Given a month number return a string with a month name.
% Input  : - Vector of month number.
% Output : - Cell vector of month full name (9 chars).
%          - Cell vector of month short name (3 chars).
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jan 2003
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [FullName,ShortName]=celestial.time.month_name([1;2])
% Reliable: 1
%--------------------------------------------------------------------------

FullName  = cell(0,0);
ShortName = cell(0,0);

for I=1:1:length(Month),
   switch Month(I)
    case 1
       FName = 'January  ';
       SName = 'Jan';
    case 2
       FName = 'Febuary  ';
       SName = 'Feb';
    case 3
       FName = 'March    ';
       SName = 'Mar';
    case 4
       FName = 'April    ';
       SName = 'Apr';
    case 5
       FName = 'May      ';
       SName = 'May';
    case 6
       FName = 'June     ';
       SName = 'Jun';
    case 7
       FName = 'July     ';
       SName = 'Jul';
    case 8
       FName = 'August   ';
       SName = 'Aug';
    case 9
       FName = 'September';
       SName = 'Sep';
    case 10
       FName = 'October  ';
       SName = 'Oct';
    case 11
       FName = 'November ';
       SName = 'Nov';
    case 12
       FName = 'December ';
       SName = 'Dec';
    otherwise
       error('Illegal month number');
   end
   FullName{I,1}  = FName;
   ShortName{I,1} = SName;
end


