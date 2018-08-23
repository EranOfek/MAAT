function [Date,JD]=read_date(DateString,FormatCell,DateFormat)
%------------------------------------------------------------------------------
% read_date function                                                     ephem
% Description: Given a string or cell array of strings containing dates and
%              time in an arbitrary format, convert it to a matrix of dates
%              and a vector of JDs.
% Input  : - String or cell array of strings containing dates and time.
%          - Cell array containing three columns:
%            {start column, end column, format string}.
%            Example: {1,5,'%d'; 7,15,'%f'}.
%            Each one of the rows correspinds to one part of the date string.
%          - Description of each one (by order) of the parts of the date/time
%            string. This is a cell array in which each cell contain
%            a string corresponding to one of the parts of the string format.
%            Options are: 'Year', 'MonthName', 'MonthNameShort', 'Month', 'Day'
%                         'Frac', 'Hour', 'Min', 'Sec'.
%            Example: {'Year','MonthNameShort','Day','Hour','Min'}.
% Output : - Matrix of dates [D M Y FracDay]
%          - Vector of JD.
% Example: [Date,JD]=read_date(Data{1},...
%          {1,4,'%d'; 6,8,'%s'; 10, 11, '%d'; 13, 14, '%d'; 16,17,'%d'},...
%          {'Year','MonthNameShort','Day','Hour','Min'});
%          % will read string with the format: '2000-Feb-01 00:00'.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    August 2008
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------------------
Month = {'January';'Febuary';'March';'April';'May';'June';'July';'August';'September';'October';'November';'December'};
MonthShort = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};

if (isstr(DateString)==1)
   % convert string to a single cell of string
   DateString{1} = DateString;
else
   % do nothing - already in a cell format
end

Nd = length(DateString);
Ncol = size(FormatCell,1);
for Id=1:1:Nd
   Line = DateString{Id};

   if (length(Line)==0)
      % ignore line
   else
      for Icol=1:1:Ncol
         switch FormatCell{Icol,3}
          case {'%s','%c'}
             Data{Id,Icol} = sscanf(Line(FormatCell{Icol,1}:FormatCell{Icol,2}),FormatCell{Icol,3});
          otherwise
             Data{Id,Icol} = Util.string.str2num_nan(Line(FormatCell{Icol,1}:FormatCell{Icol,2}));
         end
      end
   end
end

Nf = length(DateFormat);
Hour = zeros(Nd,1);
Min  = zeros(Nd,1);
Sec  = zeros(Nd,1);

for If=1:1:Nf
   switch lower(DateFormat{If})
    case 'year'
       Year = [Data{:,If}].';
    case 'month'
       Month = [Data{:,If}].';
    case 'monthname'
       Month = zeros(Nd,1);
       for Id=1:1:Nd
          I=find(Util.cell.isempty_cell(strfind(Month,Data{Id,If}))==0);
          if (isempty(I)==1)
             error('Illegal short month name');
          else
	     Month(Id) = I;
          end
       end
    case 'monthnameshort'
       Month = zeros(Nd,1);
       for Id=1:1:Nd
          I=find(Util.cell.isempty_cell(strfind(MonthShort,Data{Id,If}))==0);
          if (isempty(I)==1)
             error('Illegal short month name');
          else
	     Month(Id) = I;
          end
       end
    case 'day'
       Day = [Data{:,If}].';
    case 'hour'
       Hour = [Data{:,If}].';
    case 'min'
       Min = [Data{:,If}].';
    case 'sec'
       Sec = [Data{:,If}].';
    case {'frac','fracday'}
       FracDay = [Data{:,If}].';
    otherwise
       error('Unknown date/time type option')
   end
end

Date = [Day Month Year Hour Min Sec];
JD   = celestial.time.julday(Date).';
