function JD=year2jd(Year,Type)
% Convert year to JD
% Package: celestial.time
% Description: Return the Julian day at Jan 1 st of a given list of years.
%              See also: convert.time instead.
% Input  : - Column vector of years.
%          - Year type. Options are:
%            '1'  - January 1st of the year. Default
%            'by','b' - Bessilian year.
%            'jy','year,'yr',j' - Julian year.
%            'jd' - Input is julian days - return input as is.
% Output : - Vector of JDs.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: JD=celestial.time.year2jd(2000);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<2),
    Type = '1';
end

Ny = numel(Year);

switch lower(Type)
    case '1'
        % year at jan 1st
        JD = celestial.time.julday([ones(Ny,2), Year]);
    case {'jy','j','year','yr'}
        % Julian year
        JulYear  = 365.25;
        JD0      = celestial.time.julday([1 1 2000 0]);
        JD       = (Year - 2000).*JulYear + JD0;
        
    case {'by','b'}
        % Bessilian year
        BesYear  = 365.242189;
        JD0      = 2451544.53;
        JD       = (Year - 2000).*BesYear + JD0;
    case 'jd'
        % Input is already in JD
        JD = Year;
        
    otherwise
        error('Unknown Type option');
end
