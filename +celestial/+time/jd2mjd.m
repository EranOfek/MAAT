function MJD=jd2mjd(JD)
% Convert JD to MJD
% Package: celestial.time
% Description: Convert JD to MJD. See also: convert.time
% Input  : - JD
% Output : - MJD (i.e., JD - 2400000.5)
% See also: julday.m, jd2date.m, mjd2jd.m
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Dec 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.time.jd2mjd(2450000)
% Reliable: 1
%---------------------------------------------------------------------------

MJD = JD - 2400000.5;
