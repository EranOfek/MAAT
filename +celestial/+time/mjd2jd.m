function JD=mjd2jd(MJD)
% Convert MJD to JD
% Package: celestial.time
% Description: Convert MJD to JD.
%              See also: convert.time
% Input  : - MJD
% Output : - JD (i.e., MJD + 2400000.5)
% See also: julday.m, jd2date.m, jd2mjd.m
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Dec 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.time.mjd2jd(2450000)
% Reliable: 1
%--------------------------------------------------------------------------

JD = MJD + 2400000.5;
