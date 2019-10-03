function HelioEcLong=ecliptic2helioecliptic(Long,Date,Algo)
% Ecliptic longitude to Helio-ecliptic longitude
% Package: celestial.coo
% Description: Transform ecliptic longitude to Helio-ecliptic longitude.
% Input  : - Ecliptic longitude [radians]. See convert_dms.m for
%            additional options.
%          - Date [Day Mounth Year FracDay], or JD column vector.
%            See julday.m for options.
%          - Solar longitude algorithm:
%            'low' - low accuracy (0.01 deg). Default.
% Output : - Helio-ecliptic longitude [radians].
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: HelioEcLong=celestial.coo.ecliptic2helioecliptic(1,[21 3 2000])
% Reliable: 2
%--------------------------------------------------------------------------

Def.Algo = 'low';
if (nargin==2)
    Algo = Def.Algo;
end

Long = celestial.coo.convertdms(Long,'gH','r');

if (size(Date,2)==1)
    JD = Date;
else
    JD   = celestial.time.julday(Date);
end

switch lower(Algo)
    case 'low'
        [~,~,~,SL]=celestial.SolarSys.suncoo(JD,'j');
    otherwise
        error('Unknown Algo option');
end

HelioEcLong = Long - SL;
HelioEcLong = celestial.coo.angle_in2pi(HelioEcLong);
        
        