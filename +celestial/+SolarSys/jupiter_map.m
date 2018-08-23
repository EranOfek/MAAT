function jupiter_map(Par)
%  Plot Jupiter image as observed from Earth at a given time
% Package: celestial.SolarSys
% Description: Plot Jupiter image as observed from Earth at a given time.
% Input  : - If two elements vector then:
%            [long_of_Sys_I, long_of_Sys_II]
%            else JD (scalar), or date vector (single date;
%            see julday.m for options). Date in TT time scale.
% Output : null
% Plot   : Jupiter RGB image refer to Jupiter System II, illuminated disk.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jan 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Needed : Map of Jupiter JupiterVoyagerMap.jpg (source: Voyager)
% Web example: http://astroclub.tau.ac.il/ephem/JovMap/
% Example: celestial.SolarSys.jupiter_map(2451545);
% Reliable: 2
%--------------------------------------------------------------------------
JupiterMapFile = 'JupiterVoyagerMap.jpg';

if (length(Par)==2),
   CMi = Par;
else
   if (length(Par)>1),
      JD = julday(Par);
   else
      JD = Par;
   end
   [CMg,CMi,Ds,De]=celestial.SolarSys.jup_meridian(JD);
end


%
% Assuming the red spot is at long. 109
%
SysI  = celestial.coo.angle_in2pi(CMi(1)+180,360); % 195
SysII = celestial.coo.angle_in2pi(CMi(2)+180,360); % 195


Im = imread(JupiterMapFile);
[SizeIm] = size(Im);
Scale = SizeIm(1)./180;

axesm ortho;
geoshow(Im, [Scale 90 SysII])
axis off;
set(gca,'YDir','Reverse');

