function MagSDSS=ps1_2_sdss_mag(Band,MagPS1,gPS1,iPS1)
% PS1 magnitudes to PS1 magnitude using the Finkbeiner+2015 relations
% Package: VO
% Description: Convert PS1 magnitudes to PS1 magnitude using the
%              Finkbeiner+2015 relations.
% Input  : - PS1 Band name, 'g','r','i','z','y'
%          - PS1 magnitude
%          - PS1 g-band magnitude.
%          - PS1 i-band magnitude.
% Output : - SDSS magnitude in band.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: MagSDSS=VO.PS1.ps1_2_sdss_mag('g',20,20,21)
% Reliable: 2
%--------------------------------------------------------------------------

%     a0        a1        a2       a3
u = [ 0.04438 -2.26095 -0.13387  0.27099];
g = [-0.01808 -0.13595  0.01941 -0.00183];
r = [-0.01836 -0.03577  0.02612 -0.00558];
i = [ 0.01170 -0.00400  0.00066 -0.00058];
z = [-0.01062  0.07529 -0.03592  0.00890];
y = [ 0.08924 -0.20878  0.10360 -0.02441];

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

switch lower(Band)
    case 'g'
        a = g;
    case 'r'
        a = r;
    case 'i'
        a = i;
    case 'z'
        a = z;
    case 'y'
        a = y;
    otherwise
        error('Unknown PS1 band');
end

        
x = gPS1 - iPS1;
MagSDSS = MagPS1 - a(1) -a(2).*x -a(3).*x.^2 -a(4).*x.^3;