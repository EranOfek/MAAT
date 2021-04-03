function OutCat=nearest_bright_star(RA,Dec,Mag,N)
% Select N bright stars near a given coordinate
% Package: +celestial.stars
% Input  : - J2000.0 RA [deg] or sexagesimal string.
%          - J2000.0 Dec [deg] or sexagesimal string.
%          - Mag. limit. Default is 2.
%          - Number of stars (sorted by distance). Default is 3.
% Example : celestial.stars.nearest_bright_star(1,1,2)
%      By : Eran Ofek              Mar 2021

if nargin<4
    N = 3;
    if nargin<3
        Maf = 2;
    end
end


RAD = 180./pi;

ColRA  = 1;
ColDec = 2;
ColMag = 3;

Cat = cats.bright.mag6;

Flag = Cat(:,3)<Mag;
Cat  = Cat(Flag,:);

if ~ischar(RA)
    RA = RA./RAD;  % convert deg to rad
end
if ~ischar(Dec)
    Dec = Dec./RAD;  % convert deg to rad
end


D = celestial.coo.sphere_dist(RA,Dec,Cat(:,ColRA),Cat(:,ColDec));


[~,SI] = sort(D);

OutCat = Cat(SI(1:1:N),[ColRA, ColDec, ColMag]);
OutCat(:,[ColRA, ColDec]) = OutCat(:,[ColRA, ColDec]).*RAD;
