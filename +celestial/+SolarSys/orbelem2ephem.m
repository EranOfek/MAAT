function []=orbelem2ephem(JD_TT,OrbElem,varargin)
% SHORT DESCRIPTION HERE
% Package: celestial
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------
RAD = 180./pi;

SEC_IN_DAY = 86400;
c = constant.c./constant.au .* SEC_IN_DAY;   % speed of light [au/day]

DefV.Observer             = 'Earth';
DefV.ThreshLightTime      = 1e-5;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


%Orbital elements : [q e T i \Omega \omega], in au, days, radians.
OrbElem = [2.544709153978707, .07987906346370539, 2453193.6614275328, ...
           [10.58671483589909, 80.40846590069125,   73.1893463033331]./RAD];
JD_TT = celestial.time.julday([1 1 2005; 12 8 2019; 16 8 2019]);       
       
%  EC= .1155130281765817   QR= 2.364480109770418   TP= 2456873.8451829329      
%  OM= 121.6635322632723   W=  83.08494011660611   IN= 13.51502086244268  



% get Observer ephemerides
% currently default is Earth
% output is in [au], [au/day]
[ObserverCoo,ObserverVel] = celestial.SolarSys.calc_vsop87(JD_TT, InPar.Observer, 'e', 'E');
ObserverCoo = ObserverCoo.';
ObserverVel = ObserverVel.';

ConvergedLightTime = false;
LightTime = 0;
Ilt = 0;
while ~ConvergedLightTime
    Ilt = Ilt + 1;
    % coordinates of object
    [ObjectCoo,ObjectVel,Nu]=celestial.Kepler.elements2position(JD_TT-LightTime,OrbElem);

    % Coordinates of objects relative to observer
    Coo   = ObjectCoo - ObserverCoo;
    % Object-Observer distance
    Delta = sqrt(sum(Coo.^2,2));
    % light time correction
    PrevLightTime = LightTime;
    LightTime     = Delta./c;
    if (abs(LightTime-PrevLightTime)<InPar.ThreshLightTime)
        ConvergedLightTime = true;
    end
    %ConvergedLightTime = true;
end

% calculate RA/Dec
EcLon = atan2(Coo(:,2),Coo(:,1));
EcLat = asin(Coo(:,3)./sqrt( sum(Coo.^2,2) ));
EqCoo = celestial.coo.coco([EcLon, EcLat], 'e','j2000');

    



