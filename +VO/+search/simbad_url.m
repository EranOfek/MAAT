function [Data]=simbad_url(RA,Dec,varargin)
% Generate a SIMBAD URL for coordinates
% Package: VO.search
% Description: Generate a SIMBAD URL for coordinates
% Input  : - J2000 RA [rad]. See celestial.coo.convertdms for other options.
%            This is either a scalar or a vector of coordinates.
%          - J2000 Dec [rad]. See celestial.coo.convertdms for other options.
%            This is either a scalar or a vector of coordinates.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SearchRadius' - Default is 2.
%            'SearchRadiusUnits' - Default is 'arcmin'.
%            'CooFrame' - Default is 'FK5'.
%            'CooEpoch' - Default is 2000.
%            'CooEquinox' - Default is 2000.
% Output : - A structure arrary withe the following fields:
%            'URL' - The SIMBAD URL for this coordinates.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Data]=VO.search.simbad_url(1,1);
% Reliable: 2
%--------------------------------------------------------------------------


DefV.SearchRadius         = 2;
DefV.SearchRadiusUnits    = 'arcmin';
DefV.CooFrame             = 'FK5';
DefV.CooEpoch             = 2000;
DefV.CooEquinox           = 2000;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

%http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=01%3A47%3A02.707+%2B59%3A36%3A22.57&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=2&Radius.unit=arcmin&submit=submit+query&CoordList=
%http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=145+%2B12&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=2&Radius.unit=arcmin&submit=submit+query&CoordList=

RA  = celestial.coo.convertdms(RA,'gH','d');
Dec = celestial.coo.convertdms(Dec,'gD','d');

Ncoo = numel(RA);

for Icoo=1:1:Ncoo
    Data(Icoo).URL = sprintf('http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=%f+%f&CooFrame=%s&CooEpoch=%d&CooEqui=%d&CooDefinedFrames=none&Radius=%f&Radius.unit=%s&submit=submit+query&CoordList=',...
                              RA,Dec,InPar.CooFrame,InPar.CooEpoch,InPar.CooEquinox,...
                              InPar.SearchRadius,InPar.SearchRadiusUnits);
    %Data(Icoo).URL = urlencode(Data(Icoo).URL);
end

    