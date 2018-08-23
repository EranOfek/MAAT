function Link=decals_viewer_link(RA,Dec,varargin)
% Get link to DECaLS image viewer
% Package: VO
% Description: Get link to DECaLS image viewer
% Input  : - J2000.0 R.A. [rad, sexagesimal, or [H M S]],
%            or object name to be resolved using NameServerFun.
%          - J2000.0 Dec. [rad, sexagesimal, or [sign D M S]].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'IsDeg' - If true then assume input coordinates is in degrees.
%                      Default is false.
%            'Layer' - Layer to plot. Default is 'decals-dr5'.
%            'Zoom'  - Zoom factor. Default is 13.
%            'NameServerFun' - Name server function handle.
%                      Default is @VO.name.server_ned.
%            'BaseURL' - Viewer base URL.
%                      Default is 'http://legacysurvey.org/viewer?'
% Output : - Cell array lof links to DECaLS image viewer.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Link=VO.DECaLS.decals_viewer_link(220./RAD,-0.1./RAD)
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;

if (nargin<2)
    Dec = [];
end

DefV.IsDeg                = false;
DefV.Layer                = 'decals-dr5';
DefV.Zoom                 = 13;
DefV.NameServerFun        = @VO.name.server_ned;
DefV.BaseURL              = 'http://legacysurvey.org/viewer?';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (isempty(Dec))
    % use name server
    [RA,Dec] = InPar.NameServerFun(RA,'r');
end
if (InPar.IsDeg)
    RA  = RA./RAD;
    Dec = Dec./RAD;
else
    RA  = celestial.coo.convertdms(RA,'gH','r');
    Dec = celestial.coo.convertdms(Dec,'gD','R');
end

RA  = RA.*RAD;
Dec = Dec.*RAD;

Nlink = numel(RA);
for Ilink=1:1:Nlink

    Key.RA    = sprintf('ra=%f',RA(Ilink));
    Key.Dec   = sprintf('dec=%f',Dec(Ilink));
    Key.Zoom  = sprintf('zoom=%d',InPar.Zoom);
    Key.Layer = sprintf('layer=%s',InPar.Layer);
    FN        = fieldnames(Key);
    Pars      = '';
    for Ikey=1:1:(numel(FN)-1)
        Pars = sprintf('%s%s&',Pars,Key.(FN{Ikey}));
    end
    Pars = sprintf('%s%s',Pars,Key.(FN{end}));

    Link{Ilink} = sprintf('%s%s',InPar.BaseURL,Pars);
end
%http://legacysurvey.org/viewer?ra=220.2876&dec=-0.3548&zoom=13&layer=decals-dr5