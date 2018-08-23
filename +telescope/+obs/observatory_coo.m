function [ObsCoo,IDstring,FullName]=observatory_coo(ObsName)
% Geodetic coordinates of selected observatories
% Package: telescope.obs
% observatory_coo function                                           ephem
% Description: Return geodetic coordinates of an observatory.
% Input  : - Observatory name:
%            'Wise' - wise observatory        [  34.763, 30.596]/RAD.
%            'Kraar'- Kraar observatory (WIS) [34.812711,31.908083]./RAD
%            'KPNO' - Keat Peak National obs. [-111.60 , 31.980]/RAD.
%            'Keck' - Keck observatory.       [-155.478, 19.828]/RAD.
%            'APO'  - Apache Point obs.       [-105.82 , 32.780]/RAD.
%            'CA'   - Calar-Alto obs.         [   2.546, 37.224]/RAD.
%            'MMT'  - MMT obs. (Mt. Hopkins)  [-110.885, 31.688]/RAD.
%            'Paran'- ESO Paranal obs.        [ -70.403,-24.625]/RAD.
%            'LaPal'- La Palma Island obs.    [ -18.881, 28.760]/RAD.
%            'SAO'  - Special Astrop. obs.    [  41.442, 43.653]/RAD.
%            'Palom'- Palomar observatory     [-116.863, 33.357]/RAD.
%            'AAO'  - Anglo-Australian obs.   [ 149.067,-31.277]/RAD.
%            'CTIO' - Cerro Tololo Inter-Am.  [ -70.815,-30.165]/RAD.
%            'ESO'  - ESO La Silla obs.       [ -70.730,-29.257]/RAD.
%            'Camp' - Las Campanas obs.       [ -70.700,-29.008]/RAD.
%            'Strom'- Mount Stromlo obs.      [ 149.008,-35.320]/RAD.
%            'Lowel'- Lowell obs. Flagstaff   [-111.665, 35.203]/RAD.
% Output : - Observatory geodetic position [Longitude, Latitude]
%            in radians.
%          - Matrix of ID strings. Line per observatory.
%          - Matrix of obs. names strings. Line per observatory.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [ObsCoo]=telescope.obs.observatory_coo('ESO');
% Reliable: 2
%--------------------------------------------------------------------------
RAD      = 180./pi;

switch lower(ObsName)
    case 'wise'
        ObsLon = [  34.763]./RAD;
        ObsLat = [  30.596]./RAD;
    case 'kraar'
        ObsLon = 34.812711./RAD;
        ObsLat = 31.908083./RAD;
    case 'kpno'
        ObsLon = [-111.60 ]./RAD;
        ObsLat = [  31.980]./RAD;
    case 'keck'
        ObsLon = [-155.478]./RAD;
        ObsLat = [  19.828]./RAD;
    case 'apo'
        ObsLon = [-105.82 ]./RAD;
        ObsLat = [  32.780]./RAD;
    case 'ca'
        ObsLon = [   2.546]./RAD;
        ObsLat = [  37.224]./RAD;
     case 'mmt' 
        % MMT obs. (Mt. Hopkins)
        ObsLon = [-110.885]./RAD;
        ObsLat = [  31.688]./RAD;
    case 'paran'
         % ESO Paranal obs.   
        ObsLon = [ -70.403]./RAD;
        ObsLat = [ -24.625]./RAD;
    case 'lapal'
        % La Palma Island obs.
        ObsLon = [ -18.881]./RAD;
        ObsLat = [  28.760]./RAD;
     case 'sao'
        % Special Astrop. obs.
        ObsLon = [  41.442]./RAD;
        ObsLat = [  43.653]./RAD;
     case 'palom'
        % Palomar observatory
        ObsLon = [-116.863]./RAD;
        ObsLat = [  33.357]./RAD;
     case 'aao'
        % Anglo-Australian obs.
        ObsLon = [ 149.067]./RAD;
        ObsLat = [ -31.277]./RAD;
     case 'ctio'
        % Cerro Tololo Inter-Am.
        ObsLon = [ -70.815]./RAD;
        ObsLat = [ -30.165]./RAD;
     case 'eso'
        % ESO La Silla obs.
        ObsLon = [ -70.730]./RAD;
        ObsLat = [ -29.257]./RAD;
     case 'camp'
        % Las Campanas obs.
        ObsLon = [ -70.700]./RAD;
        ObsLat = [ -29.008]./RAD;
     case 'strom'
        %Mount Stromlo obs.
        ObsLon = [ 149.008]./RAD;
        ObsLat = [-35.320]./RAD;
     case 'lowel'
        % Lowell obs. Flagstaff
        ObsLon = [-111.665]./RAD;
        ObsLat = [  35.203]./RAD;
     otherwise
        error('Unkown observatory name');
end

ObsCoo = [ObsLon, ObsLat];

if (nargout>1)
   IDstring = strvcat('Wise',...
                      'KPNO',...
                      'Keck',...
                      'APO',...
                      'CA',...
                      'MMT',...
                      'Paran',...
                      'LaPal',...
                      'SAO',...
                      'Palom',...
                      'AAO',...
                      'CTIO',...
                      'ESO',...
                      'Camp',...
                      'Strom',...
                      'Lowel');

   FullName = strvcat('Wise observatory',...
                      'KPNO',...
                      'Keck observatory',...
                      'APO',...
                      'Calar-Alto',...
                      'MMT observatory',...
                      'ESO Paranal',...
                      'La Palma Island obs.',...
                      'Special Astrop. obs.',...
                      'Palomar observatory',...
                      'Anglo-Australian obs.',...
                      'Cerro Tololo Inter-Am.',...
                      'ESO La Silla obs.',...
                      'Las Campanas obs.',...
                      'Mount Stromlo obs',...
                      'Lowell obs. Flagstaff');

end
