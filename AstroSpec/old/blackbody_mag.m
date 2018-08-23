function [MagTemp]=blackbody_mag(Temp,Filter,System,Radius,Dist);
%-------------------------------------------------------------------------
% blackbody_mag function                                        AstroSpec
% Description: Calculate the magnitude, in a given bad, of a black body
%              given its temperature, radius and distance.
%              Work only in the range 100-10^7 K.
%              Accuracy is better than 0.015 mag for temperature range.
%              For other temperatures and additional filters
%              use: blackbody_mag_c.m
% Input  : - Vector of blackbody temperature.
%          - Filter name: 12FNUBVRIJHKugriz
%            for UW1, UW2, FUV, NUV,...,Johnson, Cousins, 2MASS, SDSS...
%          - Magnitude system: 'V' for Vega, 'A' for 'AB'.
%          - blackbody sphere radius [cm].
%          - blackbody distance from observer [pc]. 
% Output : - Vector of magnitude in each temperature.
% Tested : Matlab 7.0
%     By : Eran O. Ofek           June 2005
%    URL : http://wise-ftp.tau.ac.il/~eran/matlab.html
% Notes  : The program uses a mat file named: data_blackbody_mag.mat 
% See Also: blackbody_mag_c.m
%-------------------------------------------------------------------------
Pc           = get_constant('pc','cgs');

InterpMethod = 'linear';
%FilterList   = '12FNUBVRIJHKugriz';   % F - FUV ; N - NUV
%                                      % 1 - UW1 ; 2 - UW2


switch lower(System)
 case {'v','vega'}
    load BB_Data_Vega.mat
    MagList = BB_Mag_Vega;
    clear BB_Mag_Vega;
 case {'a','ab'}
    load BB_Data_AB.mat
    MagList = BB_Mag_AB;
    clear BB_Mag_AB;
 otherwise
    error('Unknown System option');
end

FilterColInd = findstr(FilterList,Filter);

% T is loaded from mat file

% magnitude in filter as function of T
MagInFilterT = MagList(:,FilterColInd);

MagTemp = interp1(T,MagInFilterT,Temp,InterpMethod);

%MagTemp = MagTemp - 2.5.*log10(4.*pi);   % convert to bb with radius of 1cm
MagTemp = MagTemp - 5.0.*log10(Radius);  % convert to radius
MagTemp = MagTemp + 5.0.*log10(Dist.*Pc);
%MagTemp = MagTemp + 20;    % 20 = 2.5.*log10(1e8)  Ang to cm





