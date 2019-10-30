function [L,n]=refraction_index(varargin)
% Return the refraction index as a function of wavelength
% Package: telescope.Optics
% Description: Return the refraction index as a function of wavelength
%              for various materials.
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'T' - Temperature [C]. Default is 20.
%            'Material' - Options are:
%                         'SiO2' (fused silica) default.
% Output : - Wavelength [Ang].
%          - Refraction index
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [L,n]=telescope.Optics.refraction_index('Material','SiO2')
% Reliable: 
%--------------------------------------------------------------------------



DefV.T                    = 20;
DefV.Material             = 'SiO2';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


switch lower(InPar.Material)
    case 'sio2'
        L = [0.21:0.01:6.7]';
        n = sqrt(1 + 0.6961663.*L.^2./(L.^2 - 0.0684043.^2) + ...
                0.4079426.*L.^2./(L.^2 - 0.1162414.^2) + ...
                0.8974794.*L.^2./(L.^2 - 9.896161.^2));
        L         = L.*1e4;
        WaveRange = [210, 67000];  % [Ang]
        T         = 20;  % [C]
        Reference = 'https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson';
    otherwise
        error('Unknown Material option');
end

if (InPar.T~=T)
    error('Request for illigal temperature');
end
