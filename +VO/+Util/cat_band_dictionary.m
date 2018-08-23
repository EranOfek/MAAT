function [CatBand,CatBandErr,CatColor,CatColorErr]=cat_band_dictionary(Cat,Band)
% Return the band (filter) name in a given catalog.
% Package: VO.Util
% Description: Given a catalog name (e.g., 'sdss') and a band name
%              (e.g., 'r'), return the corresponding band in the catalog 
%              (e.g., 'modelMag_g') its error and color for photometric
%              calibration (e.g., 'modelMag_g - modelMag_r').
% Input  : - A catalog name (e.g., 'sdss').
%            Available catalog names: 'sdss'
%          - A band name (e.g., 'r').
% Output : - A string of corresponding band name in the catalog.
%          - A string of corresponding band error name in the catalog.
%          - A string of color formula for photometric calibration.
%          - A string of color error formula for photometric calibration.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [CatBand,CatBandErr,CatColor,CatColorErr]=VO.Util.cat_band_dictionary('sdss','g')
% Reliable: 2
%--------------------------------------------------------------------------

Icat = 0;
Icat = Icat + 1;
Dictionary(Icat).Cat          = 'SDSSDR10';
Dictionary(Icat).MasterBand   = {'u',            'g',            'r',            'i',            'z'};
Dictionary(Icat).CatBand      = {'modelMag_u',   'modelMag_g',   'modelMag_r',   'modelMag_i',   'modelMag_z'};
Dictionary(Icat).CatBandErr   = {'modelMagErr_u','modelMagErr_g','modelMagErr_r','modelMagErr_i','modelMagErr_z'};
Dictionary(Icat).CatColor     = {'modelMag_u - modelMag_g',...
                                 'modelMag_g - modelMag_r',...
                                 'modelMag_g - modelMag_r',...
                                 'modelMag_r - modelMag_i',...
                                 'modelMag_i - modelMag_z'};
Dictionary(Icat).CatColorErr  = {'sqrt(modelMagErr_u.^2 + modelMagErr_g.^2)',...
                                 'sqrt(modelMagErr_g.^2 + modelMagErr_r.^2)',...
                                 'sqrt(modelMagErr_g.^2 + modelMagErr_r.^2)',...
                                 'sqrt(modelMagErr_r.^2 + modelMagErr_i.^2)',...
                                 'sqrt(modelMagErr_i.^2 + modelMagErr_z.^2)'};      
                             
Icat = Icat + 1;
Dictionary(Icat).Cat          = 'APASS';
Dictionary(Icat).MasterBand   = {'u',            'g',            'r',            'i'};
Dictionary(Icat).CatBand      = {'B',            'g',            'r',            'i'};
Dictionary(Icat).CatBandErr   = {'Berr',         'gerr',         'rerr',         'ierr'};
Dictionary(Icat).CatColor     = {'B - g',...
                                 'g - r',...
                                 'g - r',...
                                 'r - i'};
Dictionary(Icat).CatColorErr  = {'sqrt(Berr.^2 + gerr.^2)',...
                                 'sqrt(gerr.^2 + rerr.^2)',...
                                 'sqrt(gerr.^2 + rerr.^2)',...
                                 'sqrt(rerr.^2 + ierr.^2)'};
                                             
                             
                             
% search
I = find(strcmpi(Cat,{Dictionary.Cat}));
if (isempty(I))
    error('Cat %s not found in dictionary',Cat);
end

Ib = find(strcmpi(Band,Dictionary(I).MasterBand));
CatBand     = Dictionary(I).CatBand{Ib};
CatBandErr  = Dictionary(I).CatBandErr{Ib};
CatColor    = Dictionary(I).CatColor{Ib};
CatColorErr = Dictionary(I).CatColorErr{Ib};




