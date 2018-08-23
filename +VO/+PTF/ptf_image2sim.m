function Sim=ptf_image2sim(Files,Cat,varargin)
% Load PTF fits images and catalogs into a SIM object
% Package: VO.PTF
% Description: Load PTF fits images and catalogs into a SIM object.
%              The catalog names are constructed from the file names.
%              Either sextractor or daophot catalogs can be loaded.
% Input  : - List of fits file to load. See Util.files.create_list for
%            options.
%          - Catalog type. If empty, do not load catalog.
%            If 'daophot' then attempt to load daophot catalog.
%            If sex' then attempt to load sextractor catalog.
%            Default is empty.
% Output : - A SIM object with the images and catalogs.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=VO.PTF.ptf_image2sim('PTF*.fits','daophot')
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<2)
    Cat = [];
end


[~,ListFile] = Util.files.create_list(Files,NaN);


Sim = FITS.read2sim(ListFile);

if (~isempty(Cat))
    if (ischar(Cat))
        
        switch lower(Cat)
            case 'daophot'
                ListCat = strrep(ListFile,'_i_p_scie_t','_c_d_scie_t');
                ListCat = strrep(ListCat,'.fits','.ctlg');
            case 'sex'
                ListCat = strrep(ListFile,'_i_p_scie_t','_c_p_scie_t');
                ListCat = strrep(ListCat,'.fits','.ctlg');
            otherwise
                error('Unknown catalog option');
        end
    end
    [Out]=FITS.read_table(ListCat);
    Sim = astcat2sim(Out,Sim);
    
end
            