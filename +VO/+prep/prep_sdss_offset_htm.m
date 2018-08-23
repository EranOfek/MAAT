function prep_sdss_offset_htm


% ColCell=VO.prep.prep_sdss_offset_htm


% SELECT ra, dec, type, flags,
%   modelMag_u, modelMag_g, modelMag_r, modelMag_i, modelMag_z,
%   modelMagErr_u, modelMagErr_g, modelMagErr_r, modelMagErr_i, modelMagErr_z, 
%   offsetRa_u, offsetRa_g, offsetRa_r, offsetRa_i, offsetRa_z,
%   offsetDec_u, offsetDec_g, offsetDec_r, offsetDec_i, offsetDec_z,
%   psffwhm_u, psffwhm_g, psffwhm_r, psffwhm_i, psffwhm_z,
%   TAI_r into mydb.MyTable_1 from PhotoPrimary
%   WHERE modelMag_z<20.0

ColCell = {'RA','Dec','type','flags','modelMag_u','modelMag_g','modelMag_r','modelMag_i','modelMag_z','modelMagErr_u','modelMagErr_g','modelMagErr_r','modelMagErr_i','modelMagErr_z','offsetRa_u','offsetRa_g','offsetRa_r','offsetRa_i','offsetRa_z','offsetDec_u','offsetDec_g','offsetDec_r','offsetDec_i','offsetDec_z','psffwhm_u','psffwhm_g','psffwhm_r','psffwhm_i','psffwhm_z','TAI_r'};

RAD = 180./pi;

File = dir('S_*.mat');
Sp   = regexp({File.name},'_','split');
Nfile = numel(File);
for Ifile=1:1:Nfile
    Dec1 = Sp{Ifile}{2};
    Dec2 = Sp{Ifile}{3}(1:end-4);
    if (strcmp(Dec1(1),'m'))
        D1(Ifile) = -str2double(Dec1(2:end));
    else
        D1(Ifile) = str2double(Dec1);
    end
    if (strcmp(Dec2(1),'m'))
        D2(Ifile) = -str2double(Dec2(2:end));
    else
        D2(Ifile) = str2double(Dec2);
    end
end
[~,SI] = sort(D1);

% sort files by declination
File = File(SI);
D1 = D1(SI);
D2 = D2(SI);
    


for Ifile=1:1:Nfile
     tic;
     clear T
     T = AstCat;
     I=1;
     T(I).Cat = Util.IO.load2(File(Ifile).name);
     T(I).ColCell = ColCell;
     T = colcell2col(T);
     
    if (Ifile>1)
        I = I + 1;
        T(I).Cat = Util.IO.load2(File(Ifile).name);
        T(I).ColCell = ColCell;
        T = colcell2col(T);
    end
    if (Ifile<Nfile)
        I = I + 1;
        T(I).Cat = Util.IO.load2(File(Ifile).name);
        T(I).ColCell = ColCell;
        T = colcell2col(T);
    end
    T = merge(T);
    T.Cat(:,1:2) = T.Cat(:,1:2)./RAD;
    T = sortrows(T,2);
    %T.Cat(T.Cat<-1e8) = NaN;
    T.Cat(:,end) = T.Cat(:,end)./86400 + 2400000.5;
    
    
    VO.prep.build_htm_catalog(T.Cat,...
                'CatName','SDSSoffset',...
                'HTM_Level',8,...
                'SaveInd',false,...
                'DecRange',[D1(Ifile), D2(Ifile)]./RAD);
            toc
end
