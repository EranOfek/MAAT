function  ColCell=prep_ukidss_htm

% ColCell=VO.prep.prep_ukidss_htm

% DR info: http://wsa.roe.ac.uk/theSurveys.html
% SQL: http://wsa.roe.ac.uk:8080/wsa/SQL_form.jsp
% SELECT ra, dec, sigRA, sigDec, Epoch, pStar, pGalaxy, pNoise, 
% yPsfMag, yPsfMagErr, ySerMag2D, ySerMag2DErr, yEll, yPA,
% j_1PsfMag, j_1PsfMagErr, j_1SerMag2D, j_1SerMag2DErr, j_1Ell, j_1PA,
% j_2PsfMag, j_2PsfMagErr, j_2SerMag2D, j_2SerMag2DErr, j_2Ell, j_2PA,
% hPsfMag, hPsfMagErr, hSerMag2D, hSerMag2DErr, hEll, hPA,
% kPsfMag, kPsfMagErr, kSerMag2D, kSerMag2DErr, kEll, kPA
% FROM lasSource
% WHERE  dec>-20.0 and dec<=-10.0

RAD = 180./pi;

File = dir('ukidss_*.fits');
Sp   = regexp({File.name},'_','split');
Nfile = numel(File);
for Ifile=1:1:Nfile
    Dec1 = Sp{Ifile}{2};
    Dec2 = Sp{Ifile}{3}(1:end-5);
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
     I=1;
    T(I) = FITS.read_table(File(Ifile).name);
    if (Ifile>1)
        T(I+1) = FITS.read_table(File(Ifile-1).name);
    end
    if (Ifile<Nfile)
        T(I+1) = FITS.read_table(File(Ifile+1).name);
    end
    T = merge(T);
    T = sortrows(T,2);
    T.Cat(T.Cat<-1e8) = NaN;
    
    VO.prep.build_htm_catalog(T.Cat,...
                'CatName','UKIDSS',...
                'HTM_Level',8,...
                'SaveInd',false,...
                'DecRange',[D1(Ifile), D2(Ifile)]./RAD);
            toc
end
ColCell = T.ColCell;