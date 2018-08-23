function prep_2mass_htm
% prepare 2MASS catalog in HDF5/HTM format


% http://irsa.ipac.caltech.edu/2MASS/download/allsky/
% http://irsa.ipac.caltech.edu/2MASS/download/allsky/format_psc.html

RAD = 180./pi;

Level = 8;
[HTM,LevelH] = celestial.htm.htm_build(Level);
LevelL = LevelH(end);
Nhtm   = numel(LevelL.ptr);

ColCell = {'RA', 'Dec', 'ErrMaj', 'Mag_J', 'MagErr_J', 'Mag_H', 'MagErr_H', 'Mag_K', 'MagErr_K', 'Epoch'};
Ncol = numel(ColCell);

Format = '%f %f %f %*s %*s %*s %f %f %*s %*s  %f %f %*s %*s  %f %f %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f %*[^\n]';
Files  = dir('psc_*');
Nf     = numel(Files);
StepN  = 5;
DecLimit1 = -pi./2;
for If=1:StepN:Nf
    
    Mat = zeros(0,Ncol);
    MinI = max(If-2,1);
    MaxI = min(If+StepN+2,Nf);
    IfL = (MinI:1:MaxI);
    NL  = numel(IfL);
    for IL=(If-2):1:(If+StepN+2)
        if (IL>=1 && IL<=Nf)
            FID = fopen(Files(IL).name,'r');
            C   = textscan(FID,Format,'Delimiter','|','TreatAsEmpty','\N');
            fclose(FID);
            Mat = [Mat; cell2mat(C)];
            if (IL==(If+StepN))
                DecLimit2 = max(C{2})./RAD;
            end
        end
    end
    Mat(:,1:2) = Mat(:,1:2)./RAD;
    Mat = sortrows(Mat,2);
    
    if (If>(Nf-StepN))
        DecLimit2 = pi./2;
    end
    
    [If, DecLimit1.*RAD, DecLimit2.*RAD, min(Mat(:,2)).*RAD, max(Mat(:,2)).*RAD]
    
    VO.prep.build_htm_catalog(Mat,'CatName','TMASS','HTM_Level',Level,'DecRange',[DecLimit1 DecLimit2],'SaveInd',false);
    DecLimit1 = DecLimit2;
end
