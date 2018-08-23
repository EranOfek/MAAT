function prep_wise_htm_cat
% reformat the IRSA/WISE catalog files into an HDF5/HTM catalog

%%
CatBaseName = 'WISE';
RAD = 180./pi;

NfileInHTM = 100;

Level = 8;
[HTM,LevelH] = celestial.htm.htm_build(Level);
LevelL = LevelH(Level);
Nhtm   = numel(LevelL.ptr);

FID = fopen('colnames');
C   = textscan(FID,'%s\n');
ColCell = C{1};
fclose(FID);


% RA, Dec, ErrRA, ErrDec, Mag_W1, MagErr_W1, SNR_W1, RChi2_W1, Mag_W2, MagErr_W2, SNR_W2, RChi2_W2, Mag_W3, MagErr_W3, SNR_W3, RChi2_W3, Mag_W4, MagErr_W4, SNR_W4, RChi2_W4,
%  Sky_W1, Sky_W2, Sky_W3, Sky_W4,
%   VarL_W1, MeanMJD_W1, VarL_W2, VarL_W3, VarL_W4
%   Dist_2MASS, PA_2MASS
%   Mag_J, MagErr_J, Mag_H, MagErr_H, Mag_K, MagErr_K

%  http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec2_2a.html
Format = '%*s %f %f %f %f %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f %*s %*s %*s %*s %f %*s %*s %*s %*s %f %*s %*s %*s %*s %f %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f %*s %*s %f %*s %*s %*s %*s %*s %f %*s %*s %*s %*s %*s %*s %*s %*s %f %*s %*s %*s %*s %*s %*s %*s %*s %f %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f %f %*s %f %f %f %f %f %f %*s %*s %*s %*s %*s %*s \n';

Files = dir('wise-allsky-cat-part*');
Nf    = numel(Files);
StatusHTM = false(Nhtm,1);
Nsrc      = nan(Nhtm,1);
%%
If1 = 1;
%system(sprintf('bzip2 -d %s',Files(If1).name));

FID = fopen(Files(If1).name);
Line = fgetl(FID);
fclose(FID);
Ncol = numel(strfind(Line,'|'));

%
% read first file
FID = fopen(Files(If1).name);
C = textscan(FID,Format,'Delimiter','|');
fclose(FID);
%system(sprintf('bzip2 %s',Files(If1).name(1:end-4)));

Mat1 = cell2mat(C);
Mat1(:,1:2) = Mat1(:,1:2)./RAD;
    
for If=(If1+1):1:Nf
    If
    %system(sprintf('bzip2 -d %s',Files(If).name));
    FID = fopen(Files(If).name);
    C = textscan(FID,Format,'Delimiter','|');
    fclose(FID);
    %system(sprintf('bzip2 %s',Files(If).name(1:end-4)));
    
    Mat2 = cell2mat(C);
    Mat2(:,1:2) = Mat2(:,1:2)./RAD;
    
    Mat  = [Mat1; Mat2];
    Mat  = sortrows(Mat,2);
    Mat1 = Mat2;
    clear Mat2;
    
    MaxDec = max(Mat(:,2));
    
    for Iptr=1:1:Nhtm
        Ihtm = LevelL.ptr(Iptr);
        htmMaxDec = max(HTM(Ihtm).coo(:,2));
        if (htmMaxDec<MaxDec && ~StatusHTM(Iptr))
            
            tic;
            MeanRA  = mean(HTM(Ihtm).coo(:,1));
            MeanDec = mean(HTM(Ihtm).coo(:,2));
            D = celestial.coo.sphere_dist_fast(Mat(:,1),Mat(:,2),MeanRA,MeanDec);
            FlagD = D<(2.5./RAD);
            MatC  = Mat(FlagD,:);
            Flag = celestial.htm.in_polysphere(MatC(:,1:2),HTM(Ihtm).coo,2);
            Nsrc(Ihtm) = sum(Flag);
            
            if (Nsrc(Ihtm)>0)
                [FileName,DataName]=HDF5.get_file_var_from_htmid(CatBaseName,Ihtm,NfileInHTM);


                HDF5.save_cat(FileName,DataName,MatC(Flag,:),2,30);
                StatusHTM(Iptr) = true;
                save StatusHTM.mat StatusHTM Nsrc
                [If, Ihtm, Iptr, toc]
            end
        end
    end
    
end
