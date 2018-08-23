function Counter=gaia_dr1_build_cat
% Build the GAIA-DR1 fast access catalog
% Example: Counter=VO.prep.GAIA.gaia_dr1_build_cat

RAD = 180./pi;

Level       = 7;
DecStep     = 10;
CatColCell  = {'RA','ErrRA','Dec','ErrDec','MagG','ErrMagG','ExcessNoise','ExcessNoiseSig'};
NewCatColCell  = {'RA','Dec','ErrRA','ErrDec','MagG','ErrMagG','ExcessNoise','ExcessNoiseSig'};
OrigFiles   = 'GaiaSource*.hdf5';

Ncol        = numel(CatColCell);
CatCol      = cell2struct(num2cell(1:1:Ncol),CatColCell,2);

[HTM,LevList] = celestial.htm.htm_build(Level);
Nlev = numel(LevList(end).ptr);

DecBuffer     = 2; 

Files = dir(OrigFiles);
Nf    = numel(Files);


MainCat = zeros(Nlev,11);  % [FileInd,LevelIndex,HTMindex,MeanRA,MeanDec,RA1,Dec1,RA2,Dec2,RA3,Dec3];

K = 0;
VarNameInd = 0;
% for each declination range
VecDec = (-90:DecStep:+90);
Ndec   = numel(VecDec);
for Idec=1:1:Ndec-1
    Idec
    Dec1 = VecDec(Idec)-DecBuffer;
    Dec2 = VecDec(Idec+1)+DecBuffer;
    
    % load all HDF5 files and select object in Dec range
    RangeCat = zeros(0,Ncol);
    for If=1:1:Nf
        Cat = Util.IO.loadh(Files(If).name,'V',[]);
        Flag = Cat(:,CatCol.Dec)>Dec1 & Cat(:,CatCol.Dec)<Dec2;
        RangeCat = [RangeCat; Cat(Flag,:)];
    end
    
    
    
    % populate all the HTM in range
    
    
    for Ilev=1:1:Nlev
        %Ilev
       
        IndHTM = LevList(end).ptr(Ilev);
        CenterHTM = celestial.coo.cosined(mean(HTM(IndHTM).cosd));
        if (CenterHTM(2)>=VecDec(Idec)./RAD && CenterHTM(2)<VecDec(Idec+1)./RAD),
            % populate HTM
            Flag = celestial.htm.in_polysphere(RangeCat(:,[1 3])./RAD,HTM(IndHTM).cosd);
            
            CatHTM = RangeCat(Flag,:);
            
            CatHTM = CatHTM(:,[1 3 2 4 5 6 7 8]);
            CatHTM = sortrows(CatHTM,2);
            CatHTM(:,1:2) = CatHTM(:,1:2)./RAD;
            AstC = AstCat;
            AstC.Cat     = CatHTM;
            AstC.ColCell = NewCatColCell;
            AstC         = colcell2col(AstC);
            
            
            K = K + 1;
            Counter(K) = size(CatHTM,1);
            
            %VarNameInd  = VarNameInd + 1;
            %FileNameInd = floor(VarNameInd./100).*100;
            FileNameInd = floor(Ilev./100).*100;
            FileName    = sprintf('GAIA_DR1_HTM%06d.hdf5',FileNameInd);
            fprintf('Writing DecStrip %d HTM index %d to file %s\n',Idec,Ilev,FileName);
            Util.IO.saveh(FileName,CatHTM,sprintf('HTM%06d',Ilev));
            
            % MainCat
            % [FileInd,VarInd,MeanRA,MeanDec,RA1,Dec1,RA2,Dec2,RA3,Dec3];
            MainCat(Ilev,:) = [FileNameInd, Ilev, IndHTM, CenterHTM, HTM(IndHTM).coo(1,:),HTM(IndHTM).coo(2,:), HTM(IndHTM).coo(3,:)];
            
        end
    end
end

sum(Counter)    % 664571170

save MainCat.mat MainCat

% RA to range of 0..2pi
load MainCat.mat
MainCat(MainCat(:,4)<0,4)=MainCat(MainCat(:,4)<0,4)+2.*pi;

Main = AstCat;
Main.Cat = single(MainCat(:,[4 5 2]));
Main.ColCell = {'RA','Dec','IndHTM'};
Main = colcell2col(Main);
Main = sortrows(Main,'Dec');
saveh(Main,'GAIA_DR1_HTMindex.hdf5');

save GAIA_DR1_HTMcolcell.mat NewCatColCell


