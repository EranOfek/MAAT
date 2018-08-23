function read_ucac4
%--------------------------------------------------------------------------
% read_ucac4 function                                            Catalogue
% Description: Read UCAC4 catalog in original format, sort it into an HTM
%              catalog in HDF5 files.
%              Takes a few days to run.
%              use: get_ucac4.m to search local copy of the UCAC4 catalog.
% Input  : null
% Output : null
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: read_ucac4
% Reliable: 2
%--------------------------------------------------------------------------

%%
RAD = 180./pi;

I = 0;
I = I + 1;
ColData(I).Name = 'RA';
ColData(I).Fmt  = 4;
ColData(I).Fun  = @(x) x./(1000.*3600.*RAD);   % mas -> rad
ColData(I).Units= 'rad';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'Dec';
ColData(I).Fmt  = 4;
ColData(I).Fun  = @(x) -pi./2 + x./(1000.*3600.*RAD);   % mas -> rad
ColData(I).Units= 'rad';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'MagModel';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./1000;   % millimag -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 20000;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'MagAper';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./1000;   % millimag -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 20000;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrMag';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./100;   % millimag -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 99;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ObjType';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'DSflag';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;   
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrRA';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;   % mas   s.e. at central epoch in RA (*cos Dec)
ColData(I).Units= 'mas';
ColData(I).Null = 255;
ColData(I).Add  = 128;
I = I + 1;
ColData(I).Name = 'ErrDec';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;   % mas
ColData(I).Units= 'mas';
ColData(I).Null = 255;
ColData(I).Add  = 128;
I = I + 1;
ColData(I).Name = 'Nim';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'Nobs';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;   
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'Nep';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x; 
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'EpochMeanRA';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) 1900+x./100;   % -> yr
ColData(I).Units= 'yr';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'EpochMeanDec';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) 1900+x./100;   % -> yr
ColData(I).Units= 'yr';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'PM_RA';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./10;   % -> mas/yr    proper motion in RA*cos(Dec) 
ColData(I).Units= 'mas/yr';
ColData(I).Null = 0;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'PM_Dec';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./10;   % -> mas/yr
ColData(I).Units= 'mas/yr';
ColData(I).Null = 0;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrPM_RA';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./10;   % -> mas/yr    s.e. of pmRA * cos Dec 
ColData(I).Units= 'mas/yr';
ColData(I).Null = 500;
ColData(I).Add  = 128;
I = I + 1;
ColData(I).Name = 'ErrPM_Dec';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./10;   % -> mas/yr    s.e. of pmRA * cos Dec 
ColData(I).Units= 'mas/yr';
ColData(I).Null = 500;
ColData(I).Add  = 128;
I = I + 1;
ColData(I).Name = 'ID2MASS';
ColData(I).Fmt  = 4;
ColData(I).Fun  = @(x) x;
ColData(I).Units= '';
ColData(I).Null = 0;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'MagJ';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./1000;   % millimag -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 20000;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'MagH';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./1000;   % millimag -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 20000;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'MagK';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./1000;   % millimag -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 20000;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'FlagJ';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;   % 2MASS cc_flg*10 + ph_qual flag for J 
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'FlagH';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;   % 2MASS cc_flg*10 + ph_qual flag for H
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'FlagK';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;   % 2MASS cc_flg*10 + ph_qual flag for K
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrMagJ';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./100;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 99;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrMagH';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./100;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 99;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrMagK';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./100;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 99;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'MagB';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./1000;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 20000;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'MagV';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./1000;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 20000;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'Magg';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./1000;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 20000;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'Magr';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./1000;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 20000;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'Magi';
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x./1000;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 20000;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrMagB';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./100;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 99;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrMagV';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./100;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 99;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrMagg';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./100;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 99;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrMagr';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./100;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 99;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ErrMagi';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x./100;   % -> mag
ColData(I).Units= 'mag';
ColData(I).Null = 99;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'FlagYale';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;    %  Yale SPM g-flag*10  c-flag
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'FlagExtCat';   % external catalogs flags
ColData(I).Fmt  = 4;
ColData(I).Fun  = @(x) x;    %  
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'FlagLEDA';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;    %  
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'Flag2MASSext';
ColData(I).Fmt  = 1;
ColData(I).Fun  = @(x) x;    %  
ColData(I).Units= '';
ColData(I).Null = Inf;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ID';
ColData(I).Fmt  = 4;
ColData(I).Fun  = @(x) x;    %  
ColData(I).Units= '';
ColData(I).Null = 0;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'ZoneUCAC2';  
ColData(I).Fmt  = 2;
ColData(I).Fun  = @(x) x;    %  
ColData(I).Units= '';
ColData(I).Null = 0;
ColData(I).Add  = 0;
I = I + 1;
ColData(I).Name = 'RecordUCAC2'; 
ColData(I).Fmt  = 4;
ColData(I).Fun  = @(x) x;    %  
ColData(I).Units= '';
ColData(I).Null = 0;
ColData(I).Add  = 0;


Ncol = numel(ColData);
for Icol=1:1:Ncol,
    switch ColData(Icol).Fmt
        case 1
            ColData(Icol).Fmt = '*int8';
        case 2
            ColData(Icol).Fmt = '*int16';
        case 4
            ColData(Icol).Fmt = '*int32';
        otherwise
            error('Unknwon Fmt option');
    end
end





%%

Base = '/home/astro-sw/shared/cat/UCAC4/u4b/';

Level = 7;

[HTM,LevList]=htm_build(Level);
Leafs = LevList(Level);
HTML = HTM(Leafs.ptr);
Nhtm = numel(HTML);



Step    = 10;   % deg
Margin  = 2;    % deg
DecList = (-90:Step:90).';
Ndec    = numel(DecList) - 1;
ZoneList= (-90:0.2:90-0.2).';

Counter = 0;
OutCounter = 1;
for Idec=1:1:Ndec,
    DecL = DecList(Idec)   - Margin;
    DecH = DecList(Idec+1) + Margin;
    
    Izone = find(ZoneList>DecL & ZoneList<DecH);
    Nz    = numel(Izone);
    
    UC = AstCat;
    UC.ColCell = cell(1,Ncol);
    for Icol=1:1:Ncol,
        UC.ColCell{Icol} = ColData(Icol).Name;
    end
    UC = colcell2col(UC);

    for Iz=1:1:Nz,
        ZoneInd  = Izone(Iz);
        ZoneName = sprintf('%sz%03d',Base,ZoneInd);


        FID = fopen(ZoneName,'r');

        Isrc = 0;
        Cont = true;

        while (Cont),
            Isrc = Isrc + 1;

            for Icol=1:1:Ncol,
                Read = fread(FID,1,ColData(Icol).Fmt,'l') + ColData(Icol).Add;
                Cont = ~feof(FID);
                if (Cont)
                    if (Read==ColData(Icol).Null),
                        UC(Iz).Cat(Isrc,Icol) = NaN;
                    else
                        UC(Iz).Cat(Isrc,Icol) = ColData(Icol).Fun(double(Read));
                    end
                end
            end
        end
        fclose(FID);

        UC(Iz).Cat = UC(Iz).Cat(1:Isrc-1,:);
    end
    
    UCm = merge(UC);
    UCm = sortrows(UCm,2);

    
    for Ihtm=1:1:Nhtm,
        %HTML(Ihtm).coo
        MinDec = min(HTML(Ihtm).coo(:,2));
        if (MinDec>=DecList(Idec)./RAD && MinDec<DecList(Idec+1)./RAD),
            % populate HTM index Ihtm
            Ihtm
            kkk
            Counter = Counter + 1;
            
            Flag  = in_polysphere(UCm.Cat(:,1:2),HTML(Ihtm).coo);
            UChtm(Counter) = row_select(UCm,Flag);
            
            FileName = sprintf('UCAC4htm_%04d.hdf5',OutCounter);
            VarName  = sprintf('UCAC4_htm%07d',Ihtm);
            saveh(FileName,UChtm(Counter).Cat,VarName);
            
            if (Counter==100),
                % save 100 files
                OutCounter = OutCounter + 1;
                                
                Counter = 0;
                UChtm = AstCat;
            end
        end
    end          
end


% Ptr
Files = dir('UCAC4htm_*.hdf5');
Nf = numel(Files);
Size = 0;
C = 0;
UCAC4htm = AstCat;
UCAC4htm.Cat = zeros(Nhtm,4);
for If=1:1:Nf,
    If
    Tmp = loadh(Files(If).name);
    Fields = fieldnames(Tmp);
    Nfield = numel(Fields);
    for Ifield=1:1:Nfield,
        C = C + 1;
        Ihtm = str2double(Fields{Ifield}(10:end));
        
        %[MeanRA,MeanDec] = tri_equidist_center(HTML(Ihtm).coo(:,1).', HTML(Ihtm).coo(:,2).');
        MeanRA  = mean(HTML(Ihtm).coo(:,1).');
        MeanDec = mean(HTML(Ihtm).coo(:,2).');
        
        UCAC4htm.Cat(C,:) = [MeanRA, MeanDec, If, Ihtm];
        
    end
end


 UCAC4htm.Cat = single(UCAC4htm.Cat);
 UCAC4htm.ColCell = {'RA','Dec','FileInd','Ptr'};
 UCAC4htm = colcell2col(UCAC4htm);
 UCAC4htm.ColUnits = {'rad','rad','',''};
 UCAC4htm = sortrows(UCAC4htm,2);
 UCAC4htm.Reference = 'UCAC4 catalog of HTM grid';
 
 save UCAC4htm.mat UCAC4htm
 tic;load2('UCAC4htm.mat');toc




%%


% Files  = dir('z*');
% Nfiles = numel(Files);
% 
% %Zone = 'z001';
% 
% 
% 
% for Ifile=1:1:100,
%     %Nfiles,
%     UC(Ifile).Cat = zeros(2e5,Ncol);
%     if (Ifile>1),
%         UC(Ifile).ColCell = UC(Ifile-1).ColCell;
%         UC(Ifile).Col     = UC(Ifile-1).Col;
%     end
%     
%     FID = fopen(Files(Ifile).name,'r');
%     
%     Isrc = 0;
%     Cont = true;
% 
% 
%     while (Cont),
%         Isrc = Isrc + 1;
% 
%         for Icol=1:1:Ncol,
%    
%             Read = fread(FID,1,ColData(Icol).Fmt,'l') + ColData(Icol).Add;
% 
% 
%             Cont = ~feof(FID);
%             if (Cont)
%                 if (Read==ColData(Icol).Null),
%                     UC(Ifile).Cat(Isrc,Icol) = NaN;
%                 else
%                     UC(Ifile).Cat(Isrc,Icol) = ColData(Icol).Fun(double(Read));
%                 end
%             end
%         end
% 
%     end
%     fclose(FID);
% 
%     UC(Ifile).Cat = UC(Ifile).Cat(1:Isrc-1,:);
% end
% 
% UCm = merge(UC);
% UCm = sortrows(UCm,2);
% 
% Flag=in_polysphere(AstC,Cols,Corners,Crit);

%% check validity
Files = dir('UCAC4htm_*.hdf5');
Nf = numel(Files);
Size = 0;
for If=1:1:Nf,
    If
    Tmp = loadh(Files(If).name);
    Fields = fieldnames(Tmp);
    Nfield = numel(Fields);
    for Ifield=1:1:Nfield,
        Size = Size + size(Tmp.(Fields{Ifield}),1);
    end
end

%% construct master list


