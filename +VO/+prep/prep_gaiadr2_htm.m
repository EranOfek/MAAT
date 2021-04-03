function prep_gaiadr2_htm(varargin)
% SHORT DESCRIPTION HERE
% Package: VO
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.prep.prep_gaia_htm
% Reliable: 
%--------------------------------------------------------------------------



DefV.URL                  = 'http://o501sd3ggsiwxpbv8jpv.push-19.cdn77.com/Gaia/gdr2/gaia_source/csv/';
DefV.Dir                  = 'zz1';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

%%
% get all files for GAIA-DR2 sources
List=www.find_urls(InPar.URL,'strfind','.csv.gz');
%www.pwget(List,'',15);


%% read single file
% File = 'GaiaSource_1000172165251650944_1000424567594791808.csv';
% FID  = fopen(File,'r');
% Line = fgetl(FID);
% ColCellAll = regexp(Line,',','split');
% fclose(FID);

%%
ColCell = {'RA','Dec','Epoch','ErrRA','ErrDec','Plx','ErrPlx','PMRA','ErrPMRA','PMDec','ErrPMDec','RA_Dec_Corr',...
           'ExcessNoise','ExcessNoiseSig','MagErr_G','Mag_G','MagErr_BP','Mag_BP','MagErr_RP','Mag_RP',...
           'RV','ErrRV','VarFlag','Teff','Teff_low','Teff_high','A_G'};
    
%%
% ref_epoch  5
% ra 6
% ra_err 7
% dec 8
% dec_err 9
% parallax 10
% parallax_err 11
% pmra 13
% pmra 14
% pmdec 15
% pmdec_err 16
% ra_dec_corr 17
% astrom_excess_noise 33
% astrom_excess_noise_sig 34
% phot_g_mean_flux_over_error 50
% phot_g_mean_mag 51
% phot_bp_mean_flux_over_error 55
% phot_bp_mean_mag 56
% phot_rp_mean_flux_over_error 60
% phot_rp_mean_mag 61
% radial_velocity 67
% radial_celocity_err 68
% phot_variable_flag 73
% teff 79
% teff_low 80
% teff_high 81
% A_G 82






%%
%   %                                                       17
Format = '%*s %*s %*s %*s %f %f %f %f %f %f  %f %*s %f %f %f %f %f %*s %*s %*s  %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s  %*s %*s %f %f %*s %*s %*s %*s %*s %*s  %*s %*s %*s %*s %*s %*s %*s %*s %*s %f  %f %*s %*s %*s %f %f %*s %*s %*s %f  %f %*s %*s %*s %*s %*s %f %f %*s %*s  %*s %*s %s %*s %*s %*s %*s %*s %f %f  %f %f %*[^\n]';

%Dir = dir('GaiaSource_*.csv');

Nfile = numel(List); %numel(Dir);
SumC  = nan(Nfile,7);
K = 0;
for Ifile=1:10:Nfile
    Ifile
    %Nfile
    tic;
    Ifile2 = min(Ifile+9,Nfile);
    
    
    www.pwget(List(Ifile:1:Ifile2),'',10);
    system('gzip -d *.gz');
    Dir = dir('GaiaSource*.csv');
    Nf  = numel(Dir);

    
    Ifile
    for If=1:1:Nf
        File = Dir(If).name;
        FID = fopen(File,'r');
        C=textscan(FID,Format,'HeaderLines',1,'Delimiter',',');
        fclose(FID);

        delete(File);

        C{23} = regexprep(C{23},'NOT_AVAILABLE','-1');
        C{23} = regexprep(C{23},'VARIABLE','1');
        C{23} = regexprep(C{23},'CONSTANT','0');
        C{23}=cellfun(@str2double,C{23});

        %

        Mat = [C{:}];
        % S/N to mag err
        Mat = Mat(:,[2 4 1 3 5, 6:end]);
        Mat(:,[15 17 19]) = 1.086./Mat(:,[15 17 19]);

        K = K + 1;
        SaveFile = sprintf('H5_%d.hdf5',K);
        cd(InPar.Dir)
        HDF5.save(Mat,SaveFile);
        cd ..

        
        SumC(K,:) = [Ifile, K, min(Mat(:,1)), max(Mat(:,1)), min(Mat(:,2)), max(Mat(:,2)), size(Mat,1)];
    end
    
    save SumC.mat SumC
    toc
end

%%
cd /raid/eran/catsHTM/GAIA/DR2/
load SumC.mat

VecDec = (-90:1:90)';
Ndec   = numel(VecDec);

Nsum = size(SumC,1);
Status = false(Nsum,Ndec-1);
for Isum=1:1:Nsum
    %
    tic;
    CurInd = SumC(Isum,2);
        
    StrInd = num2str(CurInd);
    switch StrInd(1)
        case '1'
            if (numel(StrInd)==1)
                Dir = 'zz10';
            else
                Dir = sprintf('zz%s',StrInd(1:2));
            end
        otherwise
            Dir = sprintf('zz%s',StrInd(1));
    end
    %[StrInd, Dir]

    File = sprintf('H5_%s.hdf5',StrInd);
    cd(Dir)
    Data = h5read(File,'/V');
    cd ..
    
    % 
    for Idec=1:1:Ndec-1
        D1 = VecDec(Idec);
        D2 = VecDec(Idec+1);
        
        Status(Isum,Idec) = any(Data(:,2)>=D1 & Data(:,2)<D2);
    end
    [Isum, toc]
end

save Status.mat Status
%%
RAD = 180./pi;
        
load Status.mat
load SumC.mat

VecDec = (-90:1:90)';
Ndec   = numel(VecDec);
SumDec = zeros(Ndec-1,1);
for Idec=1:1:Ndec-1
   
    tic;
    
    D1 = VecDec(Idec);
    D2 = VecDec(Idec+1);
    
    Flag    = Status(:,Idec);
    IndRead = find(Flag);
    Nread   = numel(IndRead);
    
    [Idec, D1, D2, Nread]
    %Flag = (MinDec>=D1 & MinDec<=D2) | (MaxDec>=D1 & MaxDec<=D2) | (MinDec<=D1 & MaxDec>=D2);
    
    Ind = SumC(Flag,2);
    Nind = numel(Ind);
    DataDec = zeros(0,27);
    for Iread=1:1:Nread
        Ind = IndRead(Iread);
        
        StrInd = num2str(Ind);
        switch StrInd(1)
            case '1'
                if (numel(StrInd)==1)
                    Dir = 'zz10';
                else
                    Dir = sprintf('zz%s',StrInd(1:2));
                end
            otherwise
                Dir = sprintf('zz%s',StrInd(1));
        end
        %[StrInd, Dir]
        
        File = sprintf('H5_%s.hdf5',StrInd);
        cd(Dir)
        Data = h5read(File,'/V');
        cd ..
        
        FlagSrc = Data(:,2)>=D1 & Data(:,2)<D2;
        sum(FlagSrc);
        Data = Data(FlagSrc,:);
        
        All(Iread).D = Data.';
        %DataDec = [DataDec; Data];
    end
    DataDec = [All.D].';
    clear All;
    DataDec = sortrows(DataDec,2);
    SumDec(Idec) = size(DataDec,1);
    
    if (D1<0)
        D1s = 'm';
    else
        D1s = '';
    end
    if (D2<0)
        D2s = 'm';
    else
        D2s = '';
    end
    FileW = sprintf('GaiaDR2_%s%d_%s%d.hdf5',D1s,abs(D1),D2s,abs(D2));
    DataDec(:,1:2) = DataDec(:,1:2)./RAD;
    HDF5.save(DataDec,FileW,'/V');
    
    toc
end
    
