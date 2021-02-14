function prep_gaiadre3_htm(varargin)
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



DefV.URL                  = 'http://cdn.gea.esac.esa.int/Gaia/gedr3/gaia_source/';
%http://o501sd3ggsiwxpbv8jpv.push-19.cdn77.com/Gaia/gdr2/gaia_source/csv/';
DefV.Dir                  = 'zz1';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

%%
% get all files for GAIA-DR2 sources
List=www.find_urls(InPar.URL,'strfind','.csv.gz');
www.pwget(List,'',15);

% !gzip -d *.gz



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
% list of required columns and their new names:
%            GAIA name , catsHTM name, catsHTM position
ColNames = {'ref_epoch','Epoch',3;...
            'ra',       'RA',   1;...
            'ra_error', 'ErrRA',4;...
            'dec',      'Dec',  2;...
            'dec_error','ErrDec',5;...
            'parallax','Plx',    6;...
            'parallax_error','ErrPlx',7;...
            'pmra','PMRA',            8;...
            'pmra_error','ErrPMRA',   9;...
            'pmdec','PMDec',         10;...
            'pmdec_error','ErrPMDec',11;...
            'ra_dec_corr','RA_Dec_Corr',12;...
            'astrometric_n_good_obs_al','NobsAst',13;...
            'astrometric_excess_noise','ExcessNoise',14;...
            'astrometric_excess_noise_sig','ExcessNoiseSig',15;...
            'astrometric_chi2_al','Chi2Ast',16;...
            'nu_eff_used_in_astrometry','DofAst',17;...
            'phot_g_n_obs','NGphot',18;...
            'phot_g_mean_flux','Mag_G',19;...
            'phot_g_mean_flux_error','ErrMag_G',20;...
            'phot_bp_mean_flux','Mag_BP',21;...
            'phot_bp_mean_flux_error','ErrMag_BP',22;...
            'phot_rp_mean_flux','Mag_RP',23;...
            'phot_rp_mean_flux_error','ErrMag_RP',24;...
            'phot_bp_rp_excess_factor','BPRP_Excess',25;...
            'dr2_radial_velocity','RV',26;...
            'dr2_radial_velocity_error','ErrRV',27;...
            'dr2_rv_template_teff','Teff',28;...
            'dr2_rv_template_logg','LogG',29;...
            'dr2_rv_template_fe_h','FeH',30};


FID = fopen('GaiaSource_779915-780373.csv','r');
Line = fgetl(FID);
fclose(FID);
SpLine = regexp(Line,',','split');


ColCell = ColNames(:,2).';
N = size(ColNames,1);
for I=1:1:N
    ColInGAIA(I) = find(strcmp(ColNames{I},SpLine));
   
end

Format = '';
for Icol=1:1:numel(SpLine)
     if ~any(Icol==ColInGAIA)
         % ignore
         Format = sprintf('%s%%*s ',Format);
     else
         % read
         Format = sprintf('%s%%f ',Format);
     end
end

Format = [Format, '%*[^\n]'];
    
   





%%
[~,SI]=sort(cell2mat(ColNames(:,3)));

Dir = dir('GaiaSource_*.csv');

Nfile = numel(Dir);
SumC  = nan(Nfile,7);
K = 0;

for If=1:1:Nfile
    If
    tic;
    File = Dir(If).name;
    FID = fopen(File,'r');
    C=textscan(FID,Format,'HeaderLines',1,'Delimiter',',');
    fclose(FID);

    %delete(File);
    Mat = cell2mat(C);
    Mat = Mat(:,SI);
    
    % flux to AB mag
    % https://gea.esac.esa.int/archive/documentation/GEDR3/Data_processing/chap_cu5pho/cu5pho_sec_photProc/cu5pho_ssec_photCal.html#Ch5.T2
    
    ZP   = 25.8010;  % AB
    ColM = find(strcmp('Mag_G',ColCell));
    ColE = find(strcmp('ErrMag_G',ColCell));
    
    Mag = -2.5.*log10(Mat(:,ColM)) + ZP;
    Err = 1.086.* (Mat(:,ColE)./Mat(:,ColM));
    Mat(:,ColM) = Mag;
    Mat(:,ColE) = Err;
    
    ZP   = 25.3540;  % AB
    ColM = find(strcmp('Mag_BP',ColCell));
    ColE = find(strcmp('ErrMag_BP',ColCell));
    
    Mag = -2.5.*log10(Mat(:,ColM)) + ZP;
    Err = 1.086.* (Mat(:,ColE)./Mat(:,ColM));
    Mat(:,ColM) = Mag;
    Mat(:,ColE) = Err;
    
    ZP   = 25.1040;  % AB
    ColM = find(strcmp('Mag_RP',ColCell));
    ColE = find(strcmp('ErrMag_RP',ColCell));
    
    Mag = -2.5.*log10(Mat(:,ColM)) + ZP;
    Err = 1.086.* (Mat(:,ColE)./Mat(:,ColM));
    Mat(:,ColM) = Mag;
    Mat(:,ColE) = Err;
    
    
    K = K + 1;
    SaveFile = sprintf('H5_%d.hdf5',K);
    cd(InPar.Dir)
    HDF5.save(Mat,SaveFile);
    cd ..

    Ifile = 1;
    SumC(K,:) = [Ifile, K, min(Mat(:,1)), max(Mat(:,1)), min(Mat(:,2)), max(Mat(:,2)), size(Mat,1)];
    toc
end
    
save SumC.mat SumC


%%
%cd /raid/eran/catsHTM/GAIA/DRE3/
cd /data/euler/catsHTM/GAIA/DRE3/
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
%     switch StrInd(1)
%         case '1'
%             if (numel(StrInd)==1)
%                 Dir = 'zz10';
%             else
%                 Dir = sprintf('zz%s',StrInd(1:2));
%             end
%         otherwise
%             Dir = sprintf('zz%s',StrInd(1));
%     end
    %[StrInd, Dir]
    Dir = 'zz1';
    
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
%         switch StrInd(1)
%             case '1'
%                 if (numel(StrInd)==1)
%                     Dir = 'zz10';
%                 else
%                     Dir = sprintf('zz%s',StrInd(1:2));
%                 end
%             otherwise
%                 Dir = sprintf('zz%s',StrInd(1));
%         end
        Dir = 'zz1';
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
    
