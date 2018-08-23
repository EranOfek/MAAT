function prep_DECaLS_htm(varargin)
% SHORT DESCRIPTION HERE
% Package: VO.prep
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.prep.prep_DECaLS_htm
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

Fun360 = @(X) (X./360 - floor(X./360)).*360;
CatField = AstCat.CatField;

DefV.BaseURL              = 'ftp://archive.noao.edu/public/hlsp/ls/dr5/sweep/5.0/';
DefV.GetFiles             = false;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

[Link,FN] = www.find_urls_ftp(InPar.BaseURL);

B = regexp(FN,'sweep-(?<MinRA>\d\d\d)(?<SignMinDec>[pm])(?<MinDec>\d\d\d)-(?<MaxRA>\d\d\d)(?<SignMaxDec>[pm])(?<MaxDec>\d\d\d).fits','names');
B = cell2mat(B);
SignMinDec = 1-2.*strcmp({B.SignMinDec},'m');
SignMaxDec = 1-2.*strcmp({B.SignMaxDec},'m');
MinRA  = str2double({B.MinRA});
MinDec = SignMinDec.*str2double({B.MinDec});
MaxRA  = str2double({B.MaxRA});
MaxDec = SignMaxDec.*str2double({B.MaxDec});


if (InPar.GetFiles)
    www.pwget(Link,'',15);
end

TypeDic = {'PSF','REX','EXP','DEV','COMP'};
%ColCell = {'RA','Dec','ErrRA', 'ErrDec','Type','DeltaChi2','Mag_u','Mag_g','Mag_r','Mag_i','Mag_z','Mag_y','Mag_W1','Mag_W2','Mag_W3','Mag_W4',...
%                              'MagErr_u','MagErr_g','MagErr_r','MagErr_i','MagErr_z','MagErr_y','MagErr_W1','MagErr_W2','MagErr_W3','MagErr_W4'};
ColCell = {'RA','Dec','ErrRA', 'ErrDec','Type','DeltaChi2','Flux_u','Flux_g','Flux_r','Flux_i','Flux_z','Flux_y','Flux_W1','Flux_W2','Flux_W3','Flux_W4',...
                              'FluxErr_u','FluxErr_g','FluxErr_r','FluxErr_i','FluxErr_z','FluxErr_y','FluxErr_W1','FluxErr_W2','FluxErr_W3','FluxErr_W4'};
ColUnits = {'rad','rad','arcsec','arcsec','','','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies',...
                                                'nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies','nanomaggies'}
                                            

%%                          
% reformat the catalog into HDF5 (faster)
for If=1:1:Nf
    tic;
    %Nf
   
    try
        T1 = FITS.read_table(FN{If});
    catch
        % try to retrieve file
        warning('Retrieve file')
        FN{If}
        delete(FN{If});
        www.pwget(sprintf('%s%s',InPar.BaseURL,FN{If}));
        pause(1);
        T1 = FITS.read_table(FN{If});
    end
    DelChi2 = T1.(CatField)(:,20:24);
    DelChi2(DelChi2(:,1:end)==0) = Inf;
    DelChi2 = DelChi2(:,1) - min(DelChi2(:,2:end),[],2);

    TypeStr = cellstr(char(T1.(CatField)(:,12:15)));
    Ns   = numel(DelChi2);
    Type = zeros(Ns,1);
    for It=1:1:numel(TypeDic)
        F=strcmp(TypeStr,TypeDic{It});
        Type(F) = F(F).*It;
    end

    Mag =  T1.(CatField)(:,26:35); % nanomagies
    MagErr = T1.(CatField)(:,36:45);
    Mag(Mag==0) = NaN;
    MagErr(MagErr==0) = NaN;
%    
    Data1 = [T1.(CatField)(:,16)./RAD, T1.(CatField)(:,17)./RAD, ...
             3600./sqrt(T1.(CatField)(:,18)), 3600./sqrt(T1.(CatField)(:,19)), ...
             Type, ...
             DelChi2, Mag, MagErr];
             
    HDF5.save(Data1,sprintf('%s.h5',FN{If}));
    toc
    
end


%%
% save the HDF/HTM files

Level = 9;
[HTM,LevelHTM] = celestial.htm.htm_build(Level);
StepRA  = 10;
StepDec = 5;

%%
Nf = numel(FN);
for If=1:1:Nf
    tic;
    Inear = find((MinDec(If)==MinDec | (MinDec(If)-StepDec)==MinDec | (MinDec(If)+StepDec)==MinDec) & ...
                 ((MinRA(If)==MinRA) | Fun360(MinRA(If)-StepRA)==MinRA | Fun360(MinRA(If)+StepRA)==MinRA));

    Nnear = numel(Inear);
    
    for I=1:1:Nnear
        I
        File = sprintf('%s.h5',FN{Inear(I)});
        Data1 = h5read(File,'/V');
        
        clear T1;
        
        if (I==1)
            Data = Data1;
        else
            Data = [Data; Data1];
        end
    end
    Nsrc = size(Data,1);

    VO.prep.build_htm_catalog(Data,'CatName','DECaLS','HTM_Level',9,...
                        'HTM',HTM','LevelHTM',LevelHTM,...
                        'DecRange',[MinDec(If), MaxDec(If)]./RAD,...
                        'RARange',[MinRA(If), MaxRA(If)]./RAD,...
                        'SaveInd',false);
    
    toc
    
end

