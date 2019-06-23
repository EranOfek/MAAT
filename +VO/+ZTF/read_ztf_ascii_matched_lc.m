function [Data,ColCell]=read_ztf_ascii_matched_lc(File,varargin)
% Read ZTF ascii file of matched light curves
% Package: VO.ZTF
% Description: Read ZTF ascii file of matched light curves.
%              Optionally save the data and meta data in an HDF5 file.
%              The matched light curves files and description are
%              available from:
%              https://www.ztf.caltech.edu/page/dr1#12c
% Input  : - A string of the txt file name or a numeric scalar of
%            the ZTF field ID.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FileNamePrefix' - Default is 'field'.
%            'FileNameSuffix' - Default is '.txt'.
% Output : - Structure array of all sources with their light curves.
%            Available fields are:
%            'ID' - Source ID
%            'Nep' - Number of epochs
%            'FilterID' - Filter ID
%            'Field' - ZTF field index.
%            'RcID' - Readout channel ID (CCD/quad) 0 to 63
%            'RA' - J2000.0 R.A. [deg] of source in reference image.
%            'Dec' - J2000.0 Dec. [deg] of source in reference image.
%            'LC' - Array of source lightcurve [columns X times]
%          - Cell array of columns in the LC array.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Data,ColCell]=VO.ZTF.read_ztf_ascii_matched_lc(868);
%          % example: convert all txt files tp HDF5 including statistical properties
%          F=dir('field*.txt');
%          for I=56:1:numel(F),
%               I
%              FI=str2double(F(I).name(6:11));
%              FN=sprintf('ztfLCDR1_%06d.hdf5',FI)
%              [~,ColCell]=VO.ZTF.read_ztf_ascii_matched_lc(FI,'H5_FileName',FN);
%          end
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.FileNamePrefix       = 'field';
DefV.FileNameSuffix       = '.txt';
DefV.H5_FileName          = '';
DefV.H5_DataSetInd        = '/IndAllLC';  
DefV.H5_DataSetLC         = '/AllLC'; 
DefV.Nwrite               = 100000;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (ischar(File))
    FileName = File;
elseif (isnumeric(File))
    F = dir(sprintf('%s%06d*%s',InPar.FileNamePrefix,File,InPar.FileNameSuffix));
    if (isempty(F))
        error('File name constructed from field ID was not found');
    end
    if (numel(F)>1)
        warnning('Multiple files were found for field ID, reading first');
    end
    FileName = F(1).name;
else
    error('First argument File name should be string or numeric');
end
    
FID = fopen(FileName,'r');

ColCell = {'HMJD','Mag','MagErr','ColorCoef','Flags'};
Col     = cell2struct(num2cell(1:1:numel(ColCell)),ColCell,2);

ColCellInd = {'RA','Dec','I1','I2','Nep', 'ID',...
                  'FilterID','Field','RcID',...
                  'MeanMag', 'StdMag', 'RStdMag', 'MaxMag', 'MinMag',...
                  'Chi2', 'MaxPower', 'FreqMaxPower'};
NcolInd = numel(ColCellInd);

              
Ncol    = numel(ColCell);
%catflags: Photometric/image quality flags encoded as bits (Section 9b).
%            In particular, you will always want to exclude observation epochs
%            affected by clouds and/or the moon. These epochs have catflags = 32768
%            (decimal bit 15).


             
Iobj = 0;
IndAllLC = zeros(0,NcolInd);

Nlast = 0;
IndSave = 0;
IobjLastSave = 1;


while ~feof(FID)
    Iobj = Iobj + 1;
    Header = fscanf(FID,'%s %ld %d %d %d %d %f %f\n',8);
    Nep = Header(3);
    LC = fscanf(FID,'%f %f %f %f %d\n',5.*Nep);
    
    if (numel(Header)~=8)
        Header(1)
        Iobj
        error('Nh~=8');
    end
    
    Data(Iobj).ID       = Header(2);
    Data(Iobj).Nep      = Nep;
    Data(Iobj).FilterID = Header(4);
    Data(Iobj).Field    = Header(5);
    Data(Iobj).RcID     = Header(6);  % readout channel ID (CCD/quad 0..63)
    Data(Iobj).RA       = Header(7);  % RA in ref
    Data(Iobj).Dec      = Header(8);
    
    if (isempty(LC))
        Data(Iobj).LC = [];
        
        FIDp = fopen('Problems.txt','a+');
        fprintf(FIDp,'ID %d\n',Header(2));
        fclose(FIDp);
    else
        try
            Data(Iobj).LC       = reshape(LC,Ncol,Nep);
        catch
            Data(Iobj).LC = [];
        
            FIDp = fopen('Problems.txt','a+');
            fprintf(FIDp,'ID %d\n',Header(2));
            fclose(FIDp);

        
            size(LC)
            Ncol
            Nep
            Iobj
        end
    end
    if ~isempty(InPar.H5_FileName)
        if (Iobj==1)
            h5create(InPar.H5_FileName,InPar.H5_DataSetLC,[Inf 5],'ChunkSize',[5 5]);
            HDF5.writeatt(InPar.H5_FileName,InPar.H5_DataSetLC,[ColCell',num2cell(1:1:numel(ColCell))']');
            
            Istart = 1;
            
            Pool = parpool(12);

        end
        
        if (feof(FID) || Iobj./InPar.Nwrite==floor(Iobj./InPar.Nwrite))
            % write to file
            
            IndSave = IndSave + 1;
            
            LC   = [Data(IobjLastSave:Iobj).LC]';
            
            Nrow = size(LC,1);
            Iend = Istart + Nrow - 1;
            
            h5write(InPar.H5_FileName,InPar.H5_DataSetLC,LC,[Istart 1],[Nrow Ncol]);
            Istart = Istart + Nrow;
            
            % calc properties
            Prop = calc_prop(Data(IobjLastSave:Iobj),Nlast,Col);
            Nlast = Nlast + size(LC,1);
            IndAllLC = [IndAllLC; Prop];
            
            % clear LC
            for II=IobjLastSave:Iobj
                Data(II).LC = [];
            end
            
            
            IobjLastSave = Iobj+1;
        end
        
    end
    
end

fclose(FID);
if ~isempty(InPar.H5_FileName)
    

    KeysInd = [ColCellInd', num2cell(1:1:numel(ColCellInd))'];
    h5create(InPar.H5_FileName,InPar.H5_DataSetInd,size(IndAllLC));
    h5write(InPar.H5_FileName,InPar.H5_DataSetInd,IndAllLC);
    HDF5.writeatt(InPar.H5_FileName,InPar.H5_DataSetInd,KeysInd');
    
    
    delete(gcp('nocreate'))
     
end

end % main function

%%
function IndAllLC=calc_prop(Data,Nlast,Col)
    %
    
    RAD = 180./pi;

    Nobj = numel(Data);
    Nep = [Data.Nep]';

    I1  = cumsum([1;Nep(1:end-1)]);
    I2  = [I1(2:end)-1; I1(end)+Nep(end)-1];
    I1  = I1 + Nlast;
    I2  = I2 + Nlast;
    
    
    MeanMag = nan(Nobj,1);
    StdMag  = nan(Nobj,1);
    RStdMag = nan(Nobj,1);
    MaxMag  = nan(Nobj,1);
    MinMag  = nan(Nobj,1);
    Chi2    = nan(Nobj,1);
    MaxPower= nan(Nobj,1);
    FreqMaxPower= nan(Nobj,1);
    
    FreqVec = (0:1./600:4)';
    
    %tic;
    parfor IobjP=1:1:Nobj
        if (~isempty(Data(IobjP).LC))
            MeanMag(IobjP) = mean(Data(IobjP).LC(Col.Mag,:));

            StdMag(IobjP)  = std(Data(IobjP).LC(Col.Mag,:));

            RStdMag(IobjP) = Util.stat.rstd(Data(IobjP).LC(Col.Mag,:).');

            MaxMag(IobjP) = max(Data(IobjP).LC(Col.Mag,:)) - MeanMag(IobjP);
            MinMag(IobjP) = MeanMag(IobjP) - min(Data(IobjP).LC(Col.Mag,:));
            Chi2(IobjP)    = sum((Data(IobjP).LC(Col.Mag,:) - MeanMag(IobjP)).^2./Data(IobjP).LC(Col.MagErr,:).^2);

            P = timeseries.period_normnl(Data(IobjP).LC( [Col.HMJD, Col.Mag],:).',FreqVec);

            [MaxPower(IobjP),MaxInd] = max(P(:,2));
            FreqMaxPower(IobjP) = FreqVec(MaxInd);
        end
    end
    %toc
    
    
    
    %%
    % RA Dec I1, I2, Nep, ID, FilterID, Field, RcID,
    % MeanMag, StdMag, RStdMag, MaxMag, MinMag, Chi2, MaxPower,
    % FreqMaxPower

    IndAllLC = [ [Data.RA].'./RAD, [Data.Dec].'./RAD, I1, I2, [Data.Nep].',...
                 [Data.ID].', [Data.FilterID].',...
                 [Data.Field].', [Data.RcID].',...
                 MeanMag, StdMag, RStdMag, MaxMag,MinMag, Chi2, MaxPower, FreqMaxPower];
end