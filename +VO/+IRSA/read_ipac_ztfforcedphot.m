function [Table,FilterList]=read_ipac_ztfforcedphot(File,varargin)
% Read ZTF forced photometry file
% Package: VO.IRSA
% Description: Read ZTF forced photometry file generated using
%              VO.ZTF.wget_irsa_forcedphot_diff
% Input  : - File name
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - Table with data
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Table,FilterList]=VO.IRSA.read_ipac_ztfforcedphot('forcedphotometry_req00000265_lc.txt');
% Reliable: 
%--------------------------------------------------------------------------


DefV.BinSize               = 1;
DefV.BaseZP                = 25;
DefV.Max_forcediffimchisq  = 1.5;
DefV.Max_zpmaginpscirms    = 0.05;
DefV.Val_procstatus        = 0;
DefV.DelBadData            = true; 
DefV.SubMedFlux            = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

ColCellString = '# Order of columns below:';

FID = fopen(File,'r');
Ic = 0;
Iline = 0;
Itab  = 0;
while (~feof(FID))
    L = fgets(FID);
    ColStr = false;
    if (numel(L)>=numel(ColCellString))
        if strcmp(L(1:numel(ColCellString)),ColCellString)
            % read list of columns
            L = fgetl(FID);
            ColCell = regexp(L,',','split');
            ColCell = strtrim(ColCell); % remove blanks
            Col     = cell2struct(num2cell(1:1:numel(ColCell)),ColCell,2);
            ColStr  = true;
            Ncol    = numel(ColCell);
            
            Table   = array2table(zeros(0,Ncol));
            Table.Properties.VariableNames = ColCell;
            
        end
    end
            
    if (~ColStr)
        
        Iline = Iline + 1;
        if (strcmp(L(1),'#'))
            % comment
            Ic = Ic + 1;
            Comment{Ic} = L;
        else
            Itab = Itab + 1;
            C = regexp(L,'\s','split');
            C = C(2:end-1);
            
            Num = str2double(C);
            Table(Itab,:) = num2cell(Num);
            %Table(Itab,isnan(Num)) = C(isnan(Num));
            switch lower(C{Col.filter})
                case 'ztf_g'
                    Filter = 1;
                case 'ztf_r'
                    Filter = 2;
                case 'ztf_i'
                    Filter = 3;
                otherwise
                    error('Unknown filter option');
            end
            
            Table.filter(Itab) = Filter;
        end
    end
end
fclose(FID);

% add columns of normalized flux
FluxFactor = 10.^(-0.4.*Table.zpdiff) ./10.^(-0.4.*InPar.BaseZP);
Table      = addvars(Table,FluxFactor);
Flux       =  Table.forcediffimflux.*FluxFactor;
FluxErr    = Table.forcediffimfluxuncap.*FluxFactor;
Table      = addvars(Table,Flux,FluxErr);

% remove bad points

% Table.forcediffimchisq <1.5
% Table.procststus
% Table.zpmaginpscirms < 0.05

FlagGood = Table.forcediffimchisq < InPar.Max_forcediffimchisq & ...
           Table.zpmaginpscirms   < InPar.Max_zpmaginpscirms & ...
           Table.procstatus      == InPar.Val_procstatus;
      
Table      = addvars(Table,FlagGood);

%plot(Table.jd-2450000,Table.forcediffimflux,'.')  % forcediffimfluxuncap
%B = timeseries.binning([Table.jd, Table.forcediffimflux, Table.forcediffimfluxuncap],InPar.BinSize,[],{'MeanBin',@nanmean,@rstd,@numel});
%Data = [Table.jd, Table.forcediffimflux, Table.forcediffimfluxuncap];

Data = [Table.jd, Table.Flux, Table.FluxErr];
if (InPar.DelBadData)
    Data = Data(FlagGood,:);
    Table = Table(FlagGood,:);
end



UniqueFilter = unique(Table.filter);
Nfilter      = numel(UniqueFilter);
for Ifilter=1:1:Nfilter
    FilterList(Ifilter).FilterInd = UniqueFilter(Ifilter);
    FilterList(Ifilter).Ind       = find(Table.filter == UniqueFilter(Ifilter));
    
    MedianFlux = median(Table.Flux(FilterList(Ifilter).Ind));
    FilterList(Ifilter).MedianFlux = MedianFlux;
    FilterList(Ifilter).rstdFlux   = Util.stat.rstd(Table.Flux(FilterList(Ifilter).Ind));
    if (InPar.SubMedFlux)
        
        % subtract median flux
        Table.Flux(FilterList(Ifilter).Ind) = Table.Flux(FilterList(Ifilter).Ind) - MedianFlux;
    end
    

    B    = timeseries.binning(Data(FilterList(Ifilter).Ind,:),InPar.BinSize,[NaN NaN],{'MeanBin',@nanmedian,@Util.stat.rstd,@numel});
    Flag = ~isnan(B(:,1));
    B    = B(Flag,:);
    Err  = B(:,3)./sqrt(B(:,4));
    
    FilterList(Ifilter).B         = [B, Err];
    FilterList(Ifilter).B(:,1)    = FilterList(Ifilter).B(:,1);
end








