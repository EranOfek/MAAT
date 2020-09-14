function [Out]=wget_goes_data(varargin)
% wget GOES partcles flux data
% Package: ultrasat
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Out]=ultrasat.wget_goes_data;
% Reliable: 
%--------------------------------------------------------------------------

DefV.Year                 = 2014;
DefV.Month                = (1:1:12);
DefV.Sat                  = 'goes13';
DefV.FType                = {'epead_e1ew','epead_e2ew','epead_e3ew','maged_19me1','maged_19me2','maged_19me3','maged_19me4','maged_19me5'};
DefV.BaseURL              = 'https://satdat.ngdc.noaa.gov/sem/goes/data/full/%04d/%02d/%s/csv/';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


%https://satdat.ngdc.noaa.gov/sem/goes/data/full/2015/01/goes15/csv/

%FType = {'epead_e1ew','epead_e2ew','epead_e3ew','maged_19me1','maged_19me2','maged_19me3','maged_19me4','maged_19me5'};
%FType  = {'epead_e2ew'};


%Energy= [0.8 2  ?     ? ? ? (100-200 / 150) (200-350) ?
%        >0.8   >2
% e/(cm^2 s sr)



VecMonth = InPar.Month;
Nm = numel(VecMonth);

FType = InPar.FType;    
Nftype = numel(FType);

for Iftype=1:1:Nftype
   

    for Im=1:1:Nm
        Month = VecMonth(Im);

        DirURL = sprintf(InPar.BaseURL,InPar.Year,Month,InPar.Sat);


        List = www.find_urls(DirURL,'strfind','.csv');

        SpL = regexp(List,'_','split');
        JD = zeros(numel(SpL),1);
        for Ispl=1:1:numel(SpL)
            %Date = datevec(SpL{Ispl}{end}(1:end-4),'yyyymmdd')
            Date = datevec(SpL{Ispl}{end-1},'yyyymmdd');
            JD(Ispl) = celestial.time.julday(Date(:,[3 2 1]));
        end

        
        Ilist= find(~Util.cell.isempty_cell(strfind(List,FType{Iftype})));
        Nl = numel(Ilist);

        for Il=1:1:Nl
            [Iftype, Im, Il]
            
            FileTmp = tempname;
            FileTmp = websave(FileTmp,List{Ilist(Il)});

            [Mat1,ColCell]=VO.Util.read_csv_with_header(FileTmp,'OutType','table','Delimiter',',','StartKey','data:');

            delete(FileTmp)
            Mat1.time_tag = celestial.time.julday(Mat1.time_tag);

            if (Il==1 && Im==1)
                Mat = Mat1;
            else
                Mat = [Mat; Mat1];
            end

        end
        Out(Iftype).FType = FType{Iftype};
        Out(Iftype).Mat = Mat;

    end


end