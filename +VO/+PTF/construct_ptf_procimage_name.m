function FileName=construct_ptf_procimage_name(Data)
% SHORT DESCRIPTION HERE
% Package: VO.PTF
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: FileName=VO.PTF.construct_ptf_procimage_name(ProcImagesPTF.Cat)
% Reliable: 
%--------------------------------------------------------------------------
error('not working')

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if istable(Data)
    Year  = Data.Year;
    Month = Data.Month;
    Day   = Data.Day;
    Filter= Data.Filter;
    CCD   = Data.CCD;
    P     = Data.P;
    Ver   = Data.Ver;
    FN_Date = Data.FN_Date;
    FN_Time = Data.FN_Time;
    FN_Uni  = Data.FN_Uni;
    Field   = Data.Field;
elseif (AstCat.isastcat(Data))
    Year  = col_select(Data,'Year');
    Month = col_select(Data,'Month');
    Day   = col_select(Data,'Day');
    Filter= col_select(Data,'Filter');
    CCD   = col_select(Data,'CCD');
    P     = col_select(Data,'P');
    Ver   = col_select(Data,'Ver');
    FN_Date = col_select(Data,'FN_Date');
    FN_Time = col_select(Data,'FN_Time');
    FN_Uni  = col_select(Data,'FN_Uni');
    Field   = col_select(Data,'Field');
else
    Year    = Data(:,1);
    Month   = Data(:,2);
    Day     = Data(:,3);
    Filter  = Data(:,4);
    CCD     = Data(:,5);
    P       = Data(:,6);
    Ver     = Data(:,7);
    FN_Date = Data(:,8);
    FN_Time = Data(:,9);
    FN_Uni  = Data(:,10);
    Field   = Data(:,11);
end

    
%'/ptf/pos/archive/proc/2009/03/01/f1/c0/p5/v2/PTF_200903014118_i_p_scie_t095257_u012413683_f01_p000192_c00.fits'
N = numel(Year);
FileName = cell(N,1);
for I=1:1:N
    FileName{I} = sprintf('/ptf/pos/archive/proc/%04d/%02d/%02d/f%d/c%d/p%d/v%d/PTF_%d_i_p_scie_t%06d_u%09d_f%02d_p%06d_c%02d.fits',...
            Year(I), Month(I), Day(I), Filter(I), CCD(I), P(I), Ver(I),...
            FN_Date(I), FN_Time(I), FN_Uni(I), Filter(I), Field(I), CCD(I));
end