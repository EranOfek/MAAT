function Cat=wget_pulsar_cat(varargin)
% Read ATNF pulsar catalog from the web into an AstCat object.
% Package: VO.prep
% Description: Read ATNF pulsar catalog from the web into an AstCat object.
% Input  : null 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - AstCat object containing the ATNF pulsar catalog.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Pulsar = VO.prep.wget_pulsar_cat;
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;
CatField     = AstCat.CatField;
ColCellField = AstCat.ColCellField;

%DefV.Ver                  = '1.58';
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

URL = 'http://www.atnf.csiro.au/people/pulsar/psrcat/proc_form.php?version=1.58&PMRA=PMRA&PMDec=PMDec&PX=PX&PosEpoch=PosEpoch&RaJD=RaJD&DecJD=DecJD&P0=P0&P1=P1&F0=F0&F1=F1&F2=F2&PEpoch=PEpoch&DM=DM&DM1=DM1&RM=RM&W50=W50&W10=W10&Tau_sc=Tau_sc&S400=S400&S1400=S1400&S2000=S2000&Dist=Dist&Dist_DM=Dist_DM&R_lum=R_lum&R_lum14=R_lum14&Age=Age&Bsurf=Bsurf&Edot=Edot&Edotd2=Edotd2&PMtot=PMtot&VTrans=VTrans&P1_i=P1_i&Age_i=Age_i&Bsurf_i=Bsurf_i&B_LC=B_LC&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=decjd&sort_order=asc&condition=&pulsar_names=&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+csv+with+errors&no_value=*&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=55&table_bottom.y=25';
Str = webread(URL);
TS  = regexp(Str,'<pre>(?<t>.+)</pre>','names');               
TableStr = TS.t;

% #;PMRA;;;PMDEC;;;PX;;;POSEPOCH;;RAJD;DECJD;P0;;;P1;;;F0;;;F1;;;F2;;;PEPOCH;;DM;;;DM1;;;RM;;;W50;;;W10;;;TAU_SC;;;S400;;;S1400;;;S2000;;;DIST;DIST_DM;;R_LUM;R_LUM14;AGE;BSURF;EDOT;EDOTD2;PMTOT;;VTRANS;P1_I;AGE_I;BSURF_I;B_LC;
ColCell = {'Junk'  ,'PMRA','ErrPMRA' ,'Junk','PMDec' ,'ErrPMDec','Junk','Plx','ErrPlx','Junk','PosEpoch','Junk','RA','Dec','P0','Junk','Junk','P1','Junk','Junk','F0','Junk','Junk','F1','Junk','Junk','F2','Junk','Junk','PeriodEpoch','Junk','DM','Junk','Junk','DM1','Junk','Junk','RM','Junk','Junk','W50','Junk','Junk','W10','Junk','Junk','TAU_SC','Junk','Junk','S400','Junk','Junk','S1400','Junk','Junk','S2000','Junk','Junk','Dist','DistDM','Junk','R_Lum'    ,'R_Lum14'  ,'Age','Bsurf','Edot' ,'EdotD2'         ,'PM_tot','Junk','V_trans','P1_I','Age_I','Bsurf_I','B_LC'};
ColUnits= {''      ,'mas/yr','mas/yr',''    ,'mas/yr','mas/yr'  ,''    ,'mas','mas'   ,''    ,'MJD'     ,''    ,'rad','rad','s',''    ,''    ,''  ,''    ,''    ,'Hz',''    ,''    ,'s^-2',''  ,''    ,'s^-3',''  ,''    ,'MJD'        ,''    ,'pc*cm^-3','',''  ,'pc*cm^-3/yr','','','rad*m^-2','',''  ,'ms' ,''    ,''    ,'ms' ,''    ,''    ,'s'     ,''    ,''    ,'mJy' ,''    ,''    ,'mJy'  ,''    ,''    ,'mJy'  ,''    ,''    ,'kpc' ,'kpc'   ,''    ,'mJy*kpc^2','mJy*kpc^2','yr' ,'G'    ,'erg/s','erg*kpc^-2*s^-2','mas/yr',''    ,'km/s'   ,''    ,'yr'   ,'G'      ,'G'};

% ;(mas/yr);;;(mas/yr);;;(mas);;;(MJD);;(deg);(deg);(s);;;;;;(Hz);;;(s^-2);;;(s^-3);;;(MJD);;(cm^-3pc);;;(cm^-3pc/yr);;;(radm^-2);;;(ms);;;(ms);;;(s);;;(mJy);;;(mJy);;;(mJy);;;(kpc);(kpc);;(mJykpc^2);(mJykpc^2);(Yr);(G);(ergs/s);(ergs/kpc^2/s);(mas/yr);;(km/s);;(Yr);(G);(G);

% prep format string
Ncol = numel(ColCell);
Format = '';
for Icol=1:1:Ncol
    switch lower(ColCell{Icol})
        case 'junk'
            Fs = '%*s';
        otherwise
            Fs = '%f';
    end
    Format = sprintf('%s%s ',Format,Fs);
end
Format = sprintf('%s\\n',Format);
Iok = find(~strcmp(lower(ColCell),'junk'));
ColCell  = ColCell(Iok);
ColUnits = ColUnits(Iok);

C = textscan(TableStr,Format,'Delimiter',';','CommentStyle','%','HeaderLines',3,'TreatAsEmpty','*');
Cat = AstCat;
Cat.(CatField) = cell2mat(C);
NewColInd = [8,9,1:7,10:numel(ColCell)];
Cat.(CatField) = Cat.(CatField)(:,NewColInd);
Cat.(ColCellField) = ColCell(NewColInd);
Cat = colcell2col(Cat);
Cat.(CatField)(:,1:2) = Cat.(CatField)(:,1:2)./RAD;
Cat = sortrows(Cat,2);
Cat.ColUnits = ColUnits(NewColInd);
