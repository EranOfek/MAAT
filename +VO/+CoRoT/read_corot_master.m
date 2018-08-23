function [Master,Ind,FullFileName,StartJD,EndJD,TbCol]=read_corot_master
% Read CoRoT master file
% Package: VO.CoRoT
% Description: Read the CoRoT master file - this file contains a list of
%              all the CoRoT target stars with their mean properties.
% Input  : null.
% Output : - Cell array containing the master file - cell per column.
%          - Structure containing flags indicating by which run/ccd a star
%            was observed.
%          - Cell array of full file names.
%          - Start JD column.
%          - End JD column.
%          - Structure of selected table columns.
% Instellation: Requires the CoRoT master file
% Tested : Matlan 7.8
%     By : Eran O. Ofek                    Oct 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%------------------------------------------------------------------------------

TableFile = 'all_data.dat';
Dir = Util.files.which_dir('read_corot_master');
TableDir = 'LC';
DirLoc = sprintf('%s%s%s%s',Dir,filesep,TableDir,filesep);
MasterFile = sprintf('%s%s',DirLoc,TableFile);

% Generated from: http://nsted.ipac.caltech.edu

TbCol.RunCode   = 2;
TbCol.CCDhalf   = 3;
TbCol.StartDate = 6;
TbCol.EndDate   = 7;
TbCol.ExpTime   = 20;
TbCol.FlagMON   = 39;
TbCol.FlagCHR   = 40;
TbCol.FileName  = 44;


I =   1; Col(I).Name = 'corotid';       Col(I).Format = '%s';
I =   2; Col(I).Name = 'run_code';      Col(I).Format = '%s';
I =   3; Col(I).Name = 'hlfccdid';      Col(I).Format = '%s';
I =   4; Col(I).Name = 'alpha';         Col(I).Format = '%f';
I =   5; Col(I).Name = 'delta';         Col(I).Format = '%f';
I =   6; Col(I).Name = 'start_date';    Col(I).Format = '%s';
I =   7; Col(I).Name = 'end_date';      Col(I).Format = '%s';
I =   8; Col(I).Name = 'mag_b';         Col(I).Format = '%f';
I =   9; Col(I).Name = 'mag_v';         Col(I).Format = '%f';
I =  10; Col(I).Name = 'mag_r';         Col(I).Format = '%f';
I =  11; Col(I).Name = 'mag_i';         Col(I).Format = '%f';
I =  12; Col(I).Name = 'lc_mean_b';     Col(I).Format = '%f';
I =  13; Col(I).Name = 'lc_rms_b';      Col(I).Format = '%f';
I =  14; Col(I).Name = 'lc_mean_g';     Col(I).Format = '%f';
I =  15; Col(I).Name = 'lc_rms_g';      Col(I).Format = '%f';
I =  16; Col(I).Name = 'lc_mean_r';     Col(I).Format = '%f';
I =  17; Col(I).Name = 'lc_rms_r';      Col(I).Format = '%f';
I =  18; Col(I).Name = 'lc_mean_white'; Col(I).Format = '%f';
I =  19; Col(I).Name = 'lc_rms_white';  Col(I).Format = '%f';
I =  20; Col(I).Name = 'exptime';       Col(I).Format = '%f';
I =  21; Col(I).Name = 'degree_chrom';  Col(I).Format = '%f';
I =  22; Col(I).Name = 'star_temper';   Col(I).Format = '%f';
I =  23; Col(I).Name = 'contam_factor'; Col(I).Format = '%f';
I =  24; Col(I).Name = 'active_level';  Col(I).Format = '%f';
I =  25; Col(I).Name = 'lum_class';     Col(I).Format = '%s';
I =  26; Col(I).Name = 'var_class1';    Col(I).Format = '%s';
I =  27; Col(I).Name = 'prob_var1';     Col(I).Format = '%f';
I =  28; Col(I).Name = 'var_class2';    Col(I).Format = '%s';
I =  29; Col(I).Name = 'prob_var2';     Col(I).Format = '%f';
I =  30; Col(I).Name = 'var_class3';    Col(I).Format = '%s';
I =  31; Col(I).Name = 'prob_var3';     Col(I).Format = '%f';
I =  32; Col(I).Name = 'spec_type';     Col(I).Format = '%s';
I =  33; Col(I).Name = 'B-V';           Col(I).Format = '%f';
I =  34; Col(I).Name = 'num_hot_pix';   Col(I).Format = '%f';
I =  35; Col(I).Name = 'X_ccd';         Col(I).Format = '%f';
I =  36; Col(I).Name = 'Y_ccd';         Col(I).Format = '%f';
I =  37; Col(I).Name = 'X_win_ccd';     Col(I).Format = '%f';
I =  38; Col(I).Name = 'Y_win_ccd';     Col(I).Format = '%f';
I =  39; Col(I).Name = 'flag_MON';      Col(I).Format = '%f';   % flag indicating presence of exo monochromatic lc data product
I =  40; Col(I).Name = 'flag_CHR';      Col(I).Format = '%f';   % flag indicating presence of exo chromatic lc data product
I =  41; Col(I).Name = 'ewindescriptor';Col(I).Format = '%f';   % flag indicating presence of exo windescriptor data product
I =  42; Col(I).Name = 'an2_star';      Col(I).Format = '%f';   % flag indicating presence of astero lc data product
I =  43; Col(I).Name = 'awindescriptor';Col(I).Format = '%f';   % flag indicating presence of astero windescriptor data product
I =  44; Col(I).Name = 'file_name';     Col(I).Format = '%s';
I =  45; Col(I).Name = 'version_dir';   Col(I).Format = '%s';

N = length(Col);
Format = '';
for I=1:1:N
   Format = sprintf('%s %s',Format,Col(I).Format);
end
Format = sprintf('%s\\n',Format);

FID = fopen(MasterFile,'r');
Master = textscan(FID,Format,'CommentStyle','%','Delimiter','|');
fclose(FID);

N = length(Master{1});

%--- Convert start/end dates to JD ---
Date    = datevec(Master{TbCol.StartDate},'yyyy-mm-dd');
StartJD = celestial.time.julday(Date(:,[3 2 1 4 5 6]));
Date    = datevec(Master{TbCol.EndDate},'yyyy-mm-dd');
EndJD   = celestial.time.julday(Date(:,[3 2 1 4 5 6]));


%--- Seperate CCDs sections ---
Master{TbCol.CCDhalf} = deblank(Master{TbCol.CCDhalf});
[CCDs,CCDind,CCDflag] = Util.cell.group_cellstr(Master{TbCol.CCDhalf});
Nccd = length(CCDs);
Ind.CCD.junk = [];
for Iccd=1:1:Nccd
   Ind.CCD = setfield(Ind.CCD,CCDs{Iccd},CCDflag{Iccd});
end
Ind.CCD = rmfield(Ind.CCD,'junk');

%--- Seperate Runs ---
Master{TbCol.RunCode} = deblank(Master{TbCol.RunCode});
[RunCodes,RunCodesInd,RunCodesFlag]=Util.cell.group_cellstr(Master{TbCol.RunCode});
Nrun = length(RunCodes);
Ind.Run.junk = [];
for Irun=1:1:Nrun
   Ind.Run = setfield(Ind.Run,RunCodes{Irun},RunCodesFlag{Irun});
end
Ind.Run = rmfield(Ind.Run,'junk');


%--- Construct full file names ---
Str = cell(N,1);
[Str{find(Master{TbCol.FlagMON}==1)}] = deal('MON');
[Str{find(Master{TbCol.FlagCHR}==1)}] = deal('CHR');

Master{TbCol.FileName} = deblank(Master{TbCol.FileName});
FullFileName = cell(N,1);
for Is=1:1:N,
   FullFileName{Is} = sprintf('EN2_STAR_%s_%s.fits',Str{Is},Master{TbCol.FileName}{Is});
end
