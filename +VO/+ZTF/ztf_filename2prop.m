function [Prop]=ztf_filename2prop(FileName,varargin)
% Extract information from ZTF file name
% Package: +VO/+ZTF
% Description: Extract information and image properties from ZTF file name.
% Input  : - ZTF File name, or a cell array of ZTF file name.
% Output : - A structure array with the following fields:
%            'FileType'
%            'Field'
%            'Filter'
%            'CCDID'
%            'QuadID'
%            'Year'
%            'Month'
%            'Day'
%            'Frac'
%            'JD'
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Prop]=VO.ZTF.ztf_filename2prop('ztf_20180815258044_000686_zr_c01_o_q1_psfcat.fits');
% Reliable: 2
%--------------------------------------------------------------------------


DefV.FilterList           = {'zg','zr','zi'};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if ~iscell(FileName)
    FileName = {FileName};
end
Nfile = numel(FileName);


Prop = Util.struct.struct_def({'FileType','Field','Filter','FilterInd','CCDID','QuadID','Year','Month','Day','Frac','JD'},Nfile,1);
for Ifile=1:1:Nfile
    Tmp = regexp(FileName{Ifile},'_','split');
    switch lower(Tmp{end})
        case 'psfcat.fits'
            Prop(Ifile).FileType = 'psfcat';
            Prop(Ifile).Field    = str2double(Tmp{3});
            Prop(Ifile).Filter   = Tmp{4};
            Prop(Ifile).FilterInd= find(strcmp(Prop(Ifile).Filter ,InPar.FilterList));
            Prop(Ifile).CCDID    = str2double(Tmp{5}(2:3));
            Prop(Ifile).QuadID   = str2double(Tmp{7}(2));
            Date    = datevec(Tmp{2},'yyyymmdd');
            
            FracStr = Tmp{2}(9:end);
            Ndigit  = numel(FracStr);
            Frac    = str2double(FracStr)./(10.^Ndigit);
            Prop(Ifile).Year    = Date(1);
            Prop(Ifile).Month   = Date(2);
            Prop(Ifile).Day     = Date(3);
            Prop(Ifile).Frac    = Frac;
            Prop(Ifile).JD      = convert.date2jd([Prop(Ifile).Day, Prop(Ifile).Month, Prop(Ifile).Year, Frac]);
        otherwise
            error('Unsupported ztf file type');
    end
end
            
            