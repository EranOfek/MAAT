function Info=filename_info(ImList)
% Get information from PTF file name
% Package: VO.PTF
% Description: Get information from the IPAC PTF file name construct.
% Input  : - List of IPAC PTF FITS images (see create_list.m for options).
%            Default is '*.fits'.
% Output : - Info structure containing the extracted information for
%            each image. The following fields are available:
%            .Day
%            .Month
%            .Year
%            .FracDay
%            .JD
%            .ProductType
%            .ProcLevel
%            .ImType
%            .FilterID
%            .FilterName
%            .FieldID
%            .CCDID
%            .ID   - processed image ID
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    Sep 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Info=VO.PTF.filename_info('PTF_201209044888_c_p_scie_t114355_u013981212_f12_p002197_c04.ctlg');
% Reliable: 2
%-------------------------------------------------------------------------
import Util.files.*


Def.ImList = '*.fits';
if (nargin==0)
    ImList = Def.ImList;
elseif (nargin==1)
    % do nothing
else
    error('Illegal number of input arguments');
end

Map.PTF  = 1;
Map.Date = 2;
Map.ProductType = 3;   % 'i' | 'c'
Map.ProcLevel   = 4;   % 'r' | 'p' | 's' | 'e'
Map.ImType      = 5;   % 'scie' | 'mask' | ...
Map.Hour        = 6;
Map.ID          = 7;
Map.Filt        = 8;
Map.FieldID     = 9;
Map.CCDID       = 10;

FilterID  = [1 2 11 12];
FilterMap = {'g', 'R','Ha656','Ha663'};

[~,ListCell] = create_list(ImList,NaN);
Nim = length(ListCell);
Split = regexp(ListCell,'_','split');

for Iim=1:1:Nim
   Year   = str2double(Split{Iim}{Map.Date}(1:4));
   Month  = str2double(Split{Iim}{Map.Date}(5:6));
   Day    = str2double(Split{Iim}{Map.Date}(7:8));
   %FracDay= str2double(Split{Iim}{Map.Date}(9:12));
   Hour   = str2double(Split{Iim}{Map.Hour}(2:3));
   Min    = str2double(Split{Iim}{Map.Hour}(4:5));
   Sec    = str2double(Split{Iim}{Map.Hour}(6:7));
   FracDay= celestial.coo.convertdms([Hour Min Sec],'H','f');   
   JD     = celestial.time.julday([Day Month Year FracDay]);
   Info(Iim).Day     = Day;
   Info(Iim).Month   = Month;
   Info(Iim).Year    = Year;
   Info(Iim).FracDay = FracDay;
   Info(Iim).JD      = JD;
   Info(Iim).ProductType = Split{Iim}{Map.ProductType};
   Info(Iim).ProcLevel   = Split{Iim}{Map.ProcLevel};
   Info(Iim).ImType      = Split{Iim}{Map.ImType};
   Info(Iim).FilterID    = str2double(Split{Iim}{Map.Filt}(2:3));
   FilterIDID = find(Info(Iim).FilterID==FilterID);
   Info(Iim).FilterName  = FilterMap{FilterIDID};
   Info(Iim).FieldID     = str2double(Split{Iim}{Map.FieldID}(2:7));
   Info(Iim).CCDID       = str2double(Split{Iim}{Map.CCDID}(2:3));
   Info(Iim).ID          = str2double(Split{Iim}{Map.ID}(2:10));
end
   
   
   
