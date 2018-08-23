function [All,Col,JD0,ListCell,JumpInd]=read_corot_file(List,ExpTime)
% Read CoRoT file
% Package: VO.CoRoT
% Description: Read CoRoT fits files.
% Input  : - String containing a CoRoT fits file names in one of the
%            following formats:
%            (1) A cell vector containing a file name in each cell.
%            (2) A string containing wild cards (i.e., '*' or '?'),
%                in this case, the program will produce a list of
%                all matched files in the current directory.
%            (3) A string containing a file name (the file name not
%                necesserly exist).
%            (4) A string begining with '@' containing a file name
%                that contains a list of file (one per line).
%          - Exposure time in second. If given than find also
%            points in the LCs which are not seperated by
%            the exposure time+/-0.1s.
% Output : - Cell array containing the LCs. One cell per target.
%            Each cell containing 8 (for 'MON') or 14 (for 'CHR')
%            cells each contain one column.
%            fits file.
%          - Col{1} - structure containing fields in 'CHR' files.
%            Col{2} - structure containing fields in 'MON' files.
%          - JD0 to add to tabulated JDs.
%          - Cell array containing all the fits files names.
%          - Cell array of vectors. For each file indicate the indices
%            in the file in which a discontinuity was found.
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Jun 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%------------------------------------------------------------------------------
Def.ExpTime = [];
if (nargin==1)
   ExpTime = Def.ExpTime;
elseif (nargin==2)
   % do nothing
else
   error('Illegal number of input arguments');
end

[ListFileName,ListCell] = Util.files.create_list(List);


JD0 = 2451545.00025401;

Type.CHR = 1;
Type.MON = 2;
Col{Type.CHR}.Date          = 1;    % Calendar date 
Col{Type.CHR}.CorotSJD      = 2;    % date of the end of the measurements in the satel
                                    % for 32s exposure - refers to the end of the exposure.
                                    % for 512s exposure refers to 32s after begining of exposure.
Col{Type.CHR}.CorotHJD      = 3;    % date of the end of the measurements in the helio
                                    % for 32s exposure - refers to the end of the exposure.
                                    % for 512s exposure refers to 32s after begining of exposure.
Col{Type.CHR}.Status        = 4;    % flags indicating the status of measurements
Col{Type.CHR}.RedFlux       = 5;    % Red flux [electrons]
Col{Type.CHR}.StD_RedFlux   = 6;    % StD in red flux [electrons]
Col{Type.CHR}.GreenFlux     = 7;    % Green flux [electrons]
Col{Type.CHR}.StD_GreenFlux = 8;    % StD in green flux [electrons]
Col{Type.CHR}.BlueFlux      = 9;    % Blue flux [electrons]
Col{Type.CHR}.StD_BlueFlux  = 10;   % StD in blue flux [electrons]
Col{Type.CHR}.BG            = 11;   % Subtracted background level [electros/pixel]
Col{Type.CHR}.CorrRed       = 12;   % Difference between the corresponding N1 LC and
Col{Type.CHR}.CorrGreen     = 13;   % Difference between the corresponding N1 LC and
Col{Type.CHR}.CorrBlue      = 14;   % Difference between the corresponding N1 LC and


Col{Type.MON}.Date      = 1;    % Calendar date
Col{Type.MON}.CorotSJD  = 2;    % date of the end of the measurements in the satel
Col{Type.MON}.CorotHJD  = 3;    % date of the end of the measurements in the helio
Col{Type.MON}.Status    = 4;    % flags indicating the status of measurements
Col{Type.MON}.WhiteFlux = 5;    % white flux
Col{Type.MON}.StD_WF    = 6;    % Standard deviation associated with the calculat
Col{Type.MON}.Back      = 7;    % Subtracted background level
Col{Type.MON}.Corr      = 8;    % Difference between the corresponding N1 LC and


Nl = length(ListCell);
for Il=1:1:Nl
   All(Il).Table = fitsread(ListCell{Il},'BinTable');
   if (isempty(strfind('CHR',ListCell{Il}))==0)
      All(Il).Type  = Type.CHR;
   elseif (isempty(strfind('MON',ListCell{Il}))==0)
      All(Il).Type  = Type.MON;
   else
      error('Uknown file type');
   end

   if (isempty(ExpTime)==1)
      JumpInd{Il} = NaN;
   else
      JumpInd{Il} = find(abs(diff(All(Il).Table{2}.*86400)-ExpTime)>0.1);
   end
end

