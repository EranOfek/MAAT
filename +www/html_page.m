function html_page(Output,Content,varargin)
% Create an HTML file
% Package: www
% Description: Create an HTML file. The file contains the necessery
%              header and footer and a supplied content.
% Input  : - Output file name.
%          - A string containing file name, or a cell array of strings.
%            The file/strings will be added between the HTML header
%            and footer. If a cell vector is given then each element
%            is regarded as new line (i.e., seperated by '<br>').
%          * Arbitrary pairs of arguments:...keyword,value,...
%            to enable the user to control the web page appereance,
%            where keyword can be one of the following:
%            'PageTitle'  - Page title string.
%            'BgColor'    - HTML background color, defaukt is '#ffffff'.
%            'TextColor'  - HTML text color, defaukt is '#000000'.
%            'LinkColor'  - HTML link color, default is '#0000ff'.
%            'VLinkColor' - HTML vlink color, default is '#ff0000'.
%            'BodyStatment'-HTML (extra) body statment, default is ''.
% Output : null
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Apr 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: www.html_page('Example.html',{'Hello world'},'PageTitle','example');
% Reliable: 2
%--------------------------------------------------------------------------

Narg = length(varargin);
% set default values

DefV.PageTitle   = '';
DefV.BgColor     = '#ffffff';
DefV.TextColor   = '#000000';
DefV.LinkColor   = '#0000ff';
DefV.VLinkColor  = '#ff0000';
DefV.BodyStatment= '';

%InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);



FID = fopen(Output,'w');

%---------------------------
%--- Write HTML openinig ---
%---------------------------
fprintf(FID,'<HTML> \n');
fprintf(FID,'<HEAD> \n');
fprintf(FID,'<TITLE> \n');
fprintf(FID,'%s \n',InPar.PageTitle);
fprintf(FID,'</TITLE> \n');
fprintf(FID,'</HEAD> \n');
fprintf(FID,'<BODY bgcolor=%s text=%s link%s vlink=%s %s> \n',...
        InPar.BgColor,InPar.TextColor,InPar.LinkColor,InPar.VLinkColor,InPar.BodyStatment);

%fprintf(FID,'<TABLE border=%d %s> \n',TableBorder,TableStatment);
%
%fprintf(FID,'%s \n',PreTable);
%
%if (iscell(TableHead)==1),
%   fprintf(FID,'<tr> \n');
%   M= length(TableHead);
%   %--- write table header ---
%   for J=1:1:M,
%      fprintf(FID,'   <th %s> \n',HeadStatment{J});
%      fprintf(FID,'      %s \n',TableHead{J});   
%   end
%end

%--------------------------
%--- Write html content ---
%--------------------------
if (ischar(Content)==1)
   % file name
   StrContent = Util.files.file2str(Content,'str');
   fprintf(FID,'%s \n<br>\n',StrContent);
elseif (iscell(Content)==1)
   StrContent = Util.IO.fprintf_cell(FID,'%s \n<br>\n',Content);
else
   error('Content is of unknown type');
end



%--------------------------
%--- Write html closing ---
%--------------------------
fprintf(FID,'</BODY> \n');
fprintf(FID,'</HTML> \n');


fclose(FID);
