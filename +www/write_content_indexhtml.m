function write_content_indexhtml(varargin)
% SHORT DESCRIPTION HERE
% Package: www
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Files=dir('*');
Nf = numel(Files);
FID = fopen('index.html','w');
for If=1:1:Nf
    
    if ~strcmp(Files(If).name(1),'.')
        if(Files(If).isdir)
            fprintf(FID,'%s      <a href="./%s/index.html">%s</a><br>\n',Files(If).date,Files(If).name, Files(If).name);
        else
            fprintf(FID,'%s      <a href="./%s">%s</a><br>\n',Files(If).date,Files(If).name, Files(If).name);
        end
    end
end
fclose(FID);

