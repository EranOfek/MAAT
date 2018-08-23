function DirS=list_fun_in_package(Package)
% Find all functions in a matlab package.
% Package: Util.files
% Description: Find all functions in a matlab package.
% Input  : - Package name
%            E.g., '+AstroUtil/+spec' or 'AstroUtil.spec'.
% Output : - A structure containing the directory path and file content.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: DirS=Util.files.list_fun_in_package('AstroUtil.spec')
% Reliable: 2
%--------------------------------------------------------------------------

if (strfind(Package,'+'))
    % Package name has a valid format - e.g., +AstroUtil/+spec
elseif (strfind(Package,'.'))
    % package name has a . format - e.g., AstroUtil.spec
    % convert to +...
    
    Package = strrep(Package,'.','/+');
    Package = sprintf('+%s',Package);
    
end
DirS = what(Package);


