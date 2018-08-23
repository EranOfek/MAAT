function Path=which_dir(Prog)
% Return the directory in which a matlab program resides.
% Package: Util.files
% Description: Return the directory in which a matlab program resides.
%              This program is a version of "which.m" that trim the
%              program name from the full path.
% Input  : - Matlab program name.
% Output : - Directory in which program resides.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       Feb 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: which_dir('fits_xy2rd')
% Reliable: 1
%------------------------------------------------------------------------------

FullPath = which(Prog);
[Path,Name] = fileparts(FullPath); 

